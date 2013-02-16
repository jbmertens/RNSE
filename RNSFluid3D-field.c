/*
To-do's:
- take in file name as argument
- Some basic checks (T_SAMPLEINT > 0, X_SAMPLEINT > 0, etc)
*/
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> /* mkdir() */
#include <math.h>
#include <zlib.h>
#include <hdf5.h>

#include "defines.h"

// Precision/format we'd like to use for this simulation:
typedef double simType;

/* data structure for storing calculated quantities at a point */
typedef struct {
  simType fields[6];
  simType gradients[4][DOF];
  simType derivs2[4];
  simType ut;
  simType ut2;
  simType u2;
  simType uudu;
  simType ji[4];
  simType Ji[4];
} PointData;

/* functions to calculate commonly used quantities */
simType Ut(PointData *paq);
simType Ut2(PointData *paq);
simType magu2(PointData *paq);

/* sums, derivatives, such. */
simType sumvv(simType v1[4], simType v2[4]);
simType sumvt(simType v1[4], simType t1[4][DOF], int rc, int s);
simType sumvtv(simType v1[4], simType t1[4][DOF], simType v2[4]);
simType sp_tr(simType t1[4][DOF]);
simType derivative(simType *data, int ddim, int dim, int i, int j, int k);
simType derivative2(simType *data, int ddim, int dim, int i, int j, int k);

/* potential function */
simType dV(simType phi);

/* write to file stuff */
void writeinfo(char *dir, char *root);
void dumpstate(simType *fields, int fwrites, int datasize, char *dir, char *root);
void readstate(simType *fields, int fwrites);

/* evolution functions for different DOF's */
void jsource(PointData *paq);
void Jsource(PointData *paq);
simType energy_evfn(PointData *paq);
simType fluid_evfn(PointData *paq, int u);
simType field_evfn(PointData *paq);
simType ddtfield_evfn(PointData *paq);
void calculatequantities(simType *fields, PointData *paq, int i, int j, int k);
void evolve(simType *initial, simType *final, simType coeff, PointData *paq, int i, int j, int k);


int main(int argc, char *argv[])
{
  /* iterators here first. */
  int i, j, k, n, u, s;

  /* information about files */
  if(argc != 3)
  {
    fprintf(stderr, "usage: %s <data_dir> <data_name>\n", argv[0]);
    return EXIT_FAILURE;
  }
  char *data_dir = argv[1];
  char *data_name = argv[2];
  /* ensure data_dir ends with '/' */
  size_t len_dir_name = strlen(data_dir);
  if(data_dir[len_dir_name - 1] != '/')
  {
    data_dir = (char*) malloc((len_dir_name + 2) * sizeof(char));
    strcpy(data_dir, argv[1]);
    data_dir[len_dir_name] = '/';
    data_dir[len_dir_name + 1] = '\0';
  }

  /* create data_dir */
  mkdir(data_dir, 0755);

  int fwrites = 0;

  /* Storage space for data on 3 dimensional grid. */
  simType *fields, *fieldsnext, **rks;
  fields    = (simType *) malloc(STORAGE * sizeof(simType));
  fieldsnext  = (simType *) malloc(STORAGE * sizeof(simType));
  /* extra storage for intermediate RK steps */
  rks = (simType **) malloc(RK_STEPS * sizeof(simType *));
  for(n = 0; n < RK_STEPS; n++)
  {
    rks[n] = (simType *) malloc(STORAGE * sizeof(simType));
  }

  /* Preallocated storage space for calculated quantities */
  PointData paq;

  /* calculate intervals to sample at */
  /* interval of steps to sample at */
  const int T_SAMPLEINT = (int) ( (simType)STEPS / (simType)STEPS_TO_RECORD );
  /* position-interval to sample at */
  const int X_SAMPLEINT = (int) ( (simType)POINTS / (simType)POINTS_TO_SAMPLE );


  /* print out sampling information */
  printf("Starting simulation.  Storing data in %s\n", data_dir);
  printf("Will be writing %i of %i steps (sample every %i steps).\n",
    STEPS_TO_RECORD, STEPS, T_SAMPLEINT);
  printf("Writing %i along x-axis (sample every %i points on all axes).\n",
    POINTS_TO_SAMPLE, X_SAMPLEINT);
  printf("Setting w=%1.2f, R_0=%1.2f, (~%1.2f voxels from edge).\n\n",
    W_EOS, R0, POINTS/2-R0/dx);

  /* also write this information to file */
  writeinfo(data_dir, data_name);


  /* initialize data */
  for(i=0; i<POINTS; i++)
    for(j=0; j<POINTS; j++)
      for(k=0; k<POINTS; k++)
      {
        // 0-component is log of energy density
        fields[INDEX(i,j,k,0)] = log(0.1);

        // velocity pattern
        fields[INDEX(i,j,k,1)] = 0.0;
        fields[INDEX(i,j,k,2)] = 0.0;
        fields[INDEX(i,j,k,3)] = 0.0;

        // scalar field
        fields[INDEX(i,j,k,4)] = 0.999057*tanh(sqrt(LAMBDA)*ETA/2*(
                          sqrt(
                            pow( (i*1.0-1.0*POINTS/2.0)*dx , 2)
                            + pow( (j*1.0-1.0*POINTS/2.0)*dx , 2)
                            + pow( (k*1.0-1.0*POINTS/2.0)*dx , 2)
                          ) - R0
                        )
                      ) - 0.025063 ; // spherically symmetric soliton/"bubble" solution

        // time-derivative of scalar field
        fields[INDEX(i,j,k,5)] = 0;
      }


  /* Actual Evolution code */
  for (s=0; s<STEPS; s++)
  {

    // write data if necessary
    if(s % T_SAMPLEINT == 0)
    {
      // fieldsnext is no longer used, so we can re-use it to store an undersampled array of data:
      for(i=0; i<POINTS_TO_SAMPLE; i++)
        for(j=0; j<POINTS_TO_SAMPLE; j++)
          for(k=0; k<POINTS_TO_SAMPLE; k++)
            for(u=0; u<DOF; u++)
              fieldsnext[DOF*POINTS_TO_SAMPLE*POINTS_TO_SAMPLE*i + DOF*POINTS_TO_SAMPLE*j + DOF*k + u]
                = fields[INDEX(X_SAMPLEINT*i, X_SAMPLEINT*j, X_SAMPLEINT*k, u)];

      fflush(stdout);
      printf("\rWriting step %i of %i ...", s, STEPS);
      dumpstate(fieldsnext, fwrites, POINTS_TO_SAMPLE, data_dir, data_name);
      fwrites++;
    } // end write step


    // #pragma omp parallel for default(shared) private(i, j, k, paq) num_threads(2)
    for(i=0; i<POINTS; i++)
      for(j=0; j<POINTS; j++)
        for(k=0; k<POINTS; k++)
        {
          /* Work through normal Euler method. */
          evolve(fields, fieldsnext, 1.0, &paq, i, j, k);

          /* Work through 2nd-order RK method. */
          // midpoint calculation first
          //   evolve(fields, rks[0], 0.5, &paq, i, j, k);
          // final state calculation next
          //   evolve(rks[0], fieldsnext, 1.0, &paq, i, j, k);
        }


    // store new field data
    for(i=0; i<POINTS; i++)
      for(j=0; j<POINTS; j++)
        for(k=0; k<POINTS; k++)
          for(u=0; u<DOF; u++)
          {
#ifdef DEBUG
            /* check for NaN */
            if(isnan(fieldsnext[INDEX(i,j,k,u)]))
            {
              /* do something! */
            }
            else
            {
              fields[INDEX(i,j,k,u)] = fieldsnext[INDEX(i,j,k,u)];
            }
#else
              fields[INDEX(i,j,k,u)] = fieldsnext[INDEX(i,j,k,u)];
#endif
          }

  }

  // dump all data from current simulation... perhaps can read back in later?  Maybe.
  dumpstate(fields, fwrites, POINTS, data_dir, data_name);
  // done.
  printf("Simulation complete.\n");
  printf("%i steps were written to disk.\n", fwrites);
  
  return EXIT_SUCCESS;
}


/* 
 * U^2 - commonly used quantity, stored so we only need to use it once.
 */
simType magu2(PointData *paq)
{
  return paq->fields[1]*paq->fields[1]
      + paq->fields[2]*paq->fields[2]
      + paq->fields[3]*paq->fields[3];
}


/* 
 * U^t - commonly used quantity, stored so we only need to use it once.
 */
simType Ut(PointData *paq)
{
  return sqrt(Ut2(paq));
}


/* 
 * (U^t)^2 - commonly used quantity, stored so we only need to use it once.
 */
simType Ut2(PointData *paq)
{
  return 1 + magu2(paq);
}


/* 
 * Taking derivatives.  Assumes toroidial boundary conditions in each direction.
 */
simType derivative(simType *data, int ddim /* direction of derivative */, int dim, int i, int j, int k)
{
  /* taking modulo here for each point */
  if(ddim == 1) {
    return (data[INDEX((i+1)%POINTS,j,k,dim)] - data[INDEX((i+POINTS-1)%POINTS,j,k,dim)])/2/dx;
  }
  if(ddim == 2) {
    return (data[INDEX(i,(j+1)%POINTS,k,dim)] - data[INDEX(i,(j+POINTS-1)%POINTS,k,dim)])/2/dx;
  }
  if(ddim == 3) {
    return (data[INDEX(i,j,(k+1)%POINTS,dim)] - data[INDEX(i,j,(k+POINTS-1)%POINTS,dim)])/2/dx;
  }

  return 0;
}


/* 
 * Taking second derivatives.  Again, assumes toroidial boundary conditions in each direction.
 */
simType derivative2(simType *data, int ddim /* direction of derivative */, int dim, int i, int j, int k)
{
  if(ddim == 1) {
    return (data[INDEX((i+1)%POINTS,j,k,dim)] + data[INDEX((i+POINTS-1)%POINTS,j,k,dim)] - 2*data[INDEX(i,j,k,dim)])/dx/dx;
  }
  if(ddim == 2) {
    return (data[INDEX(i,(j+1)%POINTS,k,dim)] + data[INDEX(i,(j+POINTS-1)%POINTS,k,dim)] - 2*data[INDEX(i,j,k,dim)])/dx/dx;
  }
  if(ddim == 3) {
    return (data[INDEX(i,j,(k+1)%POINTS,dim)] + data[INDEX(i,j,(k+POINTS-1)%POINTS,dim)] - 2*data[INDEX(i,j,k,dim)])/dx/dx;
  }

  return 0;
}


/*
 * Sum function - spatial sum of two 4-vectors.
 */
simType sumvv(simType v1[4], simType v2[4])
{
  return v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}


/*
 * Overload to give spatial sum of 4-vector and a component of rank-2 tensor.
 */
simType sumvt(simType v1[4], simType t1[4][DOF], int rc /* Sum with row (1) or column (2)? */, int s /* row/col to sum over */)
{
  if(1 == rc) {
    return v1[1] * t1[1][s] + v1[2] * t1[2][s] + v1[3] * t1[3][s];
  }
  if(2 == rc) {
    return v1[1] * t1[s][1] + v1[2] * t1[s][2] + v1[3] * t1[s][3];
  }

  printf("Error: sum/row not specified correctly!\n");
  return 0;
}


/*
 * Overload to give full spatial inner product of 2 vectors and rank-2 tensor.
 * First vector is summed with first tensor index, second with second.
 */
simType sumvtv(simType v1[4], simType t1[4][DOF], simType v2[4])
{
  simType total = 0;
  int i;
  for(i=1; i<=3; i++)
  {
    total += v1[i] * sumvt(v2, t1, 2, i);
  }
  return total;
}

/*
 * Function to take spatial trace of rank-2 tensor.
 */
simType sp_tr(simType t1[4][DOF])
{
  return t1[1][1] + t1[2][2] + t1[3][3];
}


/* 
 * Symmetry breaking potential - \phi^4 with linear perturbation.
 */
simType dV(simType phi)
{
  // return 0;
  return LAMBDA/2*(phi*phi - ETA*ETA)*phi + EPSILON*LAMBDA*ETA*ETA*ETA;
}


/* 
 * Write simulation information to file.
 */
void writeinfo(char *dir, char *root)
{
  char *infofile;
  infofile = (char *) malloc(200 * sizeof(char));
  FILE *datafile;

  /* write simulation parameters to file */
  strcpy(infofile, dir);
  strcat(infofile, root);
  strcat(infofile, ".info");
  datafile = fopen(infofile, "w+");
  if( datafile == NULL )
  {
    printf("Error opening file: %s\n", infofile );
  }
  fprintf(datafile, "# Writing simulation.\n");
  fprintf(datafile, "#Points Sampled:\n%d\n", POINTS_TO_SAMPLE);
  fprintf(datafile, "# Steps recorded:\n%d\n", STEPS_TO_RECORD);
  fprintf(datafile, "# dx/dt (timestep):\n%f", dx/dt);
  fprintf(datafile, "# Equation of state parameter w:\n%f", W_EOS);
  fprintf(datafile, "# Physical lattice dimensions:\n%f %f %f\n",
    SIZE, SIZE, SIZE);
  fprintf(datafile, "# Total number of steps run:\n%d\n", STEPS);

  fclose(datafile);
  free(infofile);
}


/* 
 * Write current simulation state to file.  Should be able to handle arbitrary sized arrays.
 * The "datasize" parameter should be length of a cube.
 */
void dumpstate(simType *fields, int fwrites, int datasize, char *dir, char *name)
{
  char *filename, *buffer;
  filename = (char *) malloc(100 * sizeof(char));
  buffer = (char *) malloc(5 * sizeof(char));

  /* file data for files */
  sprintf(buffer, "%d", fwrites);
  strcpy(filename, dir);
  strcat(filename, name);
  strcat(filename, ".");
  strcat(filename, buffer);
  strcat(filename, ".h5.gz");

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[RANK] = {datasize, datasize, datasize, DOF},
          maxdims[RANK] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED},
          chunk[RANK] = {6, 6, 6, 6};

  file = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (RANK, dims, maxdims);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_deflate (dcpl, 9);
  status = H5Pset_chunk (dcpl, RANK, chunk);
  dset = H5Dcreate2 (file, "DS1", H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fields);

  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

  free(filename);
  free(buffer);
}


// need to write this function.
void readstate(simType *fields, int fwrites)
{

}


/*
 * "source" calculations - compute true fluid source term.
 */
void jsource(PointData *paq)
{
  simType coup = XI*(
          paq->ut * paq->fields[5]
          + sumvt(paq->fields, paq->gradients, 1, 4)
        );

  // little j is source
  paq->ji[0] = coup * paq->fields[5];
  paq->ji[1] = coup * paq->gradients[1][4];
  paq->ji[2] = coup * paq->gradients[2][4];
  paq->ji[3] = coup * paq->gradients[3][4];

  return;
}


/*
 * "source" calculations - compute modified fluid source term, for repeated use.
 */
void Jsource(PointData *paq)
{   
  simType sum_spatial = sumvv(paq->fields, paq->ji);
  simType sum_covariant = sum_spatial + paq->ut * paq->ji[0];
  simType eos_factor = 1.0 + W_EOS;
  simType source_factor = sum_covariant +
    W_EOS / (1.0 + (1.0 - W_EOS) * paq->u2) * (sum_spatial / eos_factor + paq->u2 * sum_covariant);

  int i;
  for(i=1; i<=3; i++) {
    paq->Ji[i] = (
        paq->ji[i] / eos_factor
        + paq->fields[i] * source_factor
      ) / exp(paq->fields[0]) / paq->ut;
  }

  return;
}


/*
 * DOF evolution - log(energy density) evolution.
 */
simType energy_evfn(PointData *paq)
{
  return (
    - W_EOSp1 * paq->ut / (1.0 - W_EOSm1 * paq->u2) * (
      -W_EOSm1 / W_EOSp1 * sumvt(paq->fields, paq->gradients, 1, 0)
      + sp_tr(paq->gradients)
      - paq->uudu / paq->ut2
    )
    - W_EOSp1 / paq->ut2 * sumvv(paq->fields, paq->Ji)
    - (paq->ut * paq->ji[0] + sumvv(paq->fields, paq->ji)) / exp(paq->fields[0]) / paq->ut
  );
}


/*
 * DOF evolution - fluid velocity component evolution.
 */
simType fluid_evfn(PointData *paq, int u)
{
  return (
    W_EOS * paq->ut * paq->fields[u] / (1.0 - W_EOSm1 * paq->u2 ) * (
        sp_tr(paq->gradients)
        - (
          W_EOS * sumvt(paq->fields, paq->gradients, 1, 0) / W_EOSp1
          + paq->uudu
        ) / paq->ut2
      )
    - (
      sumvt(paq->fields, paq->gradients, 1, u) + W_EOS / W_EOSp1 * paq->gradients[u][0]
      ) / paq->ut
    + paq->Ji[u]
  );
}


/*
 * DOF evolution - scalar field.
 */
simType field_evfn(PointData *paq)
{
  return paq->fields[5];
}


/*
 * DOF evolution - time derivative of field (w ~ d\phi/dt).
 */
simType ddtfield_evfn(PointData *paq)
{
  return (
    paq->derivs2[1] + paq->derivs2[2] + paq->derivs2[3]
    - XI*(
      paq->ut * paq->fields[5]
      + sumvt(paq->fields, paq->gradients, 1, 4)
    ) - dV(paq->fields[4])
  );
}


/*
 * Calculate commonly used quantities at each point in the simulation.  Quantities are contained within a
 * struct - ideally this will not be any faster or slower than explicitly writing out each quantity.
 */
void calculatequantities(simType *fields, PointData *paq, int i, int j, int k)
{
  int u, n;

  // Field data at pertinent point
  for(n=0; n<=5; n++) {
    paq->fields[n] = fields[INDEX(i,j,k,n)];
  }

  // [COMMON_QUANTITIES]
  paq->u2 = magu2(paq);
  paq->ut = Ut(paq);
  paq->ut2 = Ut2(paq);
  paq->uudu = sumvtv(paq->fields, paq->gradients, paq->fields);

  // [GRADIENTS]
  for(u=0; u<DOF; u++) {
    n = 0;
     for(n=1; n<=3; n++)
       paq->gradients[n][u] = derivative(fields, n, u, i, j, k);
  }

  paq->derivs2[1] = derivative2(fields, 1, 4, i, j, k);
  paq->derivs2[2] = derivative2(fields, 2, 4, i, j, k);
  paq->derivs2[3] = derivative2(fields, 3, 4, i, j, k);

  // compute source term information [SOURCE]
  // little j is source, big J is some wonko function of source
  jsource(paq);
  Jsource(paq);

  return;
}


/*
 * Evolution step of simulation - compute next step.
 */
void evolve(simType *initial, simType *final, simType coeff, PointData *paq, int i, int j, int k)
{
  calculatequantities(initial, paq, i, j, k);

  // energy density
   final[INDEX(i,j,k,0)] = paq->fields[0] + coeff*dt*energy_evfn(paq);
  // fluid
   final[INDEX(i,j,k,1)] = paq->fields[1] + coeff*dt*fluid_evfn(paq, 1);
   final[INDEX(i,j,k,2)] = paq->fields[2] + coeff*dt*fluid_evfn(paq, 2);
   final[INDEX(i,j,k,3)] = paq->fields[3] + coeff*dt*fluid_evfn(paq, 3);
  // field
   final[INDEX(i,j,k,4)] = paq->fields[4] + coeff*dt*field_evfn(paq);
  // field derivative
   final[INDEX(i,j,k,5)] = paq->fields[5] + coeff*dt*ddtfield_evfn(paq);
}
