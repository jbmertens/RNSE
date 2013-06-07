/* RNSE includes */
#include "defines.h"


int main(int argc, char **argv)
{
  /* iterators here first. */
  int i, j, k, n, u, s;

  /* information about files, default options */
  IOData filedata;
  filedata.fwrites = 0;
  filedata.data_dir = DEFAULT_DATA_DIR;
  filedata.data_name = DEFAULT_DATA_NAME;
  int read_initial_step = 0;
  int threads = 2;

  /* read in and process non-default options */
  /* -c coupling_xi -o output_dir -f output_filename -i initial_configuration -t num_threads */
  int c = 0;
  while(1)
  {
    static struct option long_options[] =
    {
      {"coupling",      required_argument, 0, 'c'},
      {"output-dir",    required_argument, 0, 'o'},
      {"output-file",   required_argument, 0, 'f'},
      {"initial-file",  required_argument, 0, 'i'},
      {"num-threads",   required_argument, 0, 't'},
      {"num-threads",   no_argument,       0, 'h'}
    };
    
    int option_index = 0;
    c = getopt_long(argc, argv, "c:o:f:i:t:h", long_options, &option_index);
    if(c == -1) // stop if done reading arguments
      break;

    switch(c)
    {
      case 'c':
        setXI(atof(optarg));
        break;
      case 'o':
        filedata.data_dir = optarg;
        break;
      case 'f':
        filedata.data_name = optarg;
        break;
      case 'i':
        filedata.read_data_name = optarg;
        read_initial_step = 1;
        break;
      case 't':
        threads = atof(optarg);
        break;
      case 'h':
      case '?':
        fprintf(stderr, "usage: %s -c coupling_xi -o output_dir -f output_filename -i initial_configuration -t num_threads\n", argv[0]);
      default:
        abort();
    }
  }

  /* ensure data_dir ends with '/', unless empty string is specified. */
  size_t len_dir_name = strlen(filedata.data_dir);
  if(filedata.data_dir[len_dir_name - 1] != '/' && len_dir_name != 0)
  {
    char *data_dir;
    data_dir = (char*) malloc((len_dir_name + 2) * sizeof(char));
    strcpy(data_dir, filedata.data_dir);
    filedata.data_dir = data_dir;
    filedata.data_dir[len_dir_name] = '/';
    filedata.data_dir[len_dir_name + 1] = '\0';
  }
  /* create data_dir */
  if(len_dir_name != 0)
    mkdir(filedata.data_dir, 0755);

  /* some validation for sampling rates */
  if(MAX_STEPS < STEPS_TO_SAMPLE || POINTS < POINTS_TO_SAMPLE || MAX_STEPS < STEPS_TO_DUMP)
  {
    fprintf(stderr, "# of samples should be larger than total number of steps/points.");
    return EXIT_FAILURE;
  }

  /* calculate intervals to sample at */
  int T_SAMPLEINT;
  int X_SAMPLEINT;
  int T_DUMPINT;
  if(STEPS_TO_SAMPLE == 0) {
    T_SAMPLEINT = MAX_STEPS+1;
  } else {
    T_SAMPLEINT = (int) floor( (simType)MAX_STEPS / (simType)STEPS_TO_SAMPLE );
  }
  if(POINTS_TO_SAMPLE == 0) {
    X_SAMPLEINT = POINTS+1;
  } else {
    X_SAMPLEINT = (int) floor( (simType)POINTS / (simType)POINTS_TO_SAMPLE );
  }
  if(STEPS_TO_DUMP == 0) {
    T_DUMPINT = MAX_STEPS+1;
  } else {
    T_DUMPINT = (int) floor( (simType)MAX_STEPS / (simType)STEPS_TO_DUMP );
  }

  /* print out sampling information */
  printf("\nStarting simulation.  Storing data in %s\n", filedata.data_dir);
  printf("Will be sampling every %i steps (recording about %i of %i steps).\n",
    T_SAMPLEINT, MAX_STEPS/(T_SAMPLEINT+1), MAX_STEPS);
  printf("Full dump output every %i steps (recording about %i of %i steps).\n",
    T_DUMPINT, MAX_STEPS/(T_DUMPINT+1), MAX_STEPS);
  printf("Writing %i along x-axis (sample every %i points on all axes).\n",
    POINTS_TO_SAMPLE, X_SAMPLEINT);
  printf("Setting w=%1.2f, R_0=%1.2f, (~%1.2f voxels from edge).\n\n",
    W_EOS, R0, POINTS/2-R0/dx);

  /* also write this information to file */
  writeinfo(filedata);


  /* Preallocated storage space for calculated quantities */
  PointData paq;

  /* Fluid/field storage space for data on 3 dimensional grid. */
  simType *fields, *wedge, *after;
  // actual grid
  fields     = (simType *) malloc(STORAGE * ((long long) sizeof(simType)));
  // For the "wedge", R^2 storage locations are needed.
  // wedge[i=0-2] are 'base', wedge[i=3] is 'peak'
  wedge      = (simType *) malloc(AREA_STORAGE * METHOD_ORDER * METHOD_ORDER * ((long long) sizeof(simType)));
  // Initial values calculated in the wedge can't be stored until it's come back and finished calculating
  // Basically, a 'snapshot' of the first two peaks - call this the 'afterimage'.
  after = (simType *) malloc(AREA_STORAGE * 2 * METHOD_ORDER * ((long long) sizeof(simType)));


  if(0 == read_initial_step)
  {
    /* initialize data in static bubble configuration */
    for(i=0; i<POINTS; i++)
      for(j=0; j<POINTS; j++)
        for(k=0; k<POINTS; k++)
        {
          // 0-component is log of energy density
          fields[INDEX(i,j,k,0)] = log(0.1);

          // velocity pattern - at rest
          fields[INDEX(i,j,k,1)] = 0.0;
          fields[INDEX(i,j,k,2)] = 0.0;
          fields[INDEX(i,j,k,3)] = 0.0;

          // scalar field
          /**
          fields[INDEX(i,j,k,4)] = tanh(sqrt(LAMBDA)*ETA/2*(
                sqrt(
                  // Bubble is slightly offset from the grid center here - if
                  // there are an even number of grid points. it will be lined
                  // up with a pixel.  But this is ok, even desirable if we
                  // want to look at a slice through the middle of the bubble.
                  pow( (i*1.0-POINTS/2.0)*dx , 2)
                  + pow( (j*1.0-POINTS/2.0)*dx , 2)
                  + pow( (k*1.0-POINTS/2.0)*dx , 2)
                ) - R0)
              );
          /**/
          fields[INDEX(i,j,k,4)] = 0.999057*tanh(sqrt(LAMBDA)*ETA/2*(
                sqrt(
                  // Bubble is slightly offset from the grid center here - if
                  // there are an even number of grid points. it will be lined
                  // up with a pixel.  But this is ok, even desirable if we
                  // want to look at a slice through the middle of the bubble.
                  pow( (i*1.0-POINTS/2.0)*dx , 2)
                  + pow( (j*1.0-POINTS/2.0)*dx , 2)
                  + pow( (k*1.0-POINTS/2.0)*dx , 2)
                ) - R0
              )
            ) - 0.025063; // spherically symmetric soliton/"bubble" solution
          /**/

          // time-derivative of scalar field
          fields[INDEX(i,j,k,5)] = 0;
        }
  }
  else
  {
    readstate(fields, filedata);
  }


  /* record simulation time / wall clock time */
  time_t time_start = time(NULL);

  /* Actual Evolution code */
  for (s=1; s<=MAX_STEPS; s++)
  {
/* TODO:
    // write data if necessary
    if(s % T_SAMPLEINT == 0)
    {
      // fieldsnext is not being used yet, so we can use it to store an undersampled array of data:
      LOOP4(i,j,k,u)
        fieldsnext[DOF*POINTS_TO_SAMPLE*POINTS_TO_SAMPLE*i + DOF*POINTS_TO_SAMPLE*j + DOF*k + u]
          = fields[INDEX(X_SAMPLEINT*i, X_SAMPLEINT*j, X_SAMPLEINT*k, u)];

      //printf("\rWriting step %i of %i ...", s, MAX_STEPS);
      //fflush(stdout);
      filedata.datasize = POINTS_TO_SAMPLE;
      dumpstate(fieldsnext, filedata);
      filedata.fwrites++;
    } // end write step
/**/

    if(DUMP_STRIP)
    {
      dumpstrip(fields, filedata);
    }


// TODO:
    // // write full dump and Hartley Transform data if appropriate
    // if(s % T_DUMPINT == 0)
    // {
    //    fieldsnext isn't being used, so we can pass it in as storage 
    //   //printf("\nDumping HT for step %i of %i ...\n", s, MAX_STEPS);
    //   //fflush(stdout);
    //   hartleydump(fields, fieldsnext, filedata);
    // }


/****
 * Memory-efficient midpoint implementation, hopefully at no significant speed cost.
 * Requires several area grids, forming a "wedge", and an initial snapshot of the
 * first two calculations, an "afterimage".
 * For the midpoint method:
 *   - wedge stores 4 grids (first 3 are 'base', last is the 'peak')
 *   - afterimage holds first 2 peak calculations
 * Parallelize on the j, k loops (area calculations).
 */

    // initial wedge base (centered on i=0)
      for(i=POINTS-1; i<=POINTS+1; i++) {
        // populate wedge
        #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
        LOOP2(j,k)
          g2wevolve(fields, wedge, &paq, i, j, k);
      }

    // build afterimage
      // store first peak afterimage
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
      {
        w2pevolve(wedge, &paq, 0, j, k);
        for(u=0; u<DOF; u++) {
          after[INDEX(0,j,k,u)] = wedge[INDEX(3,j,k,u)];
        }
      }
      // roll wedge base forward
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
        g2wevolve(fields, wedge, &paq, 2, j, k);
      // second afterimage point
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
      {
        w2pevolve(wedge, &paq, 1, j, k);
        for(u=0; u<DOF; u++) {
          after[INDEX(1,j,k,u)] = wedge[INDEX(3,j,k,u)];
        }
      }

    // finally: starting wedge
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
        g2wevolve(fields, wedge, &paq, 3, j, k);
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
      {
        w2pevolve(wedge, &paq, 2, j, k);
      }

    // Move along, move along home
    for(i=4; i<=POINTS; i++) // i is position of leading wedge base point; i=1 is position of peak
    {
      
      // roll base along
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
        g2wevolve(fields, wedge, &paq, i, j, k);

      // grid <- peak <- base
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
      {
        for(u=0; u<DOF; u++)
        {
          fields[INDEX(i-2,j,k,u)] = wedge[INDEX(3,j,k,u)];
        }
        w2pevolve(wedge, &paq, i-1, j, k);
      }
    }

    // Done; store the afterimage in the fields.
    #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
    LOOP2(j,k)
      for(u=0; i<DOF; u++)
      {
        fields[INDEX(0,j,k,u)] = after[INDEX(0,j,k,u)];
        fields[INDEX(1,j,k,u)] = after[INDEX(1,j,k,u)];
      }


/** End wedge method **/

    if(STOP_CELL > 0 && STOP_CELL < POINTS && fields[INDEX( POINTS/2, POINTS/2, STOP_CELL, 4)] < 0)
    {
      printf("Bubble wall has hit stop condition - ending simulation.\n");
      break;
    }

  }

  time_t time_end = time(NULL);

  // At end, dump all data from current simulation.
  filedata.datasize = POINTS;
  dumpstate(fields, filedata);
  // done.
  printf("Simulation complete.\n");
  printf("%i steps were written to disk.\n", filedata.fwrites);
  printf("Simulation took %ld seconds.\n", (long)(time_end - time_start));
  
  return EXIT_SUCCESS;
}


/*
 * "source" calculations - compute true fluid source term.
 */
void jsource(PointData *paq)
{
  simType coup = - getXI() * (
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
  simType sum_spatial = paq->srcsum;
  simType sum_covariant = sum_spatial + paq->ut * paq->ji[0];
  simType eos_factor = 1.0 + W_EOS;
  simType source_factor = sum_covariant +
    W_EOS / paq->relw * (sum_spatial / eos_factor + paq->u2 * sum_covariant);

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
static inline simType energy_evfn(PointData *paq)
{
  return (
    - W_EOSp1 * paq->ut / paq->relw * (
      -W_EOSm1 / W_EOSp1 * paq->udu
      + paq->trgrad
      - paq->uudu / paq->ut2
    )
    - W_EOSp1 / paq->ut2 * sumvv(paq->fields, paq->Ji)
    - (paq->ut * paq->ji[0] + paq->srcsum) / exp(paq->fields[0]) / paq->ut
  );
}


/*
 * DOF evolution - fluid velocity component evolution.
 */
static inline simType fluid_evfn(PointData *paq, int u)
{
  return (
    W_EOS * paq->ut * paq->fields[u] / paq->relw * (
        paq->trgrad
        - (
          W_EOS * paq->udu / W_EOSp1
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
static inline simType field_evfn(PointData *paq)
{
  return paq->fields[5];
}


/*
 * DOF evolution - time derivative of field (w ~ d\phi/dt).
 */
static inline simType ddtfield_evfn(PointData *paq)
{
  return (
    paq->lap + getXI() * (
      paq->ut * paq->fields[5]
      + sumvt(paq->fields, paq->gradients, 1, 4)
    ) - dV(paq->fields[4])
  );
}


/**
 * Note: these functions are now specific to the midpoint method.
 * Roll back to 0f66228 for old, non-wedge code.
 */

/*
 * Calculate commonly used quantities at each point in the simulation.  Quantities are contained within a
 * struct - ideally this will not be any faster or slower than explicitly writing out each quantity.
 * Following this, calculate the next step.
 */

// version evolving from grid to wedge base
void g2wevolve(simType *grid, simType *wedge, PointData *paq, int i, int j, int k)
{
  int u, n;

  // Field data at pertinent point
  for(n=0; n<=5; n++) {
    paq->fields[n] = grid[INDEX(i,j,k,n)];
  }

  // [COMMON QUANTITIES]
  paq->u2 = magu2(paq);
  paq->ut = Ut(paq);
  paq->ut2 = Ut2(paq);
  paq->relw = (1.0 - W_EOSm1 * paq->u2);

  // [GRADIENTS]
  for(u=0; u<5; u++) {
    n = 0;
    for(n=1; n<=3; n++)
      paq->gradients[n][u] = derivative(grid, n, u, i, j, k);
  }
  paq->lap = lapl(grid, 4, i, j, k);

  // [SOURCES] little j is source, big J is some wonko function of source
  jsource(paq);
  Jsource(paq);

  // [DERIVED QUANTITIES]
  paq->uudu = sumvtv(paq->fields, paq->gradients, paq->fields);
  paq->srcsum = sumvv(paq->fields, paq->ji);
  paq->trgrad = sp_tr(paq->gradients);
  paq->udu = sumvt(paq->fields, paq->gradients, 1, 0);

  // [EVOLVE GRID TO WEDGE BASE]
  // energy density
   wedge[INDEX(i%3,j,k,0)] = grid[INDEX(i,j,k,0)] + dt*energy_evfn(paq)/2.0;
  // fluid
   wedge[INDEX(i%3,j,k,1)] = grid[INDEX(i,j,k,1)] + dt*fluid_evfn(paq, 1)/2.0;
   wedge[INDEX(i%3,j,k,2)] = grid[INDEX(i,j,k,2)] + dt*fluid_evfn(paq, 2)/2.0;
   wedge[INDEX(i%3,j,k,3)] = grid[INDEX(i,j,k,3)] + dt*fluid_evfn(paq, 3)/2.0;
  // field
   wedge[INDEX(i%3,j,k,4)] = grid[INDEX(i,j,k,4)] + dt*field_evfn(paq)/2.0;
  // field derivative
   wedge[INDEX(i%3,j,k,5)] = grid[INDEX(i,j,k,5)] + dt*ddtfield_evfn(paq)/2.0;
}


// version evolving from wedge base to wedge peak
void w2pevolve(simType *wedge, PointData *paq, int i, int j, int k)
{
  int u, n;

  // Field data at pertinent point
  for(n=0; n<=5; n++) {
    paq->fields[n] = wedge[INDEX(i%3,j,k,n)];
  }

  // [COMMON QUANTITIES]
  paq->u2 = magu2(paq);
  paq->ut = Ut(paq);
  paq->ut2 = Ut2(paq);
  paq->relw = (1.0 - W_EOSm1 * paq->u2);

  // [GRADIENTS]
  for(u=0; u<5; u++) {
    n = 0;
    for(n=1; n<=3; n++)
      paq->gradients[n][u] = wderivative(wedge, n, u, i, j, k);
  }
  paq->lap = wlapl(wedge, 4, i, j, k);

  // [SOURCES] little j is source, big J is some wonko function of source
  jsource(paq);
  Jsource(paq);

  // [DERIVED QUANTITIES]
  paq->uudu = sumvtv(paq->fields, paq->gradients, paq->fields);
  paq->srcsum = sumvv(paq->fields, paq->ji);
  paq->trgrad = sp_tr(paq->gradients);
  paq->udu = sumvt(paq->fields, paq->gradients, 1, 0);

  // [EVOLVE WEDGE BASE TO PEAK]
  // energy density
   wedge[INDEX(3,j,k,0)] = wedge[INDEX(i%3,j,k,0)] + dt*energy_evfn(paq);
  // fluid
   wedge[INDEX(3,j,k,1)] = wedge[INDEX(i%3,j,k,1)] + dt*fluid_evfn(paq, 1);
   wedge[INDEX(3,j,k,2)] = wedge[INDEX(i%3,j,k,2)] + dt*fluid_evfn(paq, 2);
   wedge[INDEX(3,j,k,3)] = wedge[INDEX(i%3,j,k,3)] + dt*fluid_evfn(paq, 3);
  // field
   wedge[INDEX(3,j,k,4)] = wedge[INDEX(i%3,j,k,4)] + dt*field_evfn(paq);
  // field derivative
   wedge[INDEX(3,j,k,5)] = wedge[INDEX(i%3,j,k,5)] + dt*ddtfield_evfn(paq);

}
