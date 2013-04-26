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
        // setXI(atof(optarg));
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

  /* Preallocated storage space for calculated quantities */
  PointData paq;

  /* Fluid/field storage space for data on 3 dimensional grid. */
  simType *fields, *fieldsnext, **rks;
  fields      = (simType *) malloc(STORAGE * ((long long) sizeof(simType)));
  fieldsnext  = (simType *) malloc(STORAGE * ((long long) sizeof(simType)));
  /* extra storage for intermediate RK steps */
  rks = (simType **) malloc(RK_STEPS * sizeof(simType *));
  for(n = 0; n < RK_STEPS; n++)
  {
    rks[n] = (simType *) malloc(STORAGE * ((long long) sizeof(simType)));
  }

/* Gravitational perturbation storage */
// let's just evolve 2 degrees of freedom for now
//  simType *hij, *lij, *STTij;
//  hij = (simType *) malloc(GRID_STORAGE * 2 * sizeof(simType));
//  lij = (simType *) malloc(GRID_STORAGE * 2 * sizeof(simType));
//  STTij = (simType *) malloc(GRID_STORAGE * sizeof(simType));

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

  if(0 == read_initial_step)
  {
    /* initialize data in static bubble configuration */
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
                  // Bubble is slightly offset from the grid center here if
                  // there are an even number of grid points - it will be lined
                  // up with a pixel.  But this is ok, even desirable if we
                  // want to look at a slice through the middle of the bubble.
                  pow( (i*1.0-1.0*POINTS/2.0)*dx , 2)
                  + pow( (j*1.0-1.0*POINTS/2.0)*dx , 2)
                  + pow( (k*1.0-1.0*POINTS/2.0)*dx , 2)
                ) - R0
              )
            ) - 0.025063 ; // spherically symmetric soliton/"bubble" solution

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

    // write data if necessary
    if(s % T_SAMPLEINT == 0)
    {
      // fieldsnext is not being used yet, so we can use it to store an undersampled array of data:
      for(i=0; i<POINTS_TO_SAMPLE; i++)
        for(j=0; j<POINTS_TO_SAMPLE; j++)
          for(k=0; k<POINTS_TO_SAMPLE; k++)
            for(u=0; u<DOF; u++)
              fieldsnext[DOF*POINTS_TO_SAMPLE*POINTS_TO_SAMPLE*i + DOF*POINTS_TO_SAMPLE*j + DOF*k + u]
                = fields[INDEX(X_SAMPLEINT*i, X_SAMPLEINT*j, X_SAMPLEINT*k, u)];

      //printf("\rWriting step %i of %i ...", s, MAX_STEPS);
      //fflush(stdout);
      filedata.datasize = POINTS_TO_SAMPLE;
      dumpstate(fieldsnext, filedata);
      filedata.fwrites++;
    } // end write step

    if(DUMP_STRIP)
    {
      dumpstrip(fields, filedata);
    }

    // write full dump and Hartley Transform data if appropriate
    if(s % T_DUMPINT == 0)
    {
      /* fieldsnext isn't being used, so we can pass it in as storage */
      //printf("\nDumping HT for step %i of %i ...\n", s, MAX_STEPS);
      //fflush(stdout);
      hartleydump(fields, fieldsnext, filedata);
    }

    #pragma omp parallel for default(shared) private(i, j, k, paq) num_threads(threads)
    for(i=0; i<POINTS; i++)
    {
      for(j=0; j<POINTS; j++)
      {
        for(k=0; k<POINTS; k++)
        {
          /* Work through normal Euler method. */
          // evolve(fields, fieldsnext, fields, 1.0, &paq, i, j, k);
          // gw_evolve(hij, lij, STTij, &paq, i, j, k);

          /* Work through 2nd-order RK method. */
          // midpoint calculation first
            evolve(fields, rks[0], fields, 0.5, &paq, i, j, k);
          // final state calculation next
            evolve(fields, fieldsnext, rks[0], 1.0, &paq, i, j, k);

        }
      }
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
    paq->derivs2[1] + paq->derivs2[2] + paq->derivs2[3]
    + getXI() * (
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
  paq->relw = (1.0 - W_EOSm1 * paq->u2);

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

  paq->uudu = sumvtv(paq->fields, paq->gradients, paq->fields);
  paq->srcsum = sumvv(paq->fields, paq->ji);
  paq->trgrad = sp_tr(paq->gradients);
  paq->udu = sumvt(paq->fields, paq->gradients, 1, 0);

  return;
}


/*
 * Evolution step of simulation - compute next step.
 */
void evolve(simType *initial, simType *final, simType *intermediate, simType coeff, PointData *paq,
  int i, int j, int k)
{
  calculatequantities(intermediate, paq, i, j, k);

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
