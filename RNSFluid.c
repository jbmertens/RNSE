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
      {"alpha",         required_argument, 0, 'a'},
      {"num-threads",   no_argument,       0, 'h'}
    };
    
    int option_index = 0;
    c = getopt_long(argc, argv, "c:o:f:i:t:a:h", long_options, &option_index);
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
      case 'a':
        setALPHA(atof(optarg));
        break;
      case 'h':
      case '?':
        fprintf(stderr, "usage: %s -c coupling_xi -o output_dir -f output_filename -i initial_configuration -t num_threads\n", argv[0]);
        abort();
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

  /* print out simulation information */
  printf("\nStarting simulation.  Storing data in %s\n", filedata.data_dir);
  if(T_SAMPLEINT < MAX_STEPS) {
    printf("Will be sampling every %i steps (recording about %i of %i steps).\n",
      T_SAMPLEINT, MAX_STEPS/(T_SAMPLEINT+1), MAX_STEPS);
    printf(" - writing %i along x-axis (sample every %i points on all axes).\n",
      POINTS_TO_SAMPLE, X_SAMPLEINT);
  } else {
    printf("Will not be outputting grid snapshots.");
  }
  if(T_DUMPINT < MAX_STEPS) {
    printf("Full dump output every %i steps (recording about %i of %i steps).\n",
      T_DUMPINT, MAX_STEPS/(T_DUMPINT+1), MAX_STEPS);
  } else {
    printf("Will not be outputting full grid dumps except at end.\n"); 
  }
  printf("Parameter values: alpha=%1.2f, xi=%1.4f, w=%1.2f, R_0=%1.2f.\n",
    getALPHA(), getXI(), W_EOS, R0);
  printf("Simulation information: size=%d*R_0, points=%d, dx=%1.3f, dt/dx=%1.3f.\n\n",
    ((int) (SIZE/R0)), ((int) POINTS), dx, dt/dx);

  /* also write this information to file */
  writeinfo(filedata);


  /* Preallocated storage space for calculated quantities */
  PointData paq;

  /* Fluid/field storage space for data on 3 dimensional grid. */
  simType *fields, *wedge, *after, *transfer;
  // actual grid
  fields     = (simType *) malloc(STORAGE * ((long long) sizeof(simType)));
  // For the "wedge", R^2 storage locations are needed.
  // wedge[i=0-2] are 'base', wedge[i=3] is 'peak'
  wedge      = (simType *) malloc(AREA_STORAGE * METHOD_ORDER * METHOD_ORDER * ((long long) sizeof(simType)));
  // Overhead from transfering wedge values back to grid
  transfer   = (simType *) malloc(AREA_STORAGE * ((long long) sizeof(simType)));
  // Initial values calculated in the wedge can't be stored until it's come back and finished calculating
  // Basically, a 'snapshot' of the first two peaks - call this the 'afterimage'.
  after = (simType *) malloc(AREA_STORAGE * 2 * METHOD_ORDER * ((long long) sizeof(simType)));


  if(0 == read_initial_step)
  {
    /* initialize data in static bubble configuration */
    printf("Initializing bubble configuration using thin-wall approximation.\n");
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
          // spherically symmetric soliton/"bubble" solution - first order
          // approximation in vacuum energy difference
          fields[INDEX(i,j,k,4)] = (3.0 + sqrt(9.0 - 8.0*getALPHA()))/4.0/getALPHA()
            *(
              tanh(1.0/2.0*(
                sqrt(
                  pow( (i*1.0-((double) POINTS)/2.0)*dx , 2)
                  + pow( (j*1.0-((double) POINTS)/2.0)*dx , 2)
                  + pow( (k*1.0-((double) POINTS)/2.0)*dx , 2)
                ) - R0
              )) - 1.0
            );

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
      // populate wedge
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k) {
        g2wevolve(fields, wedge, &paq, POINTS-1, j, k);
        g2wevolve(fields, wedge, &paq, POINTS, j, k);
        g2wevolve(fields, wedge, &paq, POINTS+1, j, k);
      }

    // build afterimage
      // store first peak afterimage
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
      {
        w2pevolve(fields, wedge, &paq, 0, j, k);
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
        w2pevolve(fields, wedge, &paq, 1, j, k);
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
        w2pevolve(fields, wedge, &paq, 2, j, k);
      }

    // Move along, move along home
    for(i=4; i<=POINTS; i++) // i is position of leading wedge base point; i-1 is position of new peak
    {
      
      // roll base along
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
        g2wevolve(fields, wedge, &paq, i, j, k);

      // store old peak back in the grid
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
        for(u=0; u<DOF; u++)
          fields[INDEX(i-2,j,k,u)] = wedge[INDEX(3,j,k,u)];

      // calculate new wedge peak
      #pragma omp parallel for default(shared) private(j, k, paq) num_threads(threads)
      LOOP2(j,k)
        w2pevolve(fields, wedge, &paq, i-1, j, k);

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

    if(STOP_CELL > 0 && STOP_CELL < POINTS && fields[INDEX( POINTS/2, POINTS/2, STOP_CELL, 4)] < STOP_MAX)
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
  simType coup = -1.0*getXI() * (
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
  simType sum_covariant = paq->srcsum + paq->ut * paq->ji[0];
  simType source_factor = sum_covariant +
    W_EOS / paq->relw * (paq->srcsum / W_EOSp1 + paq->u2 * sum_covariant);

  int i;
  for(i=1; i<=3; i++) {
    paq->Ji[i] = (
        paq->ji[i] / W_EOSp1
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
      - W_EOSm1 / W_EOSp1 * paq->ude
      + paq->trgrad
      + paq->uudu / paq->ut2
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
          W_EOS * paq->ude / W_EOSp1
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
    paq->lap - getXI() * (
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
  int n;

  // Field data near working point
  for(n=0; n<DOF; n++)
    paq->fields[n] = grid[INDEX(i,j,k,n)];
  // move values from heap to stack (cut cache misses)
  for(n=0; n<DOF; n++)
    paq->adjacentFields[0][n] = grid[INDEX((i+1)%POINTS,j,k,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[3][n] = grid[INDEX((i+POINTS-1)%POINTS,j,k,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[1][n] = grid[INDEX(i,(j+1)%POINTS,k,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[4][n] = grid[INDEX(i,(j+POINTS-1)%POINTS,k,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[2][n] = grid[INDEX(i,j,(k+1)%POINTS,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[5][n] = grid[INDEX(i,j,(k+POINTS-1)%POINTS,n)];

  // Around edges (for laplacian; field values only)
  paq->adjacentEdges[0] = grid[INDEX((i+1)%POINTS,(j+1)%POINTS,k,4)];
  paq->adjacentEdges[1] = grid[INDEX((i+1)%POINTS,(j-1+POINTS)%POINTS,k,4)];
  paq->adjacentEdges[2] = grid[INDEX((i+1)%POINTS,j,(k+1)%POINTS,4)];
  paq->adjacentEdges[3] = grid[INDEX((i+1)%POINTS,j,(k-1+POINTS)%POINTS,4)];
  paq->adjacentEdges[4] = grid[INDEX((i-1+POINTS)%POINTS,(j+1)%POINTS,k,4)];
  paq->adjacentEdges[5] = grid[INDEX((i-1+POINTS)%POINTS,(j-1+POINTS)%POINTS,k,4)];
  paq->adjacentEdges[6] = grid[INDEX((i-1+POINTS)%POINTS,j,(k+1)%POINTS,4)];
  paq->adjacentEdges[7] = grid[INDEX((i-1+POINTS)%POINTS,j,(k-1+POINTS)%POINTS,4)];
  paq->adjacentEdges[8] = grid[INDEX(i,(j+1)%POINTS,(k+1)%POINTS,4)];
  paq->adjacentEdges[9] = grid[INDEX(i,(j+1)%POINTS,(k-1+POINTS)%POINTS,4)];
  paq->adjacentEdges[10] = grid[INDEX(i,(j-1+POINTS)%POINTS,(k+1)%POINTS,4)];
  paq->adjacentEdges[11] = grid[INDEX(i,(j-1+POINTS)%POINTS,(k-1+POINTS)%POINTS,4)];

  // calculate quantities used by evolution functions
  calculatequantities(paq);

  // [EVOLVE GRID TO WEDGE BASE]
  // energy density
   wedge[WINDEX(i,j,k,0)] = grid[INDEX(i,j,k,0)] + dt*0.5*energy_evfn(paq);
  // fluid
   wedge[WINDEX(i,j,k,1)] = grid[INDEX(i,j,k,1)] + dt*0.5*fluid_evfn(paq, 1);
   wedge[WINDEX(i,j,k,2)] = grid[INDEX(i,j,k,2)] + dt*0.5*fluid_evfn(paq, 2);
   wedge[WINDEX(i,j,k,3)] = grid[INDEX(i,j,k,3)] + dt*0.5*fluid_evfn(paq, 3);
  // field
   wedge[WINDEX(i,j,k,4)] = grid[INDEX(i,j,k,4)] + dt*0.5*field_evfn(paq);
  // field derivative
   wedge[WINDEX(i,j,k,5)] = grid[INDEX(i,j,k,5)] + dt*0.5*ddtfield_evfn(paq);

}


// version evolving from wedge base to wedge peak
void w2pevolve(simType *grid, simType *wedge, PointData *paq, int i, int j, int k)
{
  int n;

  // Field data near working point
  for(n=0; n<DOF; n++)
    paq->fields[n] = wedge[WINDEX(i,j,k,n)];
  // adjacent to point (cut cache misses)
  for(n=0; n<DOF; n++)
    paq->adjacentFields[0][n] = wedge[WINDEX((i+1)%POINTS,j,k,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[3][n] = wedge[WINDEX((i-1+POINTS)%POINTS,j,k,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[1][n] = wedge[WINDEX(i,(j+1)%POINTS,k,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[4][n] = wedge[WINDEX(i,(j-1+POINTS)%POINTS,k,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[2][n] = wedge[WINDEX(i,j,(k+1)%POINTS,n)];
  for(n=0; n<DOF; n++)
    paq->adjacentFields[5][n] = wedge[WINDEX(i,j,(k-1+POINTS)%POINTS,n)];

  // Around edges (for laplacian; field values only)
  paq->adjacentEdges[0] = wedge[WINDEX((i+1)%POINTS,(j+1)%POINTS,k,4)];
  paq->adjacentEdges[1] = wedge[WINDEX((i+1)%POINTS,(j-1+POINTS)%POINTS,k,4)];
  paq->adjacentEdges[2] = wedge[WINDEX((i+1)%POINTS,j,(k+1)%POINTS,4)];
  paq->adjacentEdges[3] = wedge[WINDEX((i+1)%POINTS,j,(k-1+POINTS)%POINTS,4)];
  paq->adjacentEdges[4] = wedge[WINDEX((i-1+POINTS)%POINTS,(j+1)%POINTS,k,4)];
  paq->adjacentEdges[5] = wedge[WINDEX((i-1+POINTS)%POINTS,(j-1+POINTS)%POINTS,k,4)];
  paq->adjacentEdges[6] = wedge[WINDEX((i-1+POINTS)%POINTS,j,(k+1)%POINTS,4)];
  paq->adjacentEdges[7] = wedge[WINDEX((i-1+POINTS)%POINTS,j,(k-1+POINTS)%POINTS,4)];
  paq->adjacentEdges[8] = wedge[WINDEX(i,(j+1)%POINTS,(k+1)%POINTS,4)];
  paq->adjacentEdges[9] = wedge[WINDEX(i,(j+1)%POINTS,(k-1+POINTS)%POINTS,4)];
  paq->adjacentEdges[10] = wedge[WINDEX(i,(j-1+POINTS)%POINTS,(k+1)%POINTS,4)];
  paq->adjacentEdges[11] = wedge[WINDEX(i,(j-1+POINTS)%POINTS,(k-1+POINTS)%POINTS,4)];

  // calculate quantities used by evolution functions
  calculatequantities(paq);

  // [EVOLVE WEDGE BASE TO PEAK]
  // energy density
   wedge[INDEX(3,j,k,0)] = grid[INDEX(i,j,k,0)] + dt*energy_evfn(paq);
  // fluid
   wedge[INDEX(3,j,k,1)] = grid[INDEX(i,j,k,1)] + dt*fluid_evfn(paq, 1);
   wedge[INDEX(3,j,k,2)] = grid[INDEX(i,j,k,2)] + dt*fluid_evfn(paq, 2);
   wedge[INDEX(3,j,k,3)] = grid[INDEX(i,j,k,3)] + dt*fluid_evfn(paq, 3);
  // field
   wedge[INDEX(3,j,k,4)] = grid[INDEX(i,j,k,4)] + dt*field_evfn(paq);
  // field derivative
   wedge[INDEX(3,j,k,5)] = grid[INDEX(i,j,k,5)] + dt*ddtfield_evfn(paq);

}


void calculatequantities(PointData *paq)
{
  int u, n;

  // [COMMON QUANTITIES]
  paq->u2 = magu2(paq);
  paq->ut2 = 1.0 + paq->u2;
  paq->ut = sqrt(paq->ut2);
  paq->relw = (1.0 - W_EOSm1 * paq->u2);

  // [GRADIENTS]
  for(u=0; u<5; u++) {
    for(n=1; n<=3; n++)
      paq->gradients[n][u] = derivative(paq, n, u);
  }
  paq->lap = lapl(paq);

  // [SOURCES] little j is source, big J is some wonko function of source
  jsource(paq);
  paq->srcsum = sumvv(paq->fields, paq->ji);
  Jsource(paq); // depends on having paq->srcsum calculated correctly

  // [DERIVED QUANTITIES]
  paq->uudu = sumvtv(paq->fields, paq->gradients, paq->fields);
  paq->trgrad = sp_tr(paq->gradients);
  paq->ude = sumvt(paq->fields, paq->gradients, 1, 0);

}
