#ifndef DEFINES_H
#define DEFINES_H

/* INCLUDES */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <math.h>
#include <sys/stat.h> /* mkdir() */
#include <time.h>
#include <omp.h>
#include <fftw3.h>
#include <getopt.h>

/* CONSTANTS */
/* potential parameters: must change the initial field configuration if these
   are changed! */
#define LAMBDA      0.01    /* potential height */
#define ETA         1.0     /* potential min position */
#define EPSILON     0.025   /* potential asymmetry */
#define R0          (1 / EPSILON / ETA / sqrt(LAMBDA))

#define W_EOS       (1.0 / 3.0)   /* EOS parameter */
#define W_EOSm1     (W_EOS - 1.0)
#define W_EOSp1     (W_EOS + 1.0)

/*1 (ln of) fluid density, 3 fluid, 1 d/dt scalar field, 1 scalar field */
#define DOF 6  

/* resolution parameters */
#define SIZE    (6*R0)                         /* physical size in space */
#define POINTS  ((long long) 512)              /* number of points on
                                                  lattice (each axis) */
#define dx      ( 1.0*SIZE / (1.0*POINTS) )
#define dt      (dx/20.0)

/* storage parameters */
#define RK_STEPS 0                            /* Usually RK Method Order - 1        */
#define RANK 4                                /* dimension of fields array          */
#define GRID_STORAGE (POINTS*POINTS*POINTS)   /* space requirement for just grid    */
#define STORAGE (GRID_STORAGE*DOF)            /* space requirement for fields/fluid */

/* simulation sampling information */
#define MAX_STEPS           1             /* Maximum # of steps to run */
#define STEPS_TO_SAMPLE     1             /* # of steps to record, undersampled */
#define STEPS_TO_DUMP       0             /* # of steps to give a full dump of and take DHT */
#define POINTS_TO_SAMPLE    50            /* # of points along (x-)axis to
                                               sample */

#define STOP_CELL           10            /* Check for negative field value (eg, bubble wall)
                                                this many voxels away from the boundary.
                                                If the field value in this cell becomes negative,
                                                the simulation should stop running.
                                                Use negative value for no stop. */

/* some defaults. */
#define DEFAULT_DATA_DIR      "data"
#define DEFAULT_DATA_NAME     "data"

/* array element access macro */
#define INDEX(i,j,k,l) (DOF*POINTS*POINTS*(i) + DOF*POINTS*(j) + DOF*(k) + (l))

/* TYPEDEFS */
/* Precision/format we'd like to use for this simulation: */
typedef double simType;

/* STRUCTS */
/* data structure for storing calculated quantities at a point */
typedef struct {
  simType fields[6];
  simType gradients[4][DOF];
  simType derivs2[4];
  
  /* calculated quantities */
  simType ut;
  simType ut2;
  simType u2;
  simType relw;
  
  /* Source terms */
  simType ji[4];
  simType Ji[4];

  /* sums */
  simType uudu;
  simType srcsum;
  simType trgrad;
  simType udu;

  /* S^{TT}_{ij}
    Gauge choices and coordinate choices allow us to only calculate two components.
    Of the original 9 components of h_ij, we can constrain a given number of components:
      Symmetric => 3
      Transverse => 3
      Traceless => 1
    Giving 9-3-3-1 = 2 independent degrees of freedom.
    However, degeneracies may exist around certain values.  Evolving 3 independent values
    quantities helps break this degeneracy, and also allows for consistency checks.

    Here, we evolve h_11, h_13, and h_33.  See writeup for specifics on consistency checking
    and formulation behind the evolution.
   */
  simType STT[3];

} PointData;

/* file information storage */
typedef struct {
  char *data_dir;
  char *data_name;
  char *read_data_name;
  int fwrites;
  int datasize;
} IOData;


/* Project functionality */
#include "RNSFluid.h"
#include "MetricEvolution.h"
#include "io.h"
#include "fft_util.h"
#include "math_util.h"


#endif
