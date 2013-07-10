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
#include <zlib.h>

/* CONSTANTS */
/* potential parameters: must change the initial field configuration if these
   are changed! */
#define R0         (1.0/(1.0-getALPHA()))

#define W_EOS       (1.0 / 3.0)   /* EOS parameter */
#define W_EOSm1     (W_EOS - 1.0)
#define W_EOSp1     (W_EOS + 1.0)

/*1 (ln of) fluid density, 3 fluid, 1 d/dt scalar field, 1 scalar field */
#define DOF ((long long) 6)  

/* resolution parameters */
#define SIZE    (10*R0)                         /* physical size in space */
#define POINTS  ((long long) 128)              /* number of points on
                                                  lattice (each axis) */
#define dx      (1.0 * SIZE / 1.0 / POINTS)
#define dt      (dx/25.0)

/* storage parameters */
#define METHOD_ORDER 2                        /* Order of method - assumes diagonal Butcher tableau */
#define RANK 4                                /* dimension of fields array                          */
#define GRID_STORAGE (POINTS*POINTS*POINTS)   /* space requirement for just grid                    */
#define STORAGE (GRID_STORAGE*DOF)            /* space requirement for fields/fluid                 */
#define AREA_STORAGE (POINTS*POINTS*DOF)      /* Storage needed for an array in the "wedge"         */

/* simulation sampling information */
#define MAX_STEPS           10000             /* Maximum # of steps to run */
#define STEPS_TO_SAMPLE     0             /* # of steps to record, undersampled */
#define STEPS_TO_DUMP       0             /* # of steps to give a full dump of and take DHT */
#define POINTS_TO_SAMPLE    50            /* # of points along (x-)axis to
                                               sample */

/* special features */
#define STOP_CELL           10  /* Check for negative field value (eg, bubble wall)
                                   this many voxels away from the boundary.
                                   If the field value in this cell becomes negative,
                                   the simulation should stop running.
                                   Use negative value for no stop. */
#define STOP_MAX            ((-3.0 + sqrt(9.0 - 8.0*getALPHA()))/2/getALPHA())
#define DUMP_STRIP 1            /* bool - full dump of a strip along the center of
                                   the simulation at each step? */

/* some defaults. */
#define DEFAULT_DATA_DIR      "data"
#define DEFAULT_DATA_NAME     "data"

/* array element access macros */
#define INDEX(i,j,k,l) (DOF*POINTS*POINTS*((i)%POINTS) + DOF*POINTS*(j) + DOF*(k) + (l))
#define WINDEX(i,j,k,l) (DOF*POINTS*POINTS*(((i)%POINTS)%3) + DOF*POINTS*(j) + DOF*(k) + (l))
/* Common for loop structure */
#define LOOP2(j,k) for(j=0; j<POINTS; j++) \
                   for(k=0; k<POINTS; k++)
#define LOOP3(i,j,k) for(i=0; i<POINTS; i++) LOOP2(j,k)
#define LOOP4(i,j,k,u) LOOP3(i,j,k) for(u=0; u<DOF; u++)

/* TYPEDEFS */
/* Precision/format we'd like to use for this simulation: */
typedef double simType;

/* STRUCTS */
/* data structure for storing calculated quantities at a point */
typedef struct {
  simType fields[6];
  simType adjacentFields[6][DOF];
  simType adjacentEdges[12]; // only needed for laplacian at the moment

  simType gradients[4][DOF];
  simType derivs2[4];
  simType lap;

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
  simType ude;

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
