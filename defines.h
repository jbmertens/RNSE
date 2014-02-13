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

#define BETA      0.01  /* Ratio of energy density to well height difference */
#define DELTA_V  (pow(sqrt(9.0-8.0*getALPHA())+3.0,2)*(sqrt(9.0-8.0*getALPHA())+3.0-4.0*getALPHA())/64.0/pow(getALPHA(),3))
#define LOG_E     log(DELTA_V/BETA)  /* Log of fluid energy density */
#define V_OFFSET  0.0 /* ~ VEV of true vacuum*/

/*1 (ln of) fluid density, 3 fluid, 1 d/dt scalar field, 1 scalar field */
#define DOF ((long long) 6)  

/* resolution parameters */
#define SIZE    (36*R0)                         /* physical size in space */
#define POINTS  ((long long) 256)               /* number of points on
                                                  lattice (each axis) */
#define dx      (1.0 * SIZE / 1.0 / POINTS)
#define dt      (dx/10.0)
#define gdt_dt  2                              /* GW spectrum evolves only every `gdt_dt` steps. */

/* Metric evolution - should have at least one projection here. */
#define PROJECT_L_TO_LTT 1
#define PROJECT_S_TO_STT 0

/* storage parameters */
#define METHOD_ORDER 2                        /* Order of method - assumes diagonal Butcher tableau */
#define RANK 4                                /* dimension of fields array                          */
#define GRID_STORAGE (POINTS*POINTS*POINTS)   /* space requirement for just grid                    */
#define STORAGE (GRID_STORAGE*DOF)            /* space requirement for fields/fluid                 */
#define AREA_STORAGE (POINTS*POINTS*DOF)      /* Storage needed for an array in the "wedge"         */

/* simulation sampling information */
#define MAX_STEPS           6000       /* Maximum # of steps to run */
#define STEPS_TO_SAMPLE     60         /* # of steps to record, undersampled */
#define STEPS_TO_DUMP       0          /* # of steps to give a full dump of and take DHT */
#define POINTS_TO_SAMPLE    64         /* # of points along (x-)axis to
                                          sample.  This should evenly divide POINTS.  */

/* special features */
#define STOP_CELL           -1  /* Check for negative field value (eg, bubble wall)
                                   this many voxels away from the boundary.
                                   If the field value in this cell becomes negative,
                                   the simulation should stop running.
                                   Use negative value for no stop. */
#define STOP_MAX            ((-3.0 + sqrt(9.0 - 8.0*getALPHA()))/2/getALPHA())
#define DUMP_STRIP 0            /* bool - full dump of a 1-d strip along the center of
                                   the simulation at each step? */

/* some defaults. */
#define DEFAULT_DATA_DIR      "data"
#define DEFAULT_DATA_NAME     "data"

/* array element access macros */
// fields array
#define INDEX(i,j,k,l) (DOF*POINTS*POINTS*((i+POINTS)%POINTS) + DOF*POINTS*((j+POINTS)%POINTS) + DOF*((k+POINTS)%POINTS) + (l))
// wedge array structure
#define WINDEX(i,j,k,l) (DOF*POINTS*POINTS*((i+3)%3) + DOF*POINTS*((j+POINTS)%POINTS) + DOF*((k+POINTS)%POINTS) + (l))
// stress-energy (pre-fft)
#define SINDEX(i,j,k) (POINTS*POINTS*(i) + POINTS*(j) + (k))
// stress-energy (fft-valued)
#define fSINDEX(i,j,k) ((POINTS/2+1)*POINTS*(i) + (POINTS/2+1)*(j) + (k))

/* Common for loop structure */
#define LOOP2(j,k) for(j=0; j<POINTS; j++) \
                   for(k=0; k<POINTS; k++)
#define LOOP3(i,j,k) for(i=0; i<POINTS; i++) LOOP2(j,k)
#define LOOP4(i,j,k,u) LOOP3(i,j,k) for(u=0; u<DOF; u++)

#define C_RE(c) ((c)[0])
#define C_IM(c) ((c)[1])
#define pw2(x) ((x)*(x))

#define M_PI 3.14159265358979323846

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

  simType visc[4];

  /* sums */
  simType uudu;
  simType srcsum;
  simType trgrad;
  simType ude;

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
