#ifndef DEFINES_H
#define DEFINES_H

/* INCLUDES */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <zlib.h>
#include <math.h>
#include <sys/stat.h> /* mkdir() */
#include <time.h>
#include <omp.h>
#include <fftw3.h>

/* CONSTANTS */
/* potential parameters: must change the initial field configuration if these
   are changed! */
#define LAMBDA      0.01    /* potential height */
#define ETA         1.0     /* potential min position */
#define EPSILON     0.025   /* potential asymmetry */
#define R0          ( 1/EPSILON/ETA/sqrt(LAMBDA) )

#define W_EOS       0.333   /* EOS parameter */
#define W_EOSm1     (W_EOS - 1.0)
#define W_EOSp1     (W_EOS + 1.0)

/*1 (log of) fluid density, 3 fluid, 1 d/dt scalar field, 1 scalar field */
#define DOF 6  

/* resolution parameters */
#define SIZE    (4*R0)                         /* physical size in space */
#define POINTS  200                            /* number of points on
                                                  lattice (each axis) */
#define dx      ( 1.0*SIZE / (1.0*POINTS) )
#define dt      (dx/20.0)

/* storage parameters */
#define RK_STEPS 0                         /* Usually RK Method Order - 1  */
#define RANK 4                             /* dimension of fields array    */
#define STORAGE (POINTS*POINTS*POINTS*DOF) /* space requirement            */

#define STEPS               500            /* # of steps to run */
#define STEPS_TO_SAMPLE     100            /* # of steps to record, undersampled */
#define STEPS_TO_DUMP       5              /* # of steps to give a full dump of and take DHT */
#define POINTS_TO_SAMPLE    50             /* # of points along (x-)axis to
                                               sample */

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

  /* Constants */
  double xi;
  
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

} PointData;


/* io functionality prototypes */
typedef struct {
  char *data_dir;
  char *data_name;
  int fwrites;
  int datasize;
} IOData;

void writeinfo(IOData filedata);
void dumpstate(simType *fields, IOData filedata);
void readstate(simType *fields, IOData filedata);


/* fft functionality */
void hartleydump(simType *fields, simType *storage, IOData filedata);


/* inline functions for fast math */
#include "math_util.h"


#endif
