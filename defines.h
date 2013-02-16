#ifndef DEFINES_H
#define DEFINES_H

/* potential parameters: must change the initial field configuration if these
   are changed! */
#define LAMBDA      0.01    /* potential height */
#define ETA         1.0     /* potential min position */
#define EPSILON     0.025   /* potential asymmetry */
#define R0          ( 1/EPSILON/ETA/sqrt(LAMBDA) )

#define W_EOS       0.333   /* EOS parameter */
#define W_EOSm1     (W_EOS - 1.0)
#define W_EOSp1     (W_EOS + 1.0)
#define XI          (-0.01)   /* fluid coupling */

/*1 (log of) fluid density, 3 fluid, 1 d/dt scalar field, 1 scalar field */
#define DOF 6  

/* resolution parameters */
#define SIZE    (4*R0)                         /* physical size in space */
#define POINTS  100                             /* number of points on
                                                  lattice (each axis) */
#define dx      ( 1.0*SIZE / (1.0*POINTS) )
#define dt      (dx/20.0)

/* storage parameters */
#define RK_STEPS 0                         /* Usually RK Method Order - 1  */
#define RANK 4                             /* dimension of fields array    */
#define STORAGE POINTS*POINTS*POINTS*DOF   /* space requirement            */

#define STEPS               80            /* # of steps to run */
#define STEPS_TO_RECORD     60            /* # of steps to record */
#define POINTS_TO_SAMPLE    25             /* # of points along (x-)axis to
                                               sample */

/* array element access macro */
#define INDEX(i,j,k,l) (DOF*POINTS*POINTS*(i) + DOF*POINTS*(j) + DOF*(k) + (l))

#endif
