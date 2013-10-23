#ifndef MATH_H
#define MATH_H

#include "defines.h"

/* PROTOTYPES */
/* functions to calculate commonly used quantities */
static inline simType Ut(PointData *paq);
static inline simType Ut2(PointData *paq);
static inline simType magu2(PointData *paq);

static inline int sgn(simType value);

/* sums, derivatives, such. */
static inline simType sumvv(simType v1[4], simType v2[4]);
static inline simType sumvt(simType v1[4], simType t1[4][DOF], int rc, int s);
static inline simType sumvtv(simType v1[4], simType t1[4][DOF], simType v2[4]);
static inline simType sp_tr(simType t1[4][DOF]);

static inline simType derivative(PointData *paq, int ddim, int dim);
static inline simType lapl(PointData *paq);

/* potential function */
static inline simType V(simType phi);
static inline simType dV(simType phi);

/* Bubble profile for nucleation */
simType tanhbubble(int i, int j, int k, simType xcent, simType ycent, simType zcent);

/* Data convolution - generic smoothing to reduce Gibbs oscillations */
void convolve(simType *data, simType *temp, simType coeff);

/* Interpolation function? */
// static gsl_spline getfn(IOData filedata, gsl_interp_accel acc);

/* 
 * U^2 - commonly used quantity, stored so we only need to use it once.
 */
static inline simType magu2(PointData *paq)
{
  return paq->fields[1]*paq->fields[1]
      + paq->fields[2]*paq->fields[2]
      + paq->fields[3]*paq->fields[3];
}


/* 
 * (U^t)^2 - commonly used quantity, stored so we only need to use it once.
 */
static inline simType Ut2(PointData *paq)
{
  return 1 + magu2(paq);
}


/* 
 * U^t - commonly used quantity, stored so we only need to use it once.
 */
static inline simType Ut(PointData *paq)
{
  return sqrt(Ut2(paq));
}


/* 
 * Taking derivatives.  This uses pre-stored values to avoid multiple calls to the wedge.
 * Working on the wedge base here, so derivative can only be taken at points 'inside' the base.
 */
static inline simType derivative(PointData *paq, int ddim /* direction of derivative */, int dim)
{
  switch(ddim)
  {
    case 1:
      return (paq->adjacentFields[0][dim] - paq->adjacentFields[3][dim])/2/dx;
    case 2:
      return (paq->adjacentFields[1][dim] - paq->adjacentFields[4][dim])/2/dx;
    case 3:
      return (paq->adjacentFields[2][dim] - paq->adjacentFields[5][dim])/2/dx;
  }

  // Try out alternative differences - prefferentially take a handed derivative.  This
  // is to help prevent instabilities from forming.

  // simType nvel = 0.0;
  // if(ddim <= 3) {
  //   simType mf = 10.0*paq->fields[ddim];
  //   nvel = mf / sqrt(1.0 + mf*mf);
  // } else {
  //   if(paq->fields[5] > 0) {
  //     nvel = 1.0;
  //   } else if(paq->fields[5] < 0) {
  //     nvel = -1.0;
  //   }
  // }

  // switch(ddim)
  // {
  //   case 1:
  //     return (
  //             (1.0 - nvel)*(paq->adjacentFields[0][dim] - paq->fields[dim])
  //             + (1.0 + nvel)*(paq->fields[dim] - paq->adjacentFields[3][dim])
  //           )/2/dx;
  //   case 2:
  //     return (
  //             (1.0 - nvel)*(paq->adjacentFields[1][dim] - paq->fields[dim])
  //             + (1.0 + nvel)*(paq->fields[dim] - paq->adjacentFields[4][dim])
  //           )/2/dx;
  //   case 3:
  //     return (
  //             (1.0 - nvel)*(paq->adjacentFields[2][dim] - paq->fields[dim])
  //             + (1.0 + nvel)*(paq->fields[dim] - paq->adjacentFields[5][dim])
  //           )/2/dx;
  // }

  /* XXX */
  return 0;
}


/* 
 * Taking second derivatives - only of the field.  Again, assumes toroidial boundary conditions in each direction.
 * Bit higher order scheme here than just taking 2nd derivatives.
 */
static inline simType lapl(PointData *paq)
{
  return (
    (
      // Edge-differences
      paq->adjacentEdges[0] + paq->adjacentEdges[1] + paq->adjacentEdges[2] + paq->adjacentEdges[3]
      + paq->adjacentEdges[4] + paq->adjacentEdges[5] + paq->adjacentEdges[6] + paq->adjacentEdges[7]
      + paq->adjacentEdges[8] + paq->adjacentEdges[9] + paq->adjacentEdges[10] + paq->adjacentEdges[11]
    )
    + 2.0*(
      // Center-differences
      paq->adjacentFields[0][4] + paq->adjacentFields[1][4] + paq->adjacentFields[2][4]
      + paq->adjacentFields[3][4] + paq->adjacentFields[4][4] + paq->adjacentFields[5][4]
    )
    - 24.0*(
      paq->fields[4]
    )
  )/6.0/dx/dx;
}



/*
 * Sum function - spatial sum of two 4-vectors.
 */
static inline simType sumvv(simType v1[4], simType v2[4])
{
  return v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}


/*
 * Overload to give spatial sum of 4-vector and a component of rank-2 tensor.
 */
static inline simType sumvt(simType v1[4], simType t1[4][DOF], int rc /* Sum with row (1) or column (2)? */, int s /* row/col to sum over */)
{
  switch(rc)
  {
    case 1:
      return v1[1] * t1[1][s] + v1[2] * t1[2][s] + v1[3] * t1[3][s];
    case 2:
      return v1[1] * t1[s][1] + v1[2] * t1[s][2] + v1[3] * t1[s][3];
    default:
      fprintf(stderr, "Error: sum/row not specified correctly!\n");
      return 0;
  }
}


/*
 * Overload to give full spatial inner product of 2 vectors and rank-2 tensor.
 * First vector is summed with first tensor index, second with second.
 */
static inline simType sumvtv(simType v1[4], simType t1[4][DOF], simType v2[4])
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
static inline simType sp_tr(simType t1[4][DOF])
{
  return t1[1][1] + t1[2][2] + t1[3][3];
}


/* 
 * Symmetry breaking potential - \phi^4 with linear perturbation.
 */
static inline simType dV(simType phi)
{
  // return 0;
  return phi + phi*(3.0/2.0*phi + getALPHA()/2.0*phi*phi);
}


/* 
 * Symmetry breaking potential - \phi^4 with linear perturbation.
 */
static inline simType V(simType phi)
{
  // return 0;
  return phi*(phi/2.0 + phi*phi/2.0 + getALPHA()/8.0*phi*phi*phi) + V_OFFSET;
}


/*
 * Function to get sign of number
 */
static inline int sgn(simType value)
{
  return (value > 0.0) - (value < 0.0);
}

#endif