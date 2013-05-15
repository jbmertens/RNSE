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
static inline simType derivative(simType *data, int ddim, int dim, int i, int j, int k);
static inline simType derivative2(simType *data, int ddim, int dim, int i, int j, int k);
static inline simType lapl(simType *data, int dof, int i, int j, int k);

/* potential function */
static inline simType dV(simType phi);

/* Data convolution - generic smoothing to reduce Gibbs oscillations */
static inline void convolve(simType *data, simType *temp, simType coeff);

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
 * Taking derivatives.  Assumes toroidial boundary conditions in each direction.
 */
static inline simType derivative(simType *data, int ddim /* direction of derivative */, int dim, int i, int j, int k)
{
  /* taking modulo here for each point */
  switch(ddim)
  {
    case 1:
      return (data[INDEX((i+1)%POINTS,j,k,dim)]
        - data[INDEX((i+POINTS-1)%POINTS,j,k,dim)])/2/dx;

    case 2:
      return (data[INDEX(i,(j+1)%POINTS,k,dim)]
        - data[INDEX(i,(j+POINTS-1)%POINTS,k,dim)])/2/dx;

    case 3:
      return (data[INDEX(i,j,(k+1)%POINTS,dim)]
        - data[INDEX(i,j,(k+POINTS-1)%POINTS,dim)])/2/dx;

  }

  /* XXX */
  return 0;
}


/* 
 * Taking second derivatives.  Again, assumes toroidial boundary conditions in each direction.
 */
static inline simType derivative2(simType *data, int ddim /* direction of derivative */, int dim, int i, int j, int k)
{
  switch(ddim)
  {
    case 1:
      return (
        data[INDEX((i+1)%POINTS,j,k,dim)]
        + data[INDEX((i+POINTS-1)%POINTS,j,k,dim)]
        - 2*data[INDEX(i,j,k,dim)]
        )/dx/dx;

    case 2:
      return (
        data[INDEX(i,(j+1)%POINTS,k,dim)]
        + data[INDEX(i,(j+POINTS-1)%POINTS,k,dim)]
        - 2*data[INDEX(i,j,k,dim)]
        )/dx/dx;

    case 3:
      return (
        data[INDEX(i,j,(k+1)%POINTS,dim)]
        + data[INDEX(i,j,(k+POINTS-1)%POINTS,dim)]
        - 2*data[INDEX(i,j,k,dim)]
        )/dx/dx;
  }

  /* XXX */
  return 0;
}

/* 
 * Taking second derivatives.  Again, assumes toroidial boundary conditions in each direction.
 * Bit higher order scheme here than just taking 2nd derivatives.
 */
static inline simType lapl(simType *data, int dof, int i, int j, int k)
{
  return (
    (
      data[INDEX((i+1)%POINTS,(j+1)%POINTS,k,dof)] + data[INDEX((i+1)%POINTS,(j-1+POINTS)%POINTS,k,dof)]
      + data[INDEX((i+1)%POINTS,j,(k+1)%POINTS,dof)] + data[INDEX((i+1)%POINTS,j,(k-1+POINTS)%POINTS,dof)]
      + data[INDEX((i-1+POINTS)%POINTS,(j+1)%POINTS,k,dof)] + data[INDEX((i-1+POINTS)%POINTS,(j-1+POINTS)%POINTS,k,dof)]
      + data[INDEX((i-1+POINTS)%POINTS,j,(k+1)%POINTS,dof)] + data[INDEX((i-1+POINTS)%POINTS,j,(k-1+POINTS)%POINTS,dof)]
      + data[INDEX(i,(j+1)%POINTS,(k+1)%POINTS,dof)] + data[INDEX(i,(j+1)%POINTS,(k-1+POINTS)%POINTS,dof)]
      + data[INDEX(i,(j-1+POINTS)%POINTS,(k+1)%POINTS,dof)] + data[INDEX(i,(j-1+POINTS)%POINTS,(k-1+POINTS)%POINTS,dof)]
    )
    + 2.0*(
      data[INDEX((i+1)%POINTS,j,k,dof)] + data[INDEX(i,(j+1)%POINTS,k,dof)] + data[INDEX(i,j,(k+1)%POINTS,dof)]
      + data[INDEX((i-1+POINTS)%POINTS,j,k,dof)] + data[INDEX(i,(j-1+POINTS)%POINTS,k,dof)] + data[INDEX(i,j,(k-1+POINTS)%POINTS,dof)]
    )
    - 24.0*(
      data[INDEX(i,j,k,dof)]
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
  return LAMBDA/2*(phi*phi - ETA*ETA)*phi + EPSILON*LAMBDA*ETA*ETA*ETA;
}


static inline void convolve(simType *data, simType *temp, simType coeff)
{
  int i, j, k, u;
  int x, y, z, n;

  LOOP4(i,j,k,u)
  {
    x=0; y=0; z=0;

    // check for "W" or "M" patterns to smooth
    for(n=0; n<4; n++)
    {
      // add up greater or less than's, convolved with "w" of signs.  x/y/z = +/-4 for W/M shapes.
      x += (2*(n%2)-1)
          * sgn(round(1000000*(data[INDEX((i+n-2+POINTS)%POINTS,j,k,u)] - data[INDEX((i+n-1+POINTS)%POINTS,j,k,u)])));
      y += (2*(n%2)-1)
          * sgn(round(1000000*(data[INDEX(i,(j+n-2+POINTS)%POINTS,k,u)] - data[INDEX(i,(j+n-1+POINTS)%POINTS,k,u)])));
      z += (2*(n%2)-1)
          * sgn(round(1000000*(data[INDEX(i,j,(k+n-2+POINTS)%POINTS,u)] - data[INDEX(i,j,(k+n-1+POINTS)%POINTS,u)])));
    }

    temp[INDEX(i,j,k,u)] = data[INDEX(i,j,k,u)];
    // if they exist, convolve with a smoothing kernel
    if(abs(x) >= 4) {
      temp[INDEX(i,j,k,u)] += coeff*(
        data[INDEX((i+1)%POINTS,j,k,u)]
        + data[INDEX((i-1+POINTS)%POINTS,j,k,u)]
        - 2.0*data[INDEX(i,j,k,u)]
      );
    }
    if(abs(y) >= 4) {
      temp[INDEX(i,j,k,u)] += coeff*(
        data[INDEX(i,(j+1)%POINTS,k,u)]
        + data[INDEX(i,(j-1+POINTS)%POINTS,k,u)]
        - 2.0*data[INDEX(i,j,k,u)]
      );
    }
    if(abs(z) >= 4) {
      temp[INDEX(i,j,k,u)] += coeff*(
        data[INDEX(i,j,(k+1)%POINTS,u)]
        + data[INDEX(i,j,(k-1+POINTS)%POINTS,u)]
        - 2.0*data[INDEX(i,j,k,u)]
      );
    }
  }

  // put back smoothed data
  LOOP4(i,j,k,u)
  {
    data[INDEX(i,j,k,u)] = temp[INDEX(i,j,k,u)];
  }
}

static inline int sgn(simType value)
{
  return (value > 0.0) - (value < 0.0);
}

#endif
