
#include "defines.h"

/* PROTOTYPES */
/* functions to calculate commonly used quantities */
inline simType Ut(PointData *paq);
inline simType Ut2(PointData *paq);
inline simType magu2(PointData *paq);

/* sums, derivatives, such. */
inline simType sumvv(simType v1[4], simType v2[4]);
inline simType sumvt(simType v1[4], simType t1[4][DOF], int rc, int s);
inline simType sumvtv(simType v1[4], simType t1[4][DOF], simType v2[4]);
inline simType sp_tr(simType t1[4][DOF]);
inline simType derivative(simType *data, int ddim, int dim, int i, int j, int k);
inline simType derivative2(simType *data, int ddim, int dim, int i, int j, int k);

/* potential function */
inline simType dV(simType phi);


/* 
 * U^2 - commonly used quantity, stored so we only need to use it once.
 */
inline simType magu2(PointData *paq)
{
  return paq->fields[1]*paq->fields[1]
      + paq->fields[2]*paq->fields[2]
      + paq->fields[3]*paq->fields[3];
}



/* 
 * (U^t)^2 - commonly used quantity, stored so we only need to use it once.
 */
inline simType Ut2(PointData *paq)
{
  return 1 + magu2(paq);
}

/* 
 * U^t - commonly used quantity, stored so we only need to use it once.
 */
inline simType Ut(PointData *paq)
{
  return sqrt(Ut2(paq));
}



/* 
 * Taking derivatives.  Assumes toroidial boundary conditions in each direction.
 */
inline simType derivative(simType *data, int ddim /* direction of derivative */, int dim, int i, int j, int k)
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
inline simType derivative2(simType *data, int ddim /* direction of derivative */, int dim, int i, int j, int k)
{
  switch(ddim)
  {
    case 1:
      return (data[INDEX((i+1)%POINTS,j,k,dim)]
        + data[INDEX((i+POINTS-1)%POINTS,j,k,dim)]
        - 2*data[INDEX(i,j,k,dim)])/dx/dx;

    case 2:
      return (data[INDEX(i,(j+1)%POINTS,k,dim)]
        + data[INDEX(i,(j+POINTS-1)%POINTS,k,dim)]
        - 2*data[INDEX(i,j,k,dim)])/dx/dx;

    case 3:
      return (data[INDEX(i,j,(k+1)%POINTS,dim)]
        + data[INDEX(i,j,(k+POINTS-1)%POINTS,dim)]
        - 2*data[INDEX(i,j,k,dim)])/dx/dx;
  }

  /* XXX */
  return 0;
}


/*
 * Sum function - spatial sum of two 4-vectors.
 */
inline simType sumvv(simType v1[4], simType v2[4])
{
  return v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}


/*
 * Overload to give spatial sum of 4-vector and a component of rank-2 tensor.
 */
inline simType sumvt(simType v1[4], simType t1[4][DOF], int rc /* Sum with row (1) or column (2)? */, int s /* row/col to sum over */)
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
inline simType sumvtv(simType v1[4], simType t1[4][DOF], simType v2[4])
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
inline simType sp_tr(simType t1[4][DOF])
{
  return t1[1][1] + t1[2][2] + t1[3][3];
}


/* 
 * Symmetry breaking potential - \phi^4 with linear perturbation.
 */
inline simType dV(simType phi)
{
  // return 0;
  return LAMBDA/2*(phi*phi - ETA*ETA)*phi + EPSILON*LAMBDA*ETA*ETA*ETA;
}

