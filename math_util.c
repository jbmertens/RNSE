#include "defines.h"


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


/* Interpolation function? 
static gsl_spline getfn(IOData filedata, gsl_interp_accel acc)
{

}*/

