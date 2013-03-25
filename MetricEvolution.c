/* RNSE includes */
#include "defines.h"

/*
 * For evolving h_ij, we work in fourier space.  Gauge choices reduce 6 degrees of freedom to 2,
 * so we in principle only need to evolve 2 degrees of freedom.
 * We also evolve in fourier space, since this has proven to be more stable in past work.
 */
void evolve_hij(simType *metric_pert)
{
  // transform metric into fourier space
}

/*
 * Return S^{TT}_{ij} at a point.
 */
void set_stt(PointData *paq, int i, int j, int k)
{
  // S_ij = T_ij - 1/3 \delta_ij * T_k^k.
  //      = (e+p)U^iU^j + d^ifd^jf - 1/3(\delta_ij)*( trace of previous terms )

  // store data in preallocated space.
  // first two terms:
  int a, b;
  for(a=1; a<=3; a++) {
    for(b=a; b<=3; b++) {
      paq->STT[(7-a)*a/2-4+b] = pow(paq->gradients[a][5],2)
                                  + W_EOSp1 * exp(paq->fields[0]) * paq->fields[a] * paq->fields[b];
    }
  }

  // calculate trace
  simType tr = 0;
  for(a=1; a<=3; a++) {
    tr += paq->STT[(9-a)*a/2-4];
  }

  // subtract trace term from diagonal elements
  for(a=1; a<3; a++) {
    paq->STT[(9-a)*a/2-4] -= tr/3.0
  }

  return;
}


