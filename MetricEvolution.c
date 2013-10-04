/* RNSE includes */
#include "defines.h"


/*
 * Calculate the actual gravitational wave spectrum from l
 */
void get_gws(PointData *paq, simType **l, int i, int j, int k)
{

}


/*
 * For evolving h_ij, we work in fourier space, as this has proven to be more stable in past work.
 */
void h_evolve(simType **hij, simType **lij, simType **fSTTij)
{
  int i, j, k, a;
  simType pz, px, py;
  simType dpx, dpz, dpy;

  // evolve h_ij (the fourier transform of)
    // solving (in fourier space) d_t l = k^2 h + S
    //                            d_t h = h


  // only evolve the transverse/traceless part:
  // #pragma omp parallel for default(shared) private(i, j, k) num_threads(threads)
  for(int i=0; i<POINTS; i++)
  {
    px = (simType) (i <= POINTS/2 ? i : i - POINTS);
    for(int j=0; j<POINTS; j++)
    {
      py = (simType) (j <= POINTS/2 ? j : j - POINTS);
      for(int k=0; k<N/2+1; k++)
      {
        pz = (simType) k;
        p2 = px*px + py*py + pz*pz; // p^2 (momentum)
        simType phat[3] = { px/sqrt(p2), py/sqrt(p2), pz/sqrt(p2) };

        // pre-calculate some stuff:          
          // Trace (S^a_a)
          simType trs_RE = C_RE(fSTTij[0][SINDEX(i,j,k)]) + C_RE(fSTTij[3][SINDEX(i,j,k)]) + C_RE(fSTTij[5][SINDEX(i,j,k)]);
          simType trs_IM = C_IM(fSTTij[0][SINDEX(i,j,k)]) + C_IM(fSTTij[3][SINDEX(i,j,k)]) + C_IM(fSTTij[5][SINDEX(i,j,k)]);
          // One contraction (k^m S_mi)
          simType k_mS_mi_RE[3] = { 
            phat[1]*C_RE(fSTTij[0][SINDEX(i,j,k)]) + phat[2]*C_RE(fSTTij[1][SINDEX(i,j,k)]) + phat[3]*C_RE(fSTTij[2][SINDEX(i,j,k)]),
            phat[1]*C_RE(fSTTij[1][SINDEX(i,j,k)]) + phat[2]*C_RE(fSTTij[3][SINDEX(i,j,k)]) + phat[3]*C_RE(fSTTij[4][SINDEX(i,j,k)]),
            phat[1]*C_RE(fSTTij[2][SINDEX(i,j,k)]) + phat[2]*C_RE(fSTTij[4][SINDEX(i,j,k)]) + phat[3]*C_RE(fSTTij[5][SINDEX(i,j,k)])
          };
          simType k_mS_mi_RE[3] = { 
            phat[1]*C_IM(fSTTij[0][SINDEX(i,j,k)]) + phat[2]*C_IM(fSTTij[1][SINDEX(i,j,k)]) + phat[3]*C_IM(fSTTij[2][SINDEX(i,j,k)]),
            phat[1]*C_IM(fSTTij[1][SINDEX(i,j,k)]) + phat[2]*C_IM(fSTTij[3][SINDEX(i,j,k)]) + phat[3]*C_IM(fSTTij[4][SINDEX(i,j,k)]),
            phat[1]*C_IM(fSTTij[2][SINDEX(i,j,k)]) + phat[2]*C_IM(fSTTij[4][SINDEX(i,j,k)]) + phat[3]*C_IM(fSTTij[5][SINDEX(i,j,k)])
          };
          // Two contractions (k^a k^b S_ab)
          simType kkS_RE = phat[1]*k_mS_mi_RE[1] + phat[2]*k_mS_mi_RE[2] + phat[3]*k_mS_mi_RE[3];
          simType kkS_IM = phat[1]*k_mS_mi_IM[1] + phat[2]*k_mS_mi_IM[2] + phat[3]*k_mS_mi_IM[3];

        for(a=1; a<=3; a++) {
          for(b=a; b<=3; b++) {
            // (7-a)*a/2-4+b formula maps indexes of h to sequential indices
            C_RE(fSTTij[(7-a)*a/2-4+b][SINDEX(i,j,k)]) +=
              (
                0.5*(phat[a]*phat[b] - (a==b))*trs_RE
                + 0.5*(phat[a]*phat[b] + (a==b))*kkS_RE
                - (phat[a]*k_mS_mi_RE[b] + phat[b]*k_mS_mi_RE[a])
              );
            C_IM(fSTTij[(7-a)*a/2-4+b][SINDEX(i,j,k)]) +=
              (
                0.5*(phat[a]*phat[b] - (a==b))*trs_RE
                + 0.5*(phat[a]*phat[b] + (a==b))*kkS_RE
                - (phat[a]*k_mS_mi_RE[b] + phat[b]*k_mS_mi_RE[a])
              );
          }
        }

      }
    }
  }

  for(a=0; a<12; a++) {
    // #pragma omp parallel for default(shared) private(i, j, k) num_threads(threads)
    for(int i=0; i<POINTS; i++)
    {
      px = (simType) (i <= POINTS/2 ? i : i - POINTS);
      for(int j=0; j<POINTS; j++)
      {
        py = (simType) (j <= POINTS/2 ? j : j - POINTS);
        for(int k=0; k<N/2+1; k++)
        {
          pz = (simType) k;
          p = sqrt(px*px + py*py + pz*pz); // p (momentum)
          simType h = hij[a][SINDEX(i,j,k)];
          simType l = lij[a][SINDEX(i,j,k)];

          simType S;
          if(a >= 6) {
            S = C_IM(fSTTij[a-6][SINDEX(i,j,k)]);
          } else {
            S = C_RE(fSTTij[a][SINDEX(i,j,k)]);
          }

          // This is RK4. Promise.
          hij[a][SINDEX(i,j,k)] += dt*(l*(6.0 + p*p*dt*dt) + dt*(h*k*k + S)*(12.0 + dt*dt*p*p)/4.0)/6.0;
          lij[a][SINDEX(i,j,k)] += dt*(6.0*S + 3.0*p*p*l*dt + p*p*p*p*l*dt*dt*dt/4.0 + h*p*p*(6.0 + p*p*dt*dt))/6.0;
        }
      }
    }
  }

  return;
}


/*
 * FFT the STT
 */
void fft_stt(simType **STTij, fftw_complex **fSTTij)
{
  // transform STT components into fourier space
  int a;
  fftw_plan p;
  p = fftw_plan_dft_r2c_3d(POINTS, POINTS, POINTS,
                            STTij[0], fSTTij[0],
                            FFTW_ESTIMATE);

  // fourier transform stress-energy tensor
  for(a=1; a<6; a++) {
    fftw_execute_dft_r2c(plan_fft, STTij[a], fSTTij[a]);
  }
}


/*
 * Calculate S at a point (i,j,k).
 * This does not remove the transverse/traceless part.
 */
void set_stt(PointData *paq, simType **STTij, int i, int j, int k)
{
  // S_ij = T_ij - 1/3 \delta_ij * T_k^k.
  //      = (e+p)U^iU^j + d^ifd^jf - 1/3(\delta_ij)*( trace of previous terms )
  // The above is traceless, but not yet transverse.

  int a, b;
  for(a=1; a<=3; a++) {
    for(b=a; b<=3; b++) {
      // (7-a)*a/2-4+b formula maps indexes of h to sequential indices
      STTij[(7-a)*a/2-4+b][SINDEX(i,j,k)] =
          pow(paq->gradients[a][5],2)
            + W_EOSp1 * exp(paq->fields[0]) * paq->fields[a] * paq->fields[b];
    }
  }

  // project into transverse plane, and subtract trace - this is done during actual evolution.

  return;
}
