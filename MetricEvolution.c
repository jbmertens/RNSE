/* RNSE includes */
#include "defines.h"

/*
 * Funcitonality in this file evolves the fourier trasnforlm of the
 * metric perturbations.  Te use, run (in this order):
 *  set_stt
 *  fft_stt
 *  h_evolve
 *  store_gws
 */

/*
 * Calculate the actual gravitational wave spectrum from l
 */
void store_gws(simType **lij, IOData filedata)
{
  simType array_out[(int)(1.73205*(POINTS/2))+1];
  int numpoints_gw[(int)(1.73205*(POINTS/2))+1]; // Number of points in each momentum bin
  simType p[(int)(1.73205*(POINTS/2))+1];
  simType f2_gw[(int)(1.73205*(POINTS/2))+1]; // Values for each bin: Momentum, |F-k|^2, n_k
  int numbins_gw = (int)(1.73205*(POINTS/2))+1; // Actual number of bins for the number of dimensions

  double pmagnitude_gw; // Total momentum (p) in units of lattice spacing, pmagnitude = Sqrt(px^2+py^2+pz^2). This also gives the bin index since bin spacing is set to equal lattice spacing.
  double fp2_gw;
  int i, j, k, px, py, pz; // px, py, and pz are components of momentum in units of grid spacing

  // Initial magnitude of momentum in each bin
  for(i=0; i<numbins_gw; i++) {
    f2_gw[i] = 0.0;
    numpoints_gw[i] = 0;
  }
    
  // Perform average over all angles here (~integral d\Omega).
  for(i=0; i<POINTS; i++)
  {
    px = (i<=POINTS/2 ? i : i-POINTS);
    for(j=0; j<POINTS; j++)
    {
      py = (j<=POINTS/2 ? j : j-POINTS);
      for(k=1; k<POINTS/2; k++)
      {
        pz = k;
        pmagnitude_gw = sqrt(pw2(px) + pw2(py) + pw2(pz));
        fp2_gw = pw2(lij[0][fSINDEX(i,j,k)]) + pw2(lij[6][fSINDEX(i,j,k)])
          + 2.*pw2(lij[1][fSINDEX(i,j,k)]) + 2.*pw2(lij[7][fSINDEX(i,j,k)])
          + 2.*pw2(lij[2][fSINDEX(i,j,k)]) + 2.*pw2(lij[8][fSINDEX(i,j,k)])
          + pw2(lij[3][fSINDEX(i,j,k)]) + pw2(lij[9][fSINDEX(i,j,k)])
          + 2.*pw2(lij[4][fSINDEX(i,j,k)]) + 2.*pw2(lij[10][fSINDEX(i,j,k)])
          + pw2(lij[5][fSINDEX(i,j,k)]) + pw2(lij[11][fSINDEX(i,j,k)]);
        numpoints_gw[(int)pmagnitude_gw] += 2;
        f2_gw[(int)pmagnitude_gw] += 2.*fp2_gw;
      }

      pz = 0;
      k = 0;
      pmagnitude_gw = sqrt(pw2(px) + pw2(py) + pw2( pz));
      fp2_gw = pw2(lij[0][fSINDEX(i,j,k)]) + pw2(lij[6][fSINDEX(i,j,k)])
        + 2.*pw2(lij[1][fSINDEX(i,j,k)]) + 2.*pw2(lij[7][fSINDEX(i,j,k)])
        + 2.*pw2(lij[2][fSINDEX(i,j,k)]) + 2.*pw2(lij[8][fSINDEX(i,j,k)])
        + pw2(lij[3][fSINDEX(i,j,k)]) + pw2(lij[9][fSINDEX(i,j,k)])
        + 2.*pw2(lij[4][fSINDEX(i,j,k)]) + 2.*pw2(lij[10][fSINDEX(i,j,k)])
        + pw2(lij[5][fSINDEX(i,j,k)]) + pw2(lij[11][fSINDEX(i,j,k)]);  
      numpoints_gw[(int)pmagnitude_gw] += 1;
      f2_gw[(int)pmagnitude_gw] += fp2_gw;
        
      pz = POINTS/2;
      k = POINTS/2;
      pmagnitude_gw = sqrt(pw2(px) + pw2(py) + pw2(pz));
      fp2_gw = pw2(lij[0][fSINDEX(i,j,k)]) + pw2(lij[6][fSINDEX(i,j,k)])
        + 2.*pw2(lij[1][fSINDEX(i,j,k)]) + 2.*pw2(lij[7][fSINDEX(i,j,k)])
        + 2.*pw2(lij[2][fSINDEX(i,j,k)]) + 2.*pw2(lij[8][fSINDEX(i,j,k)])
        + pw2(lij[3][fSINDEX(i,j,k)]) + pw2(lij[9][fSINDEX(i,j,k)])
        + 2.*pw2(lij[4][fSINDEX(i,j,k)]) + 2.*pw2(lij[10][fSINDEX(i,j,k)])
        + pw2(lij[5][fSINDEX(i,j,k)]) + pw2(lij[11][fSINDEX(i,j,k)]);
      numpoints_gw[(int)pmagnitude_gw] += 1;
      f2_gw[(int)pmagnitude_gw] += fp2_gw;
    }
  }
      
  // when also multiplied by 'i^3' in the following loop, this is \Omega_gw * h^2.
  int O_gw_scale = 8.0*pow(10.0,-6)*pow(2.0*M_PI/SIZE,3)*BETA*BETA/DELTA_V/DELTA_V/pow(SIZE,5);
  for(i=0; i<numbins_gw; i++)
  {
    // Converts sums to averages. (numpoints[i] should always be greater than zero.)
    if(numpoints_gw[i] > 0)
    {
      array_out[i] = i*i*i*O_gw_scale*f2_gw[i]/((double) numpoints_gw[i]);
    }
    else
    {
      array_out[i] = 0.;
    }
  }


  // write data
  char *filename, *buffer;
  gzFile *datafile;

  filename = malloc(200 * sizeof(*filename));
  strcpy(filename, filedata.data_dir);
  strcat(filename, filedata.data_name);
  strcat(filename,  ".gwspec.dat.gz");
  buffer = malloc(20 * sizeof(*buffer));

  datafile = (gzFile *)gzopen(filename, "ab");
  if(datafile == Z_NULL) {
    printf("Error opening file: %s\n", filename);
    return;
  }

  for(i=0; i<numbins_gw; i++)
  {
    // field values
    sprintf(buffer, "%g\t", array_out[i]);
    gzwrite(datafile, buffer, strlen(buffer));
  }  
  gzwrite(datafile, "\n", strlen("\n")); 
 
  gzclose(datafile);
  free(filename);
  free(buffer);

  return;
}


/*
 * For evolving h_ij, we work in fourier space, as this has proven to be more stable in past work.
 */
void h_evolve(simType **hij, simType **lij, fftw_complex **fSTTij)
{
  int i, j, k, a, b;
  simType pz, px, py, pp, p2, p;
  simType h, l, S, trs_RE, trs_IM, kkS_RE, kkS_IM;
  simType phat[3], k_mS_mi_RE[3], k_mS_mi_IM[3];

  // evolve h_ij (the fourier transform of)
    // solving (in fourier space) d_t l = k^2 h + S
    //                            d_t h = h


  // The transverse/traceless projection of S is calculated here.
  #pragma omp parallel for default(shared) \
    private(i,j,k,a,b,pz,px,py,pp,p2,p,h,l,S,trs_RE,trs_IM,kkS_RE,kkS_IM,phat,k_mS_mi_RE,k_mS_mi_IM)
  for(int i=0; i<POINTS; i++)
  {
    px = (simType) (i <= POINTS/2 ? i : i - POINTS);
    for(int j=0; j<POINTS; j++)
    {
      py = (simType) (j <= POINTS/2 ? j : j - POINTS);
      for(int k=0; k<POINTS/2+1; k++)
      {
        pz = (simType) k;

        // p (momentum)
        // nan when p=0... oh well.
        p = sqrt(px*px + py*py + pz*pz);
        phat[0] = px/p;
        phat[1] = py/p;
        phat[2] = pz/p;


        // pre-calculate some stuff:   
          // Trace (S^a_a)
          trs_RE = C_RE(fSTTij[0][fSINDEX(i,j,k)]) + C_RE(fSTTij[3][fSINDEX(i,j,k)]) + C_RE(fSTTij[5][fSINDEX(i,j,k)]);
          trs_IM = C_IM(fSTTij[0][fSINDEX(i,j,k)]) + C_IM(fSTTij[3][fSINDEX(i,j,k)]) + C_IM(fSTTij[5][fSINDEX(i,j,k)]);
          // One contraction (k^m S_mi)
          k_mS_mi_RE[0] = phat[0]*C_RE(fSTTij[0][fSINDEX(i,j,k)]) + phat[1]*C_RE(fSTTij[1][fSINDEX(i,j,k)]) + phat[2]*C_RE(fSTTij[2][fSINDEX(i,j,k)]);
          k_mS_mi_RE[1] = phat[0]*C_RE(fSTTij[1][fSINDEX(i,j,k)]) + phat[1]*C_RE(fSTTij[3][fSINDEX(i,j,k)]) + phat[2]*C_RE(fSTTij[4][fSINDEX(i,j,k)]);
          k_mS_mi_RE[2] = phat[0]*C_RE(fSTTij[2][fSINDEX(i,j,k)]) + phat[1]*C_RE(fSTTij[4][fSINDEX(i,j,k)]) + phat[2]*C_RE(fSTTij[5][fSINDEX(i,j,k)]);
          k_mS_mi_IM[0] = phat[0]*C_IM(fSTTij[0][fSINDEX(i,j,k)]) + phat[1]*C_IM(fSTTij[1][fSINDEX(i,j,k)]) + phat[2]*C_IM(fSTTij[2][fSINDEX(i,j,k)]);
          k_mS_mi_IM[1] = phat[0]*C_IM(fSTTij[1][fSINDEX(i,j,k)]) + phat[1]*C_IM(fSTTij[3][fSINDEX(i,j,k)]) + phat[2]*C_IM(fSTTij[4][fSINDEX(i,j,k)]);
          k_mS_mi_IM[2] = phat[0]*C_IM(fSTTij[2][fSINDEX(i,j,k)]) + phat[1]*C_IM(fSTTij[4][fSINDEX(i,j,k)]) + phat[2]*C_IM(fSTTij[5][fSINDEX(i,j,k)]);

          // Two contractions (k^a k^b S_ab)
          kkS_RE = phat[0]*k_mS_mi_RE[0] + phat[1]*k_mS_mi_RE[1] + phat[2]*k_mS_mi_RE[2];
          kkS_IM = phat[0]*k_mS_mi_IM[0] + phat[1]*k_mS_mi_IM[1] + phat[2]*k_mS_mi_IM[2];

        for(a=1; a<=3; a++) {
          for(b=a; b<=3; b++) {
            // (7-a)*a/2-4+b formula maps indexes of h to sequential indices
            // only evolve the transverse/traceless part:
            C_RE(fSTTij[(7-a)*a/2-4+b][fSINDEX(i,j,k)]) +=
              (
                0.5*(phat[a-1]*phat[b-1] - 1.0*(a==b))*trs_RE
                + 0.5*(phat[a-1]*phat[b-1] + 1.0*(a==b))*kkS_RE
                - (phat[a-1]*k_mS_mi_RE[b-1] + phat[b-1]*k_mS_mi_RE[a-1])
              );
            C_IM(fSTTij[(7-a)*a/2-4+b][fSINDEX(i,j,k)]) +=
              (
                0.5*(phat[a-1]*phat[b-1] - 1.0*(a==b))*trs_IM
                + 0.5*(phat[a-1]*phat[b-1] + 1.0*(a==b))*kkS_IM
                - (phat[a-1]*k_mS_mi_IM[b-1] + phat[b-1]*k_mS_mi_IM[a-1])
              );
          }
        }

      } // end k loop
    } // end j
  } // end i

  // No zero mode
  for(a=1; a<=3; a++) {
    for(b=a; b<=3; b++) {
      C_RE(fSTTij[(7-a)*a/2-4+b][fSINDEX(0,0,0)]) = 0;
      C_IM(fSTTij[(7-a)*a/2-4+b][fSINDEX(0,0,0)]) = 0;
    }
  }

  for(a=0; a<12; a++) {
    #pragma omp parallel for default(shared) \
      private(i,j,k,b,pz,px,py,pp,p2,p,h,l,S,trs_RE,trs_IM,kkS_RE,kkS_IM,phat,k_mS_mi_RE,k_mS_mi_IM)
    for(int i=0; i<POINTS; i++)
    {
      px = (simType) (i <= POINTS/2 ? i : i - POINTS);
      for(int j=0; j<POINTS; j++)
      {
        py = (simType) (j <= POINTS/2 ? j : j - POINTS);
        for(int k=0; k<POINTS/2+1; k++)
        {
          pz = (simType) k;

          pp = pow(2.0*M_PI/SIZE, 2)*(px*px + py*py + pz*pz); // k^2; k ~ 2pi/L * i|j|k
          h = hij[a][fSINDEX(i,j,k)];
          l = lij[a][fSINDEX(i,j,k)];

          if(a >= 6) {
            S = C_IM(fSTTij[a-6][fSINDEX(i,j,k)]);
          } else {
            S = C_RE(fSTTij[a][fSINDEX(i,j,k)]);
          }

          // FFT is not yet scaled by measure, so:
          S *= dx*dx*dx;

          // This is RK4. Promise.
          hij[a][fSINDEX(i,j,k)] += dt*(l*(6.0 - pp*dt*dt) + dt*(-h*pp + S)*(12.0 - dt*dt*pp)/4.0)/6.0;
          lij[a][fSINDEX(i,j,k)] += dt*(6.0*S - 3.0*pp*l*dt + pp*pp*l*dt*dt*dt/4.0 - h*pp*(6.0 - pp*dt*dt))/6.0;
        }
      }
    }
  }

  return;
}


/*
 * FFT the STT
 */
void fft_stt(simType **STTij, fftw_complex **fSTTij, fftw_plan p)
{
  // transform STT components into fourier space
  #pragma omp parallel sections num_threads(6)
  {
    
    #pragma omp section
    {
      fftw_execute_dft_r2c(p, STTij[0], fSTTij[0]);
    }
    #pragma omp section
    {
      fftw_execute_dft_r2c(p, STTij[1], fSTTij[1]);
    }
    #pragma omp section
    {
      fftw_execute_dft_r2c(p, STTij[2], fSTTij[2]);
    }
    #pragma omp section
    {
      fftw_execute_dft_r2c(p, STTij[3], fSTTij[3]);
    }
    #pragma omp section
    {
      fftw_execute_dft_r2c(p, STTij[4], fSTTij[4]);
    }
    #pragma omp section
    {
      fftw_execute_dft_r2c(p, STTij[5], fSTTij[5]);
    }
  }
}


/*
 * Calculate S at a point (i,j,k).
 * This does not remove the transverse/traceless part.
 */
void set_stt(PointData *paq, simType **STTij, int i, int j, int k)
{
  // S_ij = T_ij
  //      = (e+p)U^iU^j + d^if*d^jf
  // This calculation is not traceless or transverse - that projection will be done later.
  // However, we don't calculate clearly traceful parts (eg, X*\delta_ij pieces) as a small optimization.

  int a, b;
  for(a=1; a<=3; a++) {
    for(b=a; b<=3; b++) {
      // (7-a)*a/2-4+b formula maps indexes of h to sequential indices
      STTij[(7-a)*a/2-4+b][SINDEX(i,j,k)] =
          paq->gradients[a][4] * paq->gradients[b][4]  // (d^af * d^bf)
          + W_EOSp1 * exp(paq->fields[0]) * paq->fields[a] * paq->fields[b];  // (1+w)(e)U^aU^b
    }
  }

  // project into transverse plane, and subtract remaining trace - this is done during actual evolution.

  return;
}
