/*
This file has the functions needed to evolve the metric perturbations.  It uses a 4th order runge-kutta integrator that ties into the LatticeEasy program.
Written by John T. Giblin, jr .. Last modified 4.22.2007
*/
#include "latticeeasy.h"
#include <omp.h>
const int Block_size = N/8;


// Increments a grid location accounting for periodic wrapping
inline int INCREMENT(int i)
{
  return( (i==N-1) ? 0 : i+1 );
}


// Decrements a grid location accounting for periodic wrapping
inline int DECREMENT(int i)
{
  return( (i==0) ? N-1 : i-1 );
}


/*
 * delfld calculates the components of the gradient of the field.
 * fld corresponds to which field, i,j,k are placement on the lattice and
 * n is the direction
 * *** needs to be divided by dx to be a real gradient
 */
inline double delfld(int fld, int i, int j, int k, int n) {
  /* calculate derivative */
  if (n==1) {
    return (double)(f[fld][INCREMENT(i)][j][k]-f[fld][DECREMENT(i)][j][k])/2;
  }
  else if(n==2)
  {
    return (double)(f[fld][i][INCREMENT(j)][k]-f[fld][i][DECREMENT(j)][k])/2;
  }
  else if(n==3)
  {
    return (double)(f[fld][i][j][INCREMENT(k)]-f[fld][i][j][DECREMENT(k)])/2;
  }
  else
  {
    return 0.;
  }

  return 0.;
}


/*
 * jtg_stresstensor calculates the T_ab for any point (i,j,k).
 */
inline double jtg_stresstensor(int aa, int bb, int i, int j, int k) {
  /* return stress energy tensor at point i, j, k */
  return 0.0;
}


/*
 * put the source term in momentum space
 */
void fft_stresstensor() {

  int pz, px, py;
  double dpx, dpz, dpy;
  double momentum2;

  //need to declare the objects to use FFTW

  double (*tij0)[N][N], (*tij1)[N][N], (*tij2)[N][N];
  double (*tij3)[N][N], (*tij4)[N][N], (*tij5)[N][N];

  fftw_complex (*tij0_out)[N][N/2+1], (*tij1_out)[N][N/2+1], (*tij2_out)[N][N/2+1];
  fftw_complex (*tij3_out)[N][N/2+1], (*tij4_out)[N][N/2+1], (*tij5_out)[N][N/2+1];

  tij0 = (double (*)[N][N]) malloc(sizeof(double)*N*N*N);
  /**[CLIP]**/

  tij0_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
  /**[CLIP]**/

  int i,j,k;

  #pragma omp parallel private(i,j,k)
  {
   #pragma omp for schedule(static)
    for(i=0; i<N; i++)  // Fills the arrays with the SET 
    {                
      for (j=0; j<N; j++)  //  adds up diagonal components
      {      
        for (k=0; k<N; k++)
        {         
          tij0[i][j][k] = jtg_stresstensor(1,1,i,j,k);
          /**[CLIP]**/
        }
      }
    }
  }

  if (!plan_fft) plan_fft = fftw_plan_dft_r2c_3d(N, N, N, &tij0[0][0][0], &tij0_out[0][0][0], FFTW_MEASURE); 

  #pragma omp parallel sections num_threads(6)
  {
      #pragma omp section
      fftw_execute_dft_r2c(plan_fft, &tij0[0][0][0], &tij0_out[0][0][0]);
      /**[CLIP]**/
  } 
  //calculates the FFT
  
  //coefficient in front of the SET in the euqations of motion (minus the 8pi/3)
  double normmm = pow(dx,3)/(rescale_A*pow(a,2+2.*rescale_s+rescale_r));
  
  #pragma omp parallel private(i,j,k,px,py,pz,momentum2,dpz,dpx,dpy)
  {
    #pragma omp for schedule(static)
    for(i=0; i<N; i++)
    {              // stores them in the correct
      px = (i<=N/2 ? i : i-N);
      for (j=0; j<N; j++)
      {         // arrays
        py = (j<=N/2 ? j : j-N);
        for (k=0; k<N/2+1; k++)
        {
          pz = k;
          momentum2=pw2((double)px)+pw2((double)py)+pw2((double)pz);
          dpz = ((double)pz)/sqrt(momentum2);
          dpx = ((double)px)/sqrt(momentum2);
          dpy = ((double)py)/sqrt(momentum2);

          // CALCULATE components of T
          //T_{11}
          T_gw[0][i][j][k] = .5*pw2(1-dpx*dpx)*normmm*c_re(tij0_out[i][j][k])
            -(1-dpx*dpx)*dpx*dpy*normmm*c_re(tij1_out[i][j][k])
            -(1-dpx*dpx)*dpx*dpz*normmm*c_re(tij2_out[i][j][k])
            +pw2(dpx*dpy)*normmm*c_re(tij3_out[i][j][k])
            -.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_re(tij3_out[i][j][k])
            +(1+dpx*dpx)*dpy*dpz*normmm*c_re(tij4_out[i][j][k])
            +pw2(dpx*dpz)*normmm*c_re(tij5_out[i][j][k])
            -.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_re(tij5_out[i][j][k]);
        }
      }
    } 
  }

  //subtract zero mode (just in case)
  for(int i=0; i<12; i++)
  {
    T_gw[i][0][0][0] = 0;
  }

  free(tij0);
  /**[CLIP]**/
  free(tij0_out);
  /**[CLIP]**/
}


//  this evolves the metric perturbations using 4th order RK
void evolve_perts(double d) {

  fft_stresstensor();   //calculate the source terms in momentum space
  double source[12];
  int i,j,k;

  #pragma omp parallel private(i,j,k,source)
  {
    #pragma omp for schedule(static)
    for(i=0;i<N;i++)
    {        // EOM's
      int ii = (i<=N/2 ? i : i-N);
      for(j=0;j<N;j++)
      {
        int jj = (j<=N/2 ? j : j-N);
        for(k=0;k<N/2+1;k++)
        { //only need to do half the array

          source[0] = T_gw[0][i][j][k];      //take the correct parts
          /**[CLIP]**/

          double hd1[12], hd2[12], hd3[12], hd4[12];
          double ld1[12], ld2[12], ld3[12], ld4[12];
          
          derivs(h[0][i][j][k], /**[CLIP]**/
                 l[0][i][j][k], /**[CLIP]**/
                 ii, jj, k,
                 hd1, ld1,
                 t, source);
          
          derivs(h[0][i][j][k]+dt*hd1[0]/2, /**[CLIP]**/
                 l[0][i][j][k]+dt*ld1[0]/2, /**[CLIP]**/
                 ii, jj, k,
                 hd2, ld2, 
                 t+dt/2, source);
          
          derivs(h[0][i][j][k]+dt*hd2[0]/2, /**[CLIP]**/
                 l[0][i][j][k]+dt*ld2[0]/2, /**[CLIP]**/
                 ii, jj, k,
                 hd3, ld3, 
                 t+dt/2, source);

          derivs(h[0][i][j][k]+dt*hd3[0], /**[CLIP]**/
                 l[0][i][j][k]+dt*ld3[0], /**[CLIP]**/
                 ii, jj, k, 
                 hd4, ld4,
                 t+dt, source);

          h[0][i][j][k] += dt*(hd1[0]/6+hd2[0]/3+hd3[0]/3+hd4[0]/6);
          /**[CLIP]**/
          l[0][i][j][k] += dt*(ld1[0]/6+ld2[0]/3+ld3[0]/3+ld4[0]/6);
          /**[CLIP]**/
        }
      }
    }
  }

  // file for testing in case you want to output any variables--these are some i found usefull
  FILE *tomout;
  tomout = fopen("tomout.txt", "a");
  fprintf(tomout, "%10.10f   %10.10f \n", t, h[0][1][5][6]);

  fclose (tomout);
}


/*
 * Mutator function for calculating derivatives
 */
void derivs(double h11, 
        /**[CLIP]**/
        double li33,        int ii, int jj, int k,
        double hd[12], double ld[12], 
        double time, double source_gw[12])
{
    double norm = (pw2(2.*pi)/(pw2(L)*pow(a,2.*rescale_s+2)));
    double omega = (8.*pi);
    double kk = (pw2((double)ii) + pw2((double)jj) + pw2((double)k))*norm;

    ld[0] = - kk*h11 + 2.*omega*source_gw[0];
    /**[CLIP]**/
    hd[0]=l11;
    /**[CLIP]**/
}


/*
 * calculate GW spectrum
 */
void jtg_spectrum() {

    double array_out[(int)(1.73205*(N/2))+1];
    int numpoints_gw[(int)(1.73205*(N/2))+1]; // Number of points in each momentum bin
    double p[(int)(1.73205*(N/2))+1];
    double f2_gw[(int)(1.73205*(N/2))+1]; // Values for each bin: Momentum, |F-k|^2, n_k
    int numbins_gw=(int)(sqrt((double)NDIMS)*(N/2))+1; // Actual number of bins for the number of dimensions

    double pmagnitude_gw; // Total momentum (p) in units of lattice spacing, pmagnitude = Sqrt(px^2+py^2+pz^2). This also gives the bin index since bin spacing is set to equal lattice spacing.
    double dp_gw=2.*pi/(double)L*rescale_B*6.0e10/a/sqrt(rescale_B*pow(a,rescale_s)*hub); // Size of grid spacing in momentum space
    double fp2_gw;
    int i,j,k,px,py,pz; // px, py, and pz are components of momentum in units of grid spacing
    double norm1_gw=4.e-5/pow(100.,.333)/pow(a,2.*rescale_r)*pi/3/pw2(rescale_A)/pow(L,3.)/pw2(hub);

    // Calculate magnitude of momentum in each bin
    for(i=0;i<numbins_gw;i++) {
      p[i]=dp_gw*(double)i;
      f2_gw[i]=0.0;
      numpoints_gw[i]=0;
    }

    for(i=0;i<N;i++) {
      px = (i<=N/2 ? i : i-N);
      for(j=0;j<N;j++) {
        py = (j<=N/2 ? j : j-N);
        for(k=1;k<N/2;k++) { 
          pz = k;
          pmagnitude_gw = sqrt(pw2(px)+pw2(py)+pw2(pz));
          fp2_gw = pw2(l[0][i][j][k]-rescale_r*hub*h[0][i][j][k]) 
            + pw2(l[1][i][j][k]-rescale_r*hub*h[1][i][j][k])
            + 2.*pw2(l[2][i][j][k]-rescale_r*hub*h[2][i][j][k]) 
            + 2.*pw2(l[3][i][j][k]-rescale_r*hub*h[3][i][j][k])
            + 2.*pw2(l[4][i][j][k]-rescale_r*hub*h[4][i][j][k]) 
            + 2.*pw2(l[5][i][j][k]-rescale_r*hub*h[5][i][j][k])
            + pw2(l[6][i][j][k]-rescale_r*hub*h[6][i][j][k]) 
            + pw2(l[7][i][j][k]-rescale_r*hub*h[7][i][j][k])
            + 2.*pw2(l[8][i][j][k]-rescale_r*hub*h[8][i][j][k]) 
            + 2.*pw2(l[9][i][j][k]-rescale_r*hub*h[9][i][j][k])
            + pw2(l[10][i][j][k]-rescale_r*hub*h[10][i][j][k]) 
            + pw2(l[11][i][j][k]-rescale_r*hub*h[11][i][j][k]);
          numpoints_gw[(int)pmagnitude_gw] += 2;
          f2_gw[(int)pmagnitude_gw] += 2.*fp2_gw;
        }
        
        pz = 0;
        pmagnitude_gw=sqrt(pw2(px)+pw2(py)+pw2(pz));
        fp2_gw = pw2(l[0][i][j][k]-rescale_r*hub*h[0][i][j][k]) 
          + pw2(l[1][i][j][k]-rescale_r*hub*h[1][i][j][k])
          + 2.*pw2(l[2][i][j][k]-rescale_r*hub*h[2][i][j][k]) 
          + 2.*pw2(l[3][i][j][k]-rescale_r*hub*h[3][i][j][k])
          + 2.*pw2(l[4][i][j][k]-rescale_r*hub*h[4][i][j][k]) 
          + 2.*pw2(l[5][i][j][k]-rescale_r*hub*h[5][i][j][k])
          + pw2(l[6][i][j][k]-rescale_r*hub*h[6][i][j][k]) 
          + pw2(l[7][i][j][k]-rescale_r*hub*h[7][i][j][k])
          + 2.*pw2(l[8][i][j][k]-rescale_r*hub*h[8][i][j][k]) 
          + 2.*pw2(l[9][i][j][k]-rescale_r*hub*h[9][i][j][k])
          + pw2(l[10][i][j][k]-rescale_r*hub*h[10][i][j][k]) 
          + pw2(l[11][i][j][k]-rescale_r*hub*h[11][i][j][k]);
        
        numpoints_gw[(int)pmagnitude_gw] += 1;
        f2_gw[(int)pmagnitude_gw] += fp2_gw;
        
        pz = N/2;
        pmagnitude_gw=sqrt(pw2(px)+pw2(py)+pw2(pz));
        fp2_gw = pw2(l[0][i][j][k]-rescale_r*hub*h[0][i][j][k]) 
          + pw2(l[1][i][j][k]-rescale_r*hub*h[1][i][j][k])
          + 2.*pw2(l[2][i][j][k]-rescale_r*hub*h[2][i][j][k]) 
          + 2.*pw2(l[3][i][j][k]-rescale_r*hub*h[3][i][j][k])
          + 2.*pw2(l[4][i][j][k]-rescale_r*hub*h[4][i][j][k]) 
          + 2.*pw2(l[5][i][j][k]-rescale_r*hub*h[5][i][j][k])
          + pw2(l[6][i][j][k]-rescale_r*hub*h[6][i][j][k]) 
          + pw2(l[7][i][j][k]-rescale_r*hub*h[7][i][j][k])
          + 2.*pw2(l[8][i][j][k]-rescale_r*hub*h[8][i][j][k]) 
          + 2.*pw2(l[9][i][j][k]-rescale_r*hub*h[9][i][j][k])
          + pw2(l[10][i][j][k]-rescale_r*hub*h[10][i][j][k]) 
          + pw2(l[11][i][j][k]-rescale_r*hub*h[11][i][j][k]);

        numpoints_gw[(int)pmagnitude_gw] += 1;
        f2_gw[(int)pmagnitude_gw] += fp2_gw;
      }
    }
    
    for(i=0;i<numbins_gw;i++) {
      if(numpoints_gw[i]>0) {// Converts sums to averages. (numpoints[i] should always be greater than zero.)
          array_out[i] = norm1_gw*pow((double)i/L,3.)*f2_gw[i]/(double)numpoints_gw[i];
      }
      else {
          array_out[i] = 0.;
      }
    }
  
    
    FILE *gwspec;
    gwspec = fopen("gwspec.txt", "a");
    
    for(i=0;i<numbins_gw;i++) {
      fprintf(gwspec,"%10.10f %10.10f  %50.50f  %10.10f  %10.10f \n", t, p[i], array_out[i], ad/a, a);
    }
    fclose (gwspec);

    return;
}
