/*
This file has the functions needed to evolve the metric perturbations.  It uses a 4th order runge-kutta integrator that ties into the LatticeEasy program.
Written by John T. Giblin, jr .. Last modified 4.22.2007
*/

#include "latticeeasy.h"
#include <omp.h>
const int Block_size = N/8;


//start by defining the same wrapping as in evolution.cpp

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

// delfld calculates the components of the gradient of the field.
// fld corresponds to which field, i,j,k are placement on the lattice and
// n is the direction
// *** needs to be divided by dx to be a real gradient
inline double delfld(int fld, int i, int j, int k, int n) {
#if NDIMS==3
    if (n==1) {
	return (double)(f[fld][INCREMENT(i)][j][k]-f[fld][DECREMENT(i)][j][k])/2;
    }
    else if(n==2) {
	return (double)(f[fld][i][INCREMENT(j)][k]-f[fld][i][DECREMENT(j)][k])/2;
    }
    else if(n==3) {
	return (double)(f[fld][i][j][INCREMENT(k)]-f[fld][i][j][DECREMENT(k)])/2;
    }
    else {
	return 0.;
    }
#endif    
    return 0.;
}


// jtg_stresstensor calculates the T_ab for any point (i,j,k). It's not really the SET, but appropriate combinations of the SET to get S_{ij}--the source term of the metric perturbations
inline double jtg_stresstensor(int aa, int bb, int i, int j, int k) {
    
    double tempfld[nflds];  //to temp. store the field values for calling in 
   
    //the potential energy function
    double deln2 = (double)(1/pw2(dx));  //normalizaion for the gradient 
                                         //since they appear twice, only
                                         //one factor is saved
    
    for(int dd=0; dd<nflds; dd++)
    {
	tempfld[dd] = f[dd][i][j][k];   //save field values in temp.
    }
    
    double tempout=0.0;           //define output variable
    
#if NDIMS==3                      //calculates the correct term
    if (aa==1) {
	if (bb==1) {
	    for(int m=0; m < nflds; m++) {
		tempout += (2./3.)*deln2*delfld(m,i,j,k,1)*delfld(m,i,j,k,1)
		    - (1./3.)*deln2*delfld(m,i,j,k,2)*delfld(m,i,j,k,2)
		    - (1./3.)*deln2*delfld(m,i,j,k,3)*delfld(m,i,j,k,3);
	    }	  
	    return tempout;
	}
	else if (bb==2) {
	    for(int m=0; m < nflds; m++) {
		tempout += deln2*delfld(m,i,j,k,1)*delfld(m,i,j,k,2);
	    }
	    return tempout;
	}
	else if (bb==3) {
	    for(int m=0; m < nflds; m++) {
		tempout += (deln2*delfld(m,i,j,k,1)*delfld(m,i,j,k,3));
	    }
	    return tempout;
	}
	else 
	    printf("calling wrong SET");
	return 0.0;
    }
    else if (aa==2) {
	if (bb==2) {
	    for(int m=0; m < nflds; m++) {
		tempout += -(1./3.)*deln2*delfld(m,i,j,k,1)*delfld(m,i,j,k,1)
		    + (2./3.)*deln2*delfld(m,i,j,k,2)*delfld(m,i,j,k,2)
		    - (1./3.)*deln2*delfld(m,i,j,k,3)*delfld(m,i,j,k,3);
	    }	  
	    return tempout;
	}
	else if (bb==3) {
	    for(int m=0; m < nflds; m++) {
		tempout += (deln2*delfld(m,i,j,k,2)*delfld(m,i,j,k,3));
	    }
	    return tempout;
	}
	else
	    printf("calling wrong SET");
	return 0.0;
    }
    else if (aa==3) {
	if (bb==3) {
	    for(int m=0; m < nflds; m++) {
		tempout += -(1./3.)*deln2*delfld(m,i,j,k,1)*delfld(m,i,j,k,1)
		    - (1./3.)*deln2*delfld(m,i,j,k,2)*delfld(m,i,j,k,2)
		    + (2./3.)*deln2*delfld(m,i,j,k,3)*delfld(m,i,j,k,3);
	    }
	    return tempout;
	}
	else
	    printf("calling wrong SET\n");
	return 0.0;
    }
    else {
	printf("calling wrong SET\n");
	return 0.0;
    }
#else
    printf("calling wrong SET\n");
    return 0.;
#endif
}

//put the sourceterm in momentum space
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
  tij1 = (double (*)[N][N]) malloc(sizeof(double)*N*N*N);
  tij2 = (double (*)[N][N]) malloc(sizeof(double)*N*N*N);
  tij3 = (double (*)[N][N]) malloc(sizeof(double)*N*N*N);
  tij4 = (double (*)[N][N]) malloc(sizeof(double)*N*N*N);
  tij5 = (double (*)[N][N]) malloc(sizeof(double)*N*N*N);

  tij0_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
  tij1_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
  tij2_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
  tij3_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
  tij4_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
  tij5_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
  
int i,j,k;

#pragma omp parallel private(i,j,k)
{
 #pragma omp for schedule(static)
  for(i=0; i<N; i++)  // Fills the arrays with the SET 
  {                
   for (j=0; j<N; j++)            //  adds up diagonal components
   {      
   for (k=0; k<N; k++)
   {         
	tij0[i][j][k] = jtg_stresstensor(1,1,i,j,k);
	tij1[i][j][k] = jtg_stresstensor(1,2,i,j,k);
	tij2[i][j][k] = jtg_stresstensor(1,3,i,j,k);
	tij3[i][j][k] = jtg_stresstensor(2,2,i,j,k);
	tij4[i][j][k] = jtg_stresstensor(2,3,i,j,k);
	tij5[i][j][k] = jtg_stresstensor(3,3,i,j,k);
   }
   }
  }
}
 
  if (!plan_fft) plan_fft = fftw_plan_dft_r2c_3d(N, N, N, &tij0[0][0][0], &tij0_out[0][0][0], FFTW_MEASURE); 
  
#pragma omp parallel sections num_threads(6)
{
 #pragma omp section
  fftw_execute_dft_r2c(plan_fft, &tij0[0][0][0], &tij0_out[0][0][0]);

 #pragma omp section
  fftw_execute_dft_r2c(plan_fft, &tij1[0][0][0], &tij1_out[0][0][0]);

 #pragma omp section
  fftw_execute_dft_r2c(plan_fft, &tij2[0][0][0], &tij2_out[0][0][0]);

 #pragma omp section  
  fftw_execute_dft_r2c(plan_fft, &tij3[0][0][0], &tij3_out[0][0][0]);
 
 #pragma omp section
  fftw_execute_dft_r2c(plan_fft, &tij4[0][0][0], &tij4_out[0][0][0]);

 #pragma omp section
  fftw_execute_dft_r2c(plan_fft, &tij5[0][0][0], &tij5_out[0][0][0]);

} 
  //calculates the FFT
  
  double normmm = pow(dx,3)/(rescale_A*pow(a,2+2.*rescale_s+rescale_r));
  //    double normmm = 1/(rescale_A*pow(a,2+2.*rescale_s+rescale_r));
  //coefficient in front of the SET in the euqations of motion (minus the 
  //  8\p/3)
  
#pragma omp parallel private(i,j,k,px,py,pz,momentum2,dpz,dpx,dpy)
{
  #pragma omp for schedule(static)
  for(i=0; i<N; i++) {              // stores them in the correct
    px = (i<=N/2 ? i : i-N);
    for (j=0; j<N; j++) {         // arrays
      py = (j<=N/2 ? j : j-N);
      for (k=0; k<N/2+1; k++) {
	pz = k;
	momentum2=pw2((double)px)+pw2((double)py)+pw2((double)pz);
	dpz = ((double)pz)/sqrt(momentum2);
	dpx = ((double)px)/sqrt(momentum2);
	dpy = ((double)py)/sqrt(momentum2);

	//T_{11}
	T_gw[0][i][j][k] = .5*pw2(1-dpx*dpx)*normmm*c_re(tij0_out[i][j][k])
	  -(1-dpx*dpx)*dpx*dpy*normmm*c_re(tij1_out[i][j][k])
	  -(1-dpx*dpx)*dpx*dpz*normmm*c_re(tij2_out[i][j][k])
	  +pw2(dpx*dpy)*normmm*c_re(tij3_out[i][j][k])
	  -.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_re(tij3_out[i][j][k])
	  +(1+dpx*dpx)*dpy*dpz*normmm*c_re(tij4_out[i][j][k])
	  +pw2(dpx*dpz)*normmm*c_re(tij5_out[i][j][k])
	  -.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_re(tij5_out[i][j][k]);
	//T_{11}(i)
	T_gw[1][i][j][k] = .5*pw2(1-dpx*dpx)*normmm*c_im(tij0_out[i][j][k])
	  -(1-dpx*dpx)*dpx*dpy*normmm*c_im(tij1_out[i][j][k])
	  -(1-dpx*dpx)*dpx*dpz*normmm*c_im(tij2_out[i][j][k])
	  +pw2(dpx*dpy)*normmm*c_im(tij3_out[i][j][k])
	  -.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_im(tij3_out[i][j][k])
	  +(1+dpx*dpx)*dpy*dpz*normmm*c_im(tij4_out[i][j][k])
	  +pw2(dpx*dpz)*normmm*c_im(tij5_out[i][j][k])
	  -.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_im(tij5_out[i][j][k]);
	//T_{22}
	T_gw[6][i][j][k] = .5*pw2(1-dpy*dpy)*normmm*c_re(tij3_out[i][j][k])
	  -(1-dpy*dpy)*dpx*dpy*normmm*c_re(tij1_out[i][j][k])
	  -(1-dpy*dpy)*dpy*dpz*normmm*c_re(tij4_out[i][j][k])
	  +pw2(dpx*dpy)*normmm*c_re(tij0_out[i][j][k])
	  -.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_re(tij0_out[i][j][k])
	  +(1+dpy*dpy)*dpx*dpz*normmm*c_re(tij2_out[i][j][k])
	  +pw2(dpy*dpz)*normmm*c_re(tij5_out[i][j][k])
	  -.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_re(tij5_out[i][j][k]);
	//T_{22}(i)
	T_gw[7][i][j][k] = .5*pw2(1-dpy*dpy)*normmm*c_im(tij3_out[i][j][k])
	  -(1-dpy*dpy)*dpx*dpy*normmm*c_im(tij1_out[i][j][k])
	  -(1-dpy*dpy)*dpy*dpz*normmm*c_im(tij4_out[i][j][k])
	  +pw2(dpx*dpy)*normmm*c_im(tij0_out[i][j][k])
	  -.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_im(tij0_out[i][j][k])
	  +(1+dpy*dpy)*dpx*dpz*normmm*c_im(tij2_out[i][j][k])
	  +pw2(dpy*dpz)*normmm*c_im(tij5_out[i][j][k])
	  -.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_im(tij5_out[i][j][k]);
	//T_{33}
	T_gw[10][i][j][k] = .5*pw2(1-dpz*dpz)*normmm*c_re(tij5_out[i][j][k])
	  -(1-dpz*dpz)*dpx*dpz*normmm*c_re(tij2_out[i][j][k])
	  -(1-dpz*dpz)*dpy*dpz*normmm*c_re(tij4_out[i][j][k])
	  +pw2(dpx*dpz)*normmm*c_re(tij0_out[i][j][k])
	  -.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_re(tij0_out[i][j][k])
	  +(1+dpz*dpz)*dpx*dpy*normmm*c_re(tij1_out[i][j][k])
	  +pw2(dpy*dpz)*normmm*c_re(tij3_out[i][j][k])
	  -.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_re(tij3_out[i][j][k]);
	//T_{33}(i)
	T_gw[11][i][j][k] = .5*pw2(1-dpz*dpz)*normmm*c_im(tij5_out[i][j][k])
	  -(1-dpz*dpz)*dpx*dpz*normmm*c_im(tij2_out[i][j][k])
	  -(1-dpz*dpz)*dpy*dpz*normmm*c_im(tij4_out[i][j][k])
	  +pw2(dpx*dpz)*normmm*c_im(tij0_out[i][j][k])
	  -.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_im(tij0_out[i][j][k])
	  +(1+dpz*dpz)*dpx*dpy*normmm*c_im(tij1_out[i][j][k])
	  +pw2(dpy*dpz)*normmm*c_im(tij3_out[i][j][k])
	  -.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_im(tij3_out[i][j][k]);
	//T_{12}
	T_gw[2][i][j][k] = -.5*dpx*dpy*(1-dpx*dpx)*normmm*c_re(tij0_out[i][j][k])
	  -.5*dpx*dpy*(1-dpy*dpy)*normmm*c_re(tij3_out[i][j][k])
	  +.5*dpx*dpy*(1+dpz*dpz)*normmm*c_re(tij5_out[i][j][k])
	  +(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_re(tij1_out[i][j][k])
	  -dpy*dpz*(1-dpx*dpx)*normmm*c_re(tij2_out[i][j][k])
	  -dpx*dpz*(1-dpy*dpy)*normmm*c_re(tij4_out[i][j][k]);
	//T_{12}(i)
	T_gw[3][i][j][k] = -.5*dpx*dpy*(1-dpx*dpx)*normmm*c_im(tij0_out[i][j][k])
	  -.5*dpx*dpy*(1-dpy*dpy)*normmm*c_im(tij3_out[i][j][k])
	  +.5*dpx*dpy*(1+dpz*dpz)*normmm*c_im(tij5_out[i][j][k])
	  +(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_im(tij1_out[i][j][k])
	  -dpy*dpz*(1-dpx*dpx)*normmm*c_im(tij2_out[i][j][k])
	  -dpx*dpz*(1-dpy*dpy)*normmm*c_im(tij4_out[i][j][k]);
	//T_{13}
	T_gw[4][i][j][k] = -.5*dpx*dpz*(1-dpx*dpx)*normmm*c_re(tij0_out[i][j][k])
	  +.5*dpx*dpz*(1+dpy*dpy)*normmm*c_re(tij3_out[i][j][k])
	  -.5*dpx*dpz*(1-dpz*dpz)*normmm*c_re(tij5_out[i][j][k])
	  +(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_re(tij2_out[i][j][k])
	  -dpy*dpz*(1-dpx*dpx)*normmm*c_re(tij1_out[i][j][k])
	  -dpx*dpy*(1-dpz*dpz)*normmm*c_re(tij4_out[i][j][k]);
	//T_{13}
	T_gw[5][i][j][k] = -.5*dpx*dpz*(1-dpx*dpx)*normmm*c_im(tij0_out[i][j][k])
	  +.5*dpx*dpz*(1+dpy*dpy)*normmm*c_im(tij3_out[i][j][k])
	  -.5*dpx*dpz*(1-dpz*dpz)*normmm*c_im(tij5_out[i][j][k])
	  +(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_im(tij2_out[i][j][k])
	  -dpy*dpz*(1-dpx*dpx)*normmm*c_im(tij1_out[i][j][k])
	  -dpx*dpy*(1-dpz*dpz)*normmm*c_im(tij4_out[i][j][k]);
	//T_{23}
	T_gw[8][i][j][k] = +.5*dpy*dpz*(1+dpx*dpx)*normmm*c_re(tij0_out[i][j][k])
	  -.5*dpy*dpz*(1-dpy*dpy)*normmm*c_re(tij3_out[i][j][k])
	  -.5*dpy*dpz*(1-dpz*dpz)*normmm*c_re(tij5_out[i][j][k])
	  +(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_re(tij4_out[i][j][k])
	  -dpx*dpy*(1-dpz*dpz)*normmm*c_re(tij2_out[i][j][k])
	  -dpx*dpz*(1-dpy*dpy)*normmm*c_re(tij1_out[i][j][k]);
	//T_{23}(i)
	T_gw[9][i][j][k] = +.5*dpy*dpz*(1+dpx*dpx)*normmm*c_im(tij0_out[i][j][k])
	  -.5*dpy*dpz*(1-dpy*dpy)*normmm*c_im(tij3_out[i][j][k])
	  -.5*dpy*dpz*(1-dpz*dpz)*normmm*c_im(tij5_out[i][j][k])
	  +(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_im(tij4_out[i][j][k])
	  -dpx*dpy*(1-dpz*dpz)*normmm*c_im(tij2_out[i][j][k])
	  -dpx*dpz*(1-dpy*dpy)*normmm*c_im(tij1_out[i][j][k]);
      }
    }
  } 
}



    //subtract zero mode (just in case)
    for(int i=0; i<12; i++){
      T_gw[i][0][0][0] = 0;
    }
   
    //    FILTER(kmin_gw, kmax_gw, T_gw);   
    // allows for a filter on the set

  free(tij0);
  free(tij1);
  free(tij2);
  free(tij3);
  free(tij4);
  free(tij5);
  
  free(tij0_out);
  free(tij1_out);
  free(tij2_out);
  free(tij3_out);
  free(tij4_out);
  free(tij5_out);
}


//  this evolves the metric perturbations using 4th order RK
void evolve_perts(double d) {

    fft_stresstensor();   //calculate the source terms in momentum space
    
    double source[12];

//	printf("called fft_stresstensor \n");

int i,j,k;

#pragma omp parallel private(i,j,k,source)
{
 #pragma omp for schedule(static)
    for(i=0;i<N;i++){        // EOM's
	int ii = (i<=N/2 ? i : i-N);
	for(j=0;j<N;j++){
	    int jj = (j<=N/2 ? j : j-N);
	    for(k=0;k<N/2+1;k++){ //only need to do half the array
		
		source[0] = T_gw[0][i][j][k];      //take the correct parts
		source[1] = T_gw[1][i][j][k];      //point on the lattice
		source[2] = T_gw[2][i][j][k]; 
		source[3] = T_gw[3][i][j][k]; 
		source[4] = T_gw[4][i][j][k]; 
		source[5] = T_gw[5][i][j][k]; 
		source[6] = T_gw[6][i][j][k];    
		source[7] = T_gw[7][i][j][k];    
		source[8] = T_gw[8][i][j][k]; 
		source[9] = T_gw[9][i][j][k]; 
		source[10] = T_gw[10][i][j][k]; 
		source[11] = T_gw[11][i][j][k]; 

		double hd1[12], hd2[12], hd3[12], hd4[12];
		double ld1[12], ld2[12], ld3[12], ld4[12];
		
		derivs(h[0][i][j][k],
		       h[1][i][j][k],
		       h[2][i][j][k],
		       h[3][i][j][k],
		       h[4][i][j][k],
		       h[5][i][j][k],
		       h[6][i][j][k],
		       h[7][i][j][k],
		       h[8][i][j][k],
		       h[9][i][j][k],
		       h[10][i][j][k],
		       h[11][i][j][k],
		       l[0][i][j][k],
		       l[1][i][j][k],
		       l[2][i][j][k],
		       l[3][i][j][k],
		       l[4][i][j][k],
		       l[5][i][j][k],
		       l[6][i][j][k],
		       l[7][i][j][k],
		       l[8][i][j][k],
		       l[9][i][j][k],
		       l[10][i][j][k],
		       l[11][i][j][k],
		       ii, jj, k,
		       hd1, ld1,
		       t, source);
		
		derivs(h[0][i][j][k]+dt*hd1[0]/2,
		       h[1][i][j][k]+dt*hd1[1]/2,
		       h[2][i][j][k]+dt*hd1[2]/2,
		       h[3][i][j][k]+dt*hd1[3]/2,
		       h[4][i][j][k]+dt*hd1[4]/2,
		       h[5][i][j][k]+dt*hd1[5]/2,
		       h[6][i][j][k]+dt*hd1[6]/2,
		       h[7][i][j][k]+dt*hd1[7]/2,
		       h[8][i][j][k]+dt*hd1[8]/2,
		       h[9][i][j][k]+dt*hd1[9]/2,
		       h[10][i][j][k]+dt*hd1[10]/2,
		       h[11][i][j][k]+dt*hd1[11]/2,
		       l[0][i][j][k]+dt*ld1[0]/2,
		       l[1][i][j][k]+dt*ld1[1]/2,
		       l[2][i][j][k]+dt*ld1[2]/2,
		       l[3][i][j][k]+dt*ld1[3]/2,
		       l[4][i][j][k]+dt*ld1[4]/2,
		       l[5][i][j][k]+dt*ld1[5]/2,
		       l[6][i][j][k]+dt*ld1[6]/2,
		       l[7][i][j][k]+dt*ld1[7]/2,
		       l[8][i][j][k]+dt*ld1[8]/2,
		       l[9][i][j][k]+dt*ld1[9]/2,
		       l[10][i][j][k]+dt*ld1[10]/2,
		       l[11][i][j][k]+dt*ld1[11]/2,
		       ii, jj, k,
		       hd2, ld2, 
		       t+dt/2, source);
		
		derivs(h[0][i][j][k]+dt*hd2[0]/2,
		       h[1][i][j][k]+dt*hd2[1]/2,
		       h[2][i][j][k]+dt*hd2[2]/2,
		       h[3][i][j][k]+dt*hd2[3]/2,
		       h[4][i][j][k]+dt*hd2[4]/2,
		       h[5][i][j][k]+dt*hd2[5]/2,
		       h[6][i][j][k]+dt*hd2[6]/2,
		       h[7][i][j][k]+dt*hd2[7]/2,
		       h[8][i][j][k]+dt*hd2[8]/2,
		       h[9][i][j][k]+dt*hd2[9]/2,
		       h[10][i][j][k]+dt*hd2[10]/2,
		       h[11][i][j][k]+dt*hd2[11]/2,
		       l[0][i][j][k]+dt*ld2[0]/2,
		       l[1][i][j][k]+dt*ld2[1]/2,
		       l[2][i][j][k]+dt*ld2[2]/2,
		       l[3][i][j][k]+dt*ld2[3]/2,
		       l[4][i][j][k]+dt*ld2[4]/2,
		       l[5][i][j][k]+dt*ld2[5]/2, 
		       l[6][i][j][k]+dt*ld2[6]/2,
		       l[7][i][j][k]+dt*ld2[7]/2,
		       l[8][i][j][k]+dt*ld2[8]/2,
		       l[9][i][j][k]+dt*ld2[9]/2,
		       l[10][i][j][k]+dt*ld2[10]/2,
		       l[11][i][j][k]+dt*ld2[11]/2,  
		       ii, jj, k,
		       hd3, ld3, 
		       t+dt/2, source);

		derivs(h[0][i][j][k]+dt*hd3[0], 
		       h[1][i][j][k]+dt*hd3[1], 
		       h[2][i][j][k]+dt*hd3[2], 
		       h[3][i][j][k]+dt*hd3[3], 
		       h[4][i][j][k]+dt*hd3[4], 
		       h[5][i][j][k]+dt*hd3[5], 
		       h[6][i][j][k]+dt*hd3[6], 
		       h[7][i][j][k]+dt*hd3[7], 
		       h[8][i][j][k]+dt*hd3[8], 
		       h[9][i][j][k]+dt*hd3[9], 
		       h[10][i][j][k]+dt*hd3[10], 
		       h[11][i][j][k]+dt*hd3[11], 
		       l[0][i][j][k]+dt*ld3[0], 
		       l[1][i][j][k]+dt*ld3[1], 
		       l[2][i][j][k]+dt*ld3[2], 
		       l[3][i][j][k]+dt*ld3[3], 
		       l[4][i][j][k]+dt*ld3[4], 
		       l[5][i][j][k]+dt*ld3[5], 
		       l[6][i][j][k]+dt*ld3[6], 
		       l[7][i][j][k]+dt*ld3[7], 
		       l[8][i][j][k]+dt*ld3[8], 
		       l[9][i][j][k]+dt*ld3[9], 
		       l[10][i][j][k]+dt*ld3[10], 
		       l[11][i][j][k]+dt*ld3[11], 
		       ii, jj, k, 
		       hd4, ld4,
		       t+dt, source);

		h[0][i][j][k] += dt*(hd1[0]/6+hd2[0]/3+hd3[0]/3+hd4[0]/6);
		h[1][i][j][k] += dt*(hd1[1]/6+hd2[1]/3+hd3[1]/3+hd4[1]/6);
		h[2][i][j][k] += dt*(hd1[2]/6+hd2[2]/3+hd3[2]/3+hd4[2]/6);
		h[3][i][j][k] += dt*(hd1[3]/6+hd2[3]/3+hd3[3]/3+hd4[3]/6);
		h[4][i][j][k] += dt*(hd1[4]/6+hd2[4]/3+hd3[4]/3+hd4[4]/6);
		h[5][i][j][k] += dt*(hd1[5]/6+hd2[5]/3+hd3[5]/3+hd4[5]/6);
		h[6][i][j][k] += dt*(hd1[6]/6+hd2[6]/3+hd3[6]/3+hd4[6]/6);
		h[7][i][j][k] += dt*(hd1[7]/6+hd2[7]/3+hd3[7]/3+hd4[7]/6);
		h[8][i][j][k] += dt*(hd1[8]/6+hd2[8]/3+hd3[8]/3+hd4[8]/6);
		h[9][i][j][k] += dt*(hd1[9]/6+hd2[9]/3+hd3[9]/3+hd4[9]/6);
		h[10][i][j][k] += dt*(hd1[10]/6+hd2[10]/3+hd3[10]/3+hd4[10]/6);
		h[11][i][j][k] += dt*(hd1[11]/6+hd2[11]/3+hd3[11]/3+hd4[11]/6);

		l[0][i][j][k] += dt*(ld1[0]/6+ld2[0]/3+ld3[0]/3+ld4[0]/6);
		l[1][i][j][k] += dt*(ld1[1]/6+ld2[1]/3+ld3[1]/3+ld4[1]/6);
		l[2][i][j][k] += dt*(ld1[2]/6+ld2[2]/3+ld3[2]/3+ld4[2]/6);
		l[3][i][j][k] += dt*(ld1[3]/6+ld2[3]/3+ld3[3]/3+ld4[3]/6);
		l[4][i][j][k] += dt*(ld1[4]/6+ld2[4]/3+ld3[4]/3+ld4[4]/6);
		l[5][i][j][k] += dt*(ld1[5]/6+ld2[5]/3+ld3[5]/3+ld4[5]/6);
		l[6][i][j][k] += dt*(ld1[6]/6+ld2[6]/3+ld3[6]/3+ld4[6]/6);
		l[7][i][j][k] += dt*(ld1[7]/6+ld2[7]/3+ld3[7]/3+ld4[7]/6);
		l[8][i][j][k] += dt*(ld1[8]/6+ld2[8]/3+ld3[8]/3+ld4[8]/6);
		l[9][i][j][k] += dt*(ld1[9]/6+ld2[9]/3+ld3[9]/3+ld4[9]/6);
		l[10][i][j][k] += dt*(ld1[10]/6+ld2[10]/3+ld3[10]/3+ld4[10]/6);
		l[11][i][j][k] += dt*(ld1[11]/6+ld2[11]/3+ld3[11]/3+ld4[11]/6);

		//R-K routine

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


//function for calculating the derivatives
void derivs(double h11, 
	    double hi11, 
	    double h12, 
	    double hi12, 
	    double h13, 
	    double hi13, 
	    double h22, 
	    double hi22, 
	    double h23, 
	    double hi23, 
	    double h33, 
	    double hi33,
	    double l11, 
	    double li11, 
	    double l12, 
	    double li12, 
	    double l13, 
	    double li13, 
	    double l22, 
	    double li22, 
	    double l23, 
	    double li23, 
	    double l33, 
	    double li33,	    int ii, int jj, int k,
	    double hd[12], double ld[12], 
	    double time, double source_gw[12]) {
    
    double norm = (pw2(2.*pi)/(pw2(L)*pow(a,2.*rescale_s+2)));
    double omega = (8.*pi);
    double kk = (pw2((double)ii) + pw2((double)jj) + pw2((double)k))*norm;

    ld[0] = - kk*h11 + 2.*omega*source_gw[0];
    ld[1] = - kk*hi11 + 2.*omega*source_gw[1];
    ld[2] = - kk*h12 + 2.*omega*source_gw[2];
    ld[3] = - kk*hi12 + 2.*omega*source_gw[3];
    ld[4] = - kk*h13 + 2.*omega*source_gw[4];
    ld[5] = - kk*hi13 + 2.*omega*source_gw[5];
    ld[6] = - kk*h22 + 2.*omega*source_gw[6];
    ld[7] = - kk*hi22 + 2.*omega*source_gw[7];
    ld[8] = - kk*h23 + 2.*omega*source_gw[8];
    ld[9] = - kk*hi23 + 2.*omega*source_gw[9];
    ld[10] = - kk*h33 + 2.*omega*source_gw[10];
    ld[11] = - kk*hi33 + 2.*omega*source_gw[11];

    hd[0]=l11;
    hd[1]=li11;
    hd[2]=l12;
    hd[3]=li12;
    hd[4]=l13;
    hd[5]=li13;
    hd[6]=l22;
    hd[7]=li22;
    hd[8]=l23;
    hd[9]=li23;
    hd[10]=l33;
    hd[11]=li33;

}


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
	    
	    pz=0;
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
