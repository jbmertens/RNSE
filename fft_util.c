#include "defines.h"


void fftdump(simType **STTij, fftw_complex **fSTTij, IOData filedata)
{
  // fft_stt(STTij, fSTTij);

  IOData dumpfile;
  dumpfile.data_dir = filedata.data_dir;
  dumpfile.data_name = "FFT_DUMP";
  dumpfile.fwrites = filedata.fwrites;
  dumpfile.datasize = POINTS;
  // dumpstate(storage, dumpfile);

  return;
}

void powerdump(simType *in, fftw_complex *out, fftw_plan plan, IOData filedata)
{

  // Transform input array
  fftw_execute_dft_r2c(plan, in, out);


  // average power over angles
  simType array_out[(int)(1.73205*(POINTS/2))+1];
  int numpoints[(int)(1.73205*(POINTS/2))+1]; // Number of points in each momentum bin
  simType p[(int)(1.73205*(POINTS/2))+1];
  simType f2[(int)(1.73205*(POINTS/2))+1]; // Values for each bin: Momentum, |F-k|^2, n_k
  int numbins = (int)(1.73205*(POINTS/2))+1; // Actual number of bins for the number of dimensions

  double pmagnitude; // Total momentum (p) in units of lattice spacing, pmagnitude = Sqrt(px^2+py^2+pz^2).
                     // This also gives the bin index since bin spacing is set to equal lattice spacing.
  double fp2;
  int i, j, k, px, py, pz; // px, py, and pz are components of momentum in units of grid spacing

  // Initial magnitude of momentum in each bin
  for(i=0; i<numbins; i++) {
    f2[i] = 0.0;
    numpoints[i] = 0;
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
        pmagnitude = sqrt(pw2(px) + pw2(py) + pw2(pz));
        fp2 = pw2(C_RE(out[fSINDEX(i,j,k)])) + pw2(C_IM(out[fSINDEX(i,j,k)]));
        numpoints[(int)pmagnitude] += 2;
        f2[(int)pmagnitude] += 2.*fp2;
      }

      pz = 0;
      k = 0;
      pmagnitude = sqrt(pw2(px) + pw2(py) + pw2( pz));
      fp2 = pw2(C_RE(out[fSINDEX(i,j,k)])) + pw2(C_IM(out[fSINDEX(i,j,k)]));
      numpoints[(int)pmagnitude] += 1;
      f2[(int)pmagnitude] += fp2;
        
      pz = POINTS/2;
      k = POINTS/2;
      pmagnitude = sqrt(pw2(px) + pw2(py) + pw2(pz));
      fp2 = pw2(C_RE(out[fSINDEX(i,j,k)])) + pw2(C_IM(out[fSINDEX(i,j,k)]));
      numpoints[(int)pmagnitude] += 1;
      f2[(int)pmagnitude] += fp2;
    }
  }

  for(i=0; i<numbins; i++)
  {
    // Converts sums to averages. (numpoints[i] should always be greater than zero.)
    if(numpoints[i] > 0)
    {
      array_out[i] = f2[i]/((double) numpoints[i]);
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
  strcat(filename,  ".power.dat.gz");
  buffer = malloc(20 * sizeof(*buffer));

  datafile = (gzFile *)gzopen(filename, "ab");
  if(datafile == Z_NULL) {
    printf("Error opening file: %s\n", filename);
    return;
  }

  for(i=0; i<numbins; i++)
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
