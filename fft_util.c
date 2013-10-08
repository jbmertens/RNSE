#include "defines.h"


void fftdump(simType **STTij, fftw_complex **fSTTij, IOData filedata)
{
  fft_stt(STTij, fSTTij);

  IOData dumpfile;
  dumpfile.data_dir = filedata.data_dir;
  dumpfile.data_name = "FFT_DUMP";
  dumpfile.fwrites = filedata.fwrites;
  dumpfile.datasize = POINTS;
  // dumpstate(storage, dumpfile);

  return;
}
