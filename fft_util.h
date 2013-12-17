#ifndef FFT_H
#define FFT_H

void fftdump(simType **STTij, fftw_complex **fSTTij, IOData filedata);

void powerdump(simType *in, fftw_complex *out, fftw_plan plan, IOData filedata);

#endif
