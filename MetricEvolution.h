#ifndef METRICEVOLUTION_H
#define METRICEVOLUTION_H

/* Stress Energy Tensor (SET) functionality */
void h_evolve(simType **hij, simType **lij, fftw_complex **fSTTij);
void get_gws(simType **lij);
void fft_stt(simType **STTij, fftw_complex **fSTTij);
void set_stt(PointData *paq, simType **STTij, int i, int j, int k);

#endif
