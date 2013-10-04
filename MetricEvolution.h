#ifndef METRICEVOLUTION_H
#define METRICEVOLUTION_H

/* Stress Energy Tensor (SET) functionality */
void h_evolve(simType **hij, simType **lij, simType **STTij);
void get_gws(PointData *paq, simType **l, int i, int j, int k);
void fft_stt(simType **STTij);
void set_stt(PointData *paq, simType **STTij, int i, int j, int k);

#endif
