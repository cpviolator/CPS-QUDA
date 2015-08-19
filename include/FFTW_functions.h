#include <fftw3.h>

#ifndef _FFTW_FUNCTIONS_H_
#define _FFTW_FUNCTIONS_H_

void print_FFTW_arr(fftw_complex *FFTW_arr, int Lmax);
void FFT_F(int dim, int Lmax, fftw_complex *arr);
void FFT_B(int dim, int Lmax, fftw_complex *arr);
void FFT_wtf_ZYX(fftw_complex *FT_9d, int diag, int SINPz_Pz, int SINPxy_Pxy,
	int NMOM, int Lmax, int T, int sweep_counter, int t);
void FFT_wtf_ZX(fftw_complex *FT_6d, int diag, int SINPz_Pz, int SINPxy_Pxy,
	int NMOM, int Lmax, int T, int sweep_counter, int t);
void FFT_wtf_XY(fftw_complex *FT_6d, int diag, int SINPz_Pz, int SINPxy_Pxy,
	int NMOM, int Lmax, int T, int sweep_counter, int t);

#endif
