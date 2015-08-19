#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<string>

#include <fftw3.h>

#include "mom1D.h"

#define DATAPATH "../data/"

using namespace std;

void FFT_wtf_ZYX(fftw_complex *FT_9d, int diag, int SINPz_Pz, int SINPxy_Pxy,
int NMOM, int Lmax, int T, int sweep_counter, int t)
{

  char file[256];
  sprintf(file, DATAPATH"D%d_msFT_%d_%d-%d_%d_%d.dat", diag, NMOM, Lmax, T, sweep_counter, t);
  FILE *fp_9d   = fopen(file, "a");

  //Normalised momentum selection: sin(|p|)/|p| > SINP_P/100

  double SIN_Z_cutoff  = SINPz_Pz/(1.0*100);
  double SIN_XY_cutoff = SINPxy_Pxy/(1.0*100);

  mom1D z(Lmax, SIN_Z_cutoff);
  mom1D y(Lmax, SIN_XY_cutoff);
  mom1D x(Lmax, SIN_XY_cutoff);

  int mom_arr_idx = 0;

  for (z.P[2]=0; z.P[2] < Lmax; z.P[2]++)
    for (z.P[1]=0; z.P[1] < Lmax; z.P[1]++)
      for (z.P[0]=0; z.P[0] < Lmax; z.P[0]++)

  for (y.P[2]=0; y.P[2] < Lmax; y.P[2]++)
    for (y.P[1]=0; y.P[1] < Lmax; y.P[1]++)
      for (y.P[0]=0; y.P[0] < Lmax; y.P[0]++)

  for (x.P[2]=0; x.P[2] < Lmax; x.P[2]++)
    for (x.P[1]=0; x.P[1] < Lmax; x.P[1]++)
      for (x.P[0]=0; x.P[0] < Lmax; x.P[0]++) {

        mom_arr_idx = x.index() + pow(Lmax,3)*y.index() + pow(Lmax,6)*z.index();

          if( (sin(z.mod())/z.mod() > SIN_Z_cutoff || z.mod() == 0) &&
          (sin(y.mod())/y.mod() > SIN_XY_cutoff || y.mod() == 0)    &&
          (sin(x.mod())/x.mod() > SIN_XY_cutoff || x.mod() == 0)) {

            fprintf(fp_9d, "%d %d %d %d %d %d %d %d %d %d %d %d %.16e %.16e\n", sweep_counter, t,
            mom_arr_idx, z.P[2], z.P[1], z.P[0], y.P[2], y.P[1], y.P[0], x.P[2], x.P[1], x.P[0],
            FT_9d[mom_arr_idx][0]/NMOM, FT_9d[mom_arr_idx][1]/NMOM);

          }
        }
  fclose(fp_9d);
}


void FFT_wtf_ZX(fftw_complex *FT_6d, int diag, int SINPz_Pz, int SINPxy_Pxy,
 int NMOM, int Lmax, int T, int sweep_counter, int t){
  
  char file[256];
  sprintf(file, DATA_PATH"D%d_msFT_%d_%d-%d_%d_%d.dat", diag, NMOM, Lmax, T, sweep_counter, t);
  FILE *fp_6d   = fopen(file, "a");
  
  //Normalised momentum selection: sin(|p|)/|p| > SINP_P/100
  
  double SIN_Z_cutoff  = SINPz_Pz/(1.0*100);
  double SIN_XY_cutoff = SINPxy_Pxy/(1.0*100);
 
  mom1D z(Lmax, SIN_Z_cutoff);
  mom1D x(Lmax, SIN_XY_cutoff);
 
  int mom_arr_idx = 0;
  
  for (z.P[2]=0; z.P[2] < Lmax; z.P[2]++)
    for (z.P[1]=0; z.P[1] < Lmax; z.P[1]++)
      for (z.P[0]=0; z.P[0] < Lmax; z.P[0]++)

  for (x.P[2]=0; x.P[2] < Lmax; x.P[2]++)
    for (x.P[1]=0; x.P[1] < Lmax; x.P[1]++)
      for (x.P[0]=0; x.P[0] < Lmax; x.P[0]++) {

        mom_arr_idx = x.index() + pow(Lmax,3)*z.index();
	      
	if( (sin(z.mod())/z.mod() > SIN_Z_cutoff || z.mod() == 0) &&
	(sin(x.mod())/x.mod() > SIN_XY_cutoff || x.mod() == 0)) {
		
	  fprintf(fp_6d, "%d %d %d %d %d %d %d %d %d %.16e %.16e\n", sweep_counter, t,
	  mom_arr_idx, z.P[2], z.P[1], z.P[0], x.P[2], x.P[1], x.P[0],
  	  FT_6d[mom_arr_idx][0]/NMOM, FT_6d[mom_arr_idx][1]/NMOM);

	}
      }
  fclose(fp_6d);
}

void FFT_wtf_XY(fftw_complex *FT_6d, int diag, int SINPz_Pz, int SINPxy_Pxy,
int NMOM, int Lmax, int T, int sweep_counter, int t){
  
  char file[256];
  sprintf(file, DATA_PATH"D%d_msFT_%d_%d-%d_%d_%d.dat", diag, NMOM, Lmax, T, sweep_counter, t);
  FILE *fp_6d   = fopen(file, "a");
  
  //Normalised momentum selection: sin(|p|)/|p| > SINP_P/100
  
  //double SIN_Z_cutoff  = SINPz_Pz/(1.0*100);
  double SIN_XY_cutoff = SINPxy_Pxy/(1.0*100);

  mom1D y(Lmax, SIN_XY_cutoff);
  mom1D x(Lmax, SIN_XY_cutoff);

  int mom_arr_idx = 0;

  for (y.P[2]=0; y.P[2] < Lmax; y.P[2]++)
    for (y.P[1]=0; y.P[1] < Lmax; y.P[1]++)
      for (y.P[0]=0; y.P[0] < Lmax; y.P[0]++)

  for (x.P[2]=0; x.P[2] < Lmax; x.P[2]++)
    for (x.P[1]=0; x.P[1] < Lmax; x.P[1]++)
      for (x.P[0]=0; x.P[0] < Lmax; x.P[0]++) {
  
	mom_arr_idx = x.index() + pow(Lmax,3)*y.index();
	      
	if( (sin(y.mod())/y.mod() > SIN_XY_cutoff || y.mod() == 0) &&
	(sin(x.mod())/x.mod() > SIN_XY_cutoff || x.mod() == 0)) {
		
	  fprintf(fp_6d, "%d %d %d %d %d %d %d %d %d %.16e %.16e\n", sweep_counter, t, 
	  mom_arr_idx, y.P[2], y.P[1], y.P[0], x.P[2], x.P[1], x.P[0],
	  FT_6d[mom_arr_idx][0]/NMOM, FT_6d[mom_arr_idx][1]/NMOM);

	}
      }
  fclose(fp_6d);
}

void print_FFTW_arr(fftw_complex *arr, int Lmax) {
  for (int i = 0; i < pow(Lmax, 3); i++) {
    printf("%d: %.6e  %.6e\n", i, arr[i][0], arr[i][1]);
  }
}

//backwards FFT
void FFT_B(int dim, int Lmax, fftw_complex *arr) {

  int N[dim];
  for (int a=0; a<dim; a++) N[a] = Lmax;  

  //Create plan
  fftw_plan plan = fftw_plan_dft(dim, N, arr, arr, FFTW_BACKWARD, FFTW_MEASURE);
  
  //print_FFTW_arr(arr, L);
  fftw_execute(plan);
  //print_FFTW_arr(arr, L);

  fftw_destroy_plan(plan);
}

//forwards FFT
void FFT_F(int dim, int L, fftw_complex *arr) {
  
  int N[dim];
  for (int a=0; a<dim; a++) N[a] = L;  

  //Create plan
  fftw_plan plan = fftw_plan_dft(dim, N, arr, arr, FFTW_FORWARD, FFTW_MEASURE);
  
  //print_FFTW_arr(arr, L);
  fftw_execute(plan);
  //print_FFTW_arr(arr, L);

  fftw_destroy_plan(plan);
}
