// includes CUDA Runtime
#include <cuda.h>
#include <cuda_runtime.h>
//#include <cuda_functions.h>
//#include <cuda_runtime_api.h>

//input A,B; output Tr(A*Bdag)
__global__ void tr_ABdag(double *re_tr_arr, double *A, double *Bdag) {
  
  int snk_idx_vol = threadIdx.x;  //Sink index
  
  double A_arr[288];
  for(int a=0; a<288; a++) A_arr[a] = A[288*(snk_idx_vol) + a];
  

  // Perform trace sum.
  //
  // Tr(ABdag) = A_ab * Bdag_ba 
  //
  // N.B. Bdag enters the function as B. We perform
  // conjugation by manipulating B's row and column
  // indices and multiplying the imaginary elements
  // of B by -1.0 as required.
  // This has the neat result:
  // Tr(AG5_BdagG5) = sum_n sum_m a_nm*b_nm
  
  int s1 = 0;
  int c1 = 0;
  int s2 = 0;
  int c2 = 0;
  int sc_idx = 0;
  
  for(s1=0;s1<4;s1++)
    for(c1=0;c1<3;c1++)
      for(s2=0;s2<4;s2++)
 	for(c2=0;c2<3;c2++) {
 	  //REAL[(reA + im1A*i)(re2B - imB*i)] = (reAreB + imAimB)                                             
 	  sc_idx = 2*(c2 + (3*s2) + (12*c1) + (36*s1));
	  
	  re_tr_arr[snk_idx_vol] += A_arr[sc_idx] * A_arr[sc_idx] + A_arr[sc_idx+1]*A_arr[sc_idx+1];
 	}
}
