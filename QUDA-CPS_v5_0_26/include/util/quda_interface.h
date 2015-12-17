//Begin QUDA_CPS
#ifndef _QUDA_INTERFACE_H_
#define _QUDA_INTERFACE_H_

USING_NAMESPACE_CPS

void fillClover(Lattice &lat, double *h_quda_clover, ChkbType chkb, 
		   DiracOpClover dirac);

void fillCloverInv(Lattice &lat, double *h_quda_clover_inv, ChkbType chkb, 
		   DiracOpClover dirac);

void quda_clover_interface(double *h_quda_clover,  double *h_cps_clover);

void inversion_clover(Lattice &lat, double *h_quda_clover, 
		      double *h_quda_clover_inv, Vector *f_in, Vector *f_out);

void set_quda_params(CgArg *cg_arg);

#endif

//End QUDA_CPS
