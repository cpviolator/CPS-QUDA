//Begin QUDA_CPS
#ifndef _QUDA_INTERFACE_H_
#define _QUDA_INTERFACE_H_

USING_NAMESPACE_CPS

void fillCloverAll(Lattice &lat, double *h_quda_clover, 
		   double *h_quda_clover_inv, ChkbType chkb, 
		   CgArg *cg_arg, CnvFrmType cnv_frm, DiracOpClover dirac);

void quda_clover_interface(double *h_quda_clover,  double *h_cps_clover);

//Mx=y type probems (correlation functions)
int inversion_clover(Lattice &lat, double *h_quda_clover, 
		     double *h_quda_clover_inv, Vector *f_in, Vector *f_out,
		     Float *QUDA_true_res);
int inversion_wilson(Lattice &lat, Vector *f_in, Vector *f_out,
		     Float *QUDA_true_res);
void set_quda_params(CgArg *cg_arg, int WilClo);

//M^{\dagger}Mx=y type probems (HMC evolution)
int inversion_clover_HMC(Lattice &lat, double *h_quda_clover, 
			 double *h_quda_clover_inv, Vector *f_in, 
			 Vector *f_out, Float *QUDA_true_res);
int inversion_wilson_HMC(Lattice &lat, Vector *f_in, Vector *f_out, 
			 Float *QUDA_true_res);
void set_quda_params_HMC(CgArg *cg_arg, int WilClo);

#endif
//End QUDA_CPS
