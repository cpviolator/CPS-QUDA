//Begin QUDA_CPS
#ifdef USEQUDA
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>

enum {
  CLOVER_MAT_SIZE = 72,   // COLORS^2 * lat.SpinComponents()^2 / 2;
  HALF_CLOVER_MAT_SIZE = 36
};

#include <config.h>
#include <util/quda_interface.h>
#include <quda.h>

extern int gaugecounter;

QudaInvertParam inv_param;
QudaGaugeParam param;
QudaPrecision cpu_prec;

void fillCloverAll(Lattice &lat, double *h_quda_clover, 
		   double *h_quda_clover_inv, ChkbType chkb, CgArg *cg_arg, 
		   CnvFrmType cnv_frm, PreserveType prs_f_in, 
		   DiracOpClover dirac)
{
  double *ptr_clover_quda_inv = h_quda_clover_inv;
  double *ptr_clover_quda = h_quda_clover;
  
  dirac.CloverMatChkb(chkb, 0);
  double *ptr_clover_cps = (double *)(chkb==CHKB_EVEN ? 
				      lat.Aux0Ptr() : lat.Aux1Ptr());
  for(int i=0;i<GJP.VolNodeSites()/2;i++){
    quda_clover_interface(ptr_clover_quda, ptr_clover_cps);
    ptr_clover_quda += CLOVER_MAT_SIZE;
    ptr_clover_cps += CLOVER_MAT_SIZE;
  }

  ChkbType chkb_opp = (chkb == CHKB_EVEN ? CHKB_ODD : CHKB_EVEN );
  dirac.CloverMatChkb(chkb_opp, 0);
  ptr_clover_cps = (double *)(chkb_opp == CHKB_EVEN ? 
			      lat.Aux0Ptr() : lat.Aux1Ptr());
  // ptr_clover_quda starts the loop pointing at the 
  // second half of the array.
  for(int i=0;i<GJP.VolNodeSites()/2;i++){
    quda_clover_interface(ptr_clover_quda, ptr_clover_cps);
    ptr_clover_quda += CLOVER_MAT_SIZE;
    ptr_clover_cps += CLOVER_MAT_SIZE;
  }

  dirac.CloverMatChkb(chkb, 1);
  double *ptr_clover_cps_inv = (double *)(chkb==CHKB_EVEN ? 
					  lat.Aux0Ptr() : lat.Aux1Ptr());
  for(int i=0;i<GJP.VolNodeSites()/2;i++){
    quda_clover_interface(ptr_clover_quda_inv, ptr_clover_cps_inv);
    ptr_clover_quda_inv += CLOVER_MAT_SIZE;
    ptr_clover_cps_inv += CLOVER_MAT_SIZE;
  }

  dirac.CloverMatChkb(chkb_opp, 1);
  ptr_clover_cps_inv = (double *)(chkb_opp == CHKB_EVEN ? 
				  lat.Aux0Ptr() : lat.Aux1Ptr());
  // ptr_clover_quda_inv starts the loop pointing at the 
  // second half of the array.
  for(int i=0;i<GJP.VolNodeSites()/2;i++){
    quda_clover_interface(ptr_clover_quda_inv, ptr_clover_cps_inv);
    ptr_clover_quda_inv += CLOVER_MAT_SIZE;
    ptr_clover_cps_inv += CLOVER_MAT_SIZE;
  }
}

void inversion_clover(Lattice &lat, double *h_quda_clover, 
		      double *h_quda_clover_inv, Vector *f_in, Vector *f_out)
{ 
  int f_size_cb = GJP.VolNodeSites() * lat.FsiteSize();
  void *h_gauge=(void*)lat.GaugeField();
  double *gauge_reord=(double*)smalloc(4*18*GJP.VolNodeSites()*sizeof(double));
  
  reorder_gauge(h_gauge, gauge_reord);
  
  // 1 = initialize, else just calculate.
  if (gaugecounter == 1){
    freeGaugeQuda();
    loadGaugeQuda(gauge_reord, &param);
    freeCloverQuda();
    loadCloverQuda(h_quda_clover, h_quda_clover_inv, &inv_param);
    gaugecounter=0;
  }
  
  double *bQUDA=(double*)smalloc(f_size_cb*sizeof(double));
  double *pb=(double*)f_in;
  for (int i=0;i<(f_size_cb/2);i++) {
    bQUDA[i]=pb[i+f_size_cb/2];
    bQUDA[i+f_size_cb/2]=pb[i];
  }
  
  double *xQUDA=(double*)smalloc(f_size_cb*sizeof(double));
  double *px=(double*)f_out;
  for (int i=0;i<(f_size_cb/2);i++) {
    xQUDA[i]=px[i+f_size_cb/2];
    xQUDA[i+f_size_cb/2]=px[i];
  }
  invertQuda((void*)xQUDA, (void*)bQUDA, &inv_param);
  for (int i=0;i<(f_size_cb/2);i++) {
    px[i]=xQUDA[i+f_size_cb/2];
    px[i+f_size_cb/2]=xQUDA[i];
  }

  sfree(gauge_reord);
  sfree(bQUDA);
  sfree(xQUDA);
  
}

void reorder_gauge(void *h_gauge, double *gauge_reord)
{
  int g_size = 4*GJP.VolNodeSites();
  
  //Transpose
  double *gauge_tran=(double*)smalloc(4*18*GJP.VolNodeSites()*sizeof(double));
  double *ptr_orig=(double*)h_gauge;
  double *ptr_tran=gauge_tran;
  for(int i=0;i<(18*g_size);i=i+18){
    for(int j=0;j<2;j++){
      ptr_tran[i+j]    = ptr_orig[i+j];
      ptr_tran[i+2+j]  = ptr_orig[i+6+j];
      ptr_tran[i+4+j]  = ptr_orig[i+12+j];
      ptr_tran[i+6+j]  = ptr_orig[i+2+j];
      ptr_tran[i+8+j]  = ptr_orig[i+8+j];
      ptr_tran[i+10+j] = ptr_orig[i+14+j];
      ptr_tran[i+12+j] = ptr_orig[i+4+j];
      ptr_tran[i+14+j] = ptr_orig[i+10+j];
      ptr_tran[i+16+j] = ptr_orig[i+16+j];
    }
  }

  //Hermitean Conjugate
  double *gauge_herm=(double*)smalloc(4*18*GJP.VolNodeSites()*sizeof(double));
  double *ptr_herm=gauge_herm;  
  for(int i=0;i<(18*g_size);i++){
    if(!(i%2)) ptr_herm[i]=ptr_tran[i];
    else ptr_herm[i]=(-1)*ptr_tran[i];
  }

  //Reorder the gauge array
  double *ptr_reord = gauge_reord;
  int g_even = 0;
  int g_odd  = 18*4*GJP.VolNodeSites()/2;

  int T = 0;
  int Z = 0;
  int Y = 0;
  int X = 0;
  int index = 0;

  for(T=0 ; T<GJP.TnodeSites();T++)
    for(Z=0 ; Z<GJP.ZnodeSites();Z++)
      for(Y=0 ; Y<GJP.YnodeSites();Y++)
        for(X=0 ; X<GJP.XnodeSites();X++){
          index = (X+GJP.XnodeSites()*(Y+GJP.YnodeSites()
				       *(Z+GJP.ZnodeSites()*T)));
          if(!((X+Y+Z+T)%2)){
            for(int l=0;l<18*4;l++)
	      ptr_reord[g_even+l] = ptr_tran[72*index+l];	    
            g_even += 18*4;
          }
          else{
            for(int l=0;l<18*4;l++)
	      ptr_reord[g_odd+l] = ptr_tran[72*index+l];
            g_odd += 18*4;
          }
        }
  
  sfree(gauge_tran);
  sfree(gauge_herm);
  
}

void quda_clover_interface(double *h_quda_clover, double *h_cps_clover)
{
  h_quda_clover[0]=h_cps_clover[0];    // c00_00_re = C0.x, A0
  h_quda_clover[1]=h_cps_clover[3];    // c01_01_re = C0.y, A0
  h_quda_clover[2]=h_cps_clover[8];    // c02_02_re = C0.z, A0
  h_quda_clover[3]=h_cps_clover[15];   // c10_10_re = C0.w, A2
  h_quda_clover[4]=h_cps_clover[24];   // c11_11_re = C1.x, A2
  h_quda_clover[5]=h_cps_clover[35];   // c12_12_re = C1.y, A2
  h_quda_clover[6]=h_cps_clover[1];    // c01_00_re = C1.z, A0
  h_quda_clover[7]=h_cps_clover[2];    // c01_00_im = C1.w, A0
  h_quda_clover[8]=h_cps_clover[4];    // c02_00_re = C2.x, A0
  h_quda_clover[9]=h_cps_clover[5];    // c02_00_im = C2.y, A0
  h_quda_clover[10]=h_cps_clover[9];   // c10_00_re = C2.z, A1
  h_quda_clover[11]=h_cps_clover[10];  // c10_00_im = C2.w, A1
  h_quda_clover[12]=h_cps_clover[16];  // c11_00_re = C3.x, A1
  h_quda_clover[13]=h_cps_clover[17];  // c11_00_im = C3.y, A1
  h_quda_clover[14]=h_cps_clover[25];  // c12_00_re = C3.z, A1
  h_quda_clover[15]=h_cps_clover[26];  // c12_00_im = C3.w, A1
  h_quda_clover[16]=h_cps_clover[6];   // c02_01_re = C4.x, A0
  h_quda_clover[17]=h_cps_clover[7];   // c02_01_im = C4.y, A0
  h_quda_clover[18]=h_cps_clover[11];  // c10_01_re = C4.z, A1
  h_quda_clover[19]=h_cps_clover[12];  // c10_01_im = C4.w, A1
  h_quda_clover[20]=h_cps_clover[18];  // c11_01_re = C5.x, A1
  h_quda_clover[21]=h_cps_clover[19];  // c11_01_im = C5.y, A1
  h_quda_clover[22]=h_cps_clover[27];  // c12_01_re = C5.z, A1
  h_quda_clover[23]=h_cps_clover[28];  // c12_01_im = C5.w, A1
  h_quda_clover[24]=h_cps_clover[13];  // c10_02_re = C6.x, A1
  h_quda_clover[25]=h_cps_clover[14];  // c10_02_im = C6.y, A1
  h_quda_clover[26]=h_cps_clover[20];  // c11_02_re = C6.z, A1
  h_quda_clover[27]=h_cps_clover[21];  // c11_02_im = C6.w, A1
  h_quda_clover[28]=h_cps_clover[29];  // c12_02_re = C7.x, A1
  h_quda_clover[29]=h_cps_clover[30];  // c12_02_im = C7.y, A1
  h_quda_clover[30]=h_cps_clover[22];  // c11_10_re = C7.z, A2
  h_quda_clover[31]=h_cps_clover[23];  // c11_10_im = C7.w, A2
  h_quda_clover[32]=h_cps_clover[31];  // c12_10_re = C8.x, A2
  h_quda_clover[33]=h_cps_clover[32];  // c12_10_im = C8.y, A2
  h_quda_clover[34]=h_cps_clover[33];  // c12_11_re = C8.z, A2
  h_quda_clover[35]=h_cps_clover[34];  // c12_11_im = C8.w, A2
  // Simple offset for second chiral block:
  h_quda_clover[0+36]=h_cps_clover[0+36];    // c20_20_re = C0.x, A3
  h_quda_clover[1+36]=h_cps_clover[3+36];    // c21_21_re = C0.y, A3
  h_quda_clover[2+36]=h_cps_clover[8+36];    // c22_22_re = C0.z, A3
  h_quda_clover[3+36]=h_cps_clover[15+36];   // c30_30_re = C0.w, A5
  h_quda_clover[4+36]=h_cps_clover[24+36];   // c31_31_re = C1.x, A5
  h_quda_clover[5+36]=h_cps_clover[35+36];   // c32_32_re = C1.y, A5
  h_quda_clover[6+36]=h_cps_clover[1+36];    // c21_20_re = C1.z, A3
  h_quda_clover[7+36]=h_cps_clover[2+36];    // c21_20_im = C1.w, A3
  h_quda_clover[8+36]=h_cps_clover[4+36];    // c22_20_re = C2.x, A3
  h_quda_clover[9+36]=h_cps_clover[5+36];    // c22_20_im = C2.y, A3
  h_quda_clover[10+36]=h_cps_clover[9+36];   // c30_20_re = C2.z, A4
  h_quda_clover[11+36]=h_cps_clover[10+36];  // c30_20_im = C2.w, A4
  h_quda_clover[12+36]=h_cps_clover[16+36];  // c31_20_re = C3.x, A4
  h_quda_clover[13+36]=h_cps_clover[17+36];  // c31_20_im = C3.y, A4
  h_quda_clover[14+36]=h_cps_clover[25+36];  // c32_20_re = C3.z, A4
  h_quda_clover[15+36]=h_cps_clover[26+36];  // c32_20_im = C3.w, A4
  h_quda_clover[16+36]=h_cps_clover[6+36];   // c22_21_re = C4.x, A3
  h_quda_clover[17+36]=h_cps_clover[7+36];   // c22_21_im = C4.y, A3
  h_quda_clover[18+36]=h_cps_clover[11+36];  // c30_21_re = C4.z, A4
  h_quda_clover[19+36]=h_cps_clover[12+36];  // c30_21_im = C4.w, A4
  h_quda_clover[20+36]=h_cps_clover[18+36];  // c31_21_re = C5.x, A4
  h_quda_clover[21+36]=h_cps_clover[19+36];  // c31_21_im = C5.y, A4
  h_quda_clover[22+36]=h_cps_clover[27+36];  // c32_21_re = C5.z, A4
  h_quda_clover[23+36]=h_cps_clover[28+36];  // c32_21_im = C5.w, A4
  h_quda_clover[24+36]=h_cps_clover[13+36];  // c30_22_re = C6.x, A4
  h_quda_clover[25+36]=h_cps_clover[14+36];  // c30_22_im = C6.y, A4
  h_quda_clover[26+36]=h_cps_clover[20+36];  // c31_22_re = C6.z, A4
  h_quda_clover[27+36]=h_cps_clover[21+36];  // c31_22_im = C6.w, A4
  h_quda_clover[28+36]=h_cps_clover[29+36];  // c32_22_re = C7.x, A4
  h_quda_clover[29+36]=h_cps_clover[30+36];  // c32_22_im = C7.y, A4
  h_quda_clover[30+36]=h_cps_clover[22+36];  // c31_20_re = C7.z, A5
  h_quda_clover[31+36]=h_cps_clover[23+36];  // c31_20_im = C7.w, A5
  h_quda_clover[32+36]=h_cps_clover[31+36];  // c32_30_re = C8.x, A5
  h_quda_clover[33+36]=h_cps_clover[32+36];  // c32_30_im = C8.y, A5
  h_quda_clover[34+36]=h_cps_clover[33+36];  // c32_31_re = C8.z, A5
  h_quda_clover[35+36]=h_cps_clover[34+36];  // c32_31_im = C8.w, A5
}

void set_quda_params(CgArg *cg_arg)
{
  ///////////////////////////////////////////////////////////////////
  // For the reconstruct type we must be careful. We must use 
  // QUDA_RECONSTRUCT_12 until we better understand the 
  // unpacking algorithm.
  ///////////////////////////////////////////////////////////////////

  param = newQudaGaugeParam();

  // set the CUDA precisions
  param.reconstruct = QUDA_RECONSTRUCT_12;
  param.cuda_prec = QUDA_DOUBLE_PRECISION;

  // set the CUDA sloppy precisions  
  param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;
  param.cuda_prec_sloppy = QUDA_DOUBLE_PRECISION;

  // set the CUDA precondition precisions
  param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
  param.reconstruct_precondition = QUDA_RECONSTRUCT_12;
  param.cpu_prec = QUDA_DOUBLE_PRECISION; 
  
  param.X[0] = GJP.XnodeSites();
  param.X[1] = GJP.YnodeSites();
  param.X[2] = GJP.ZnodeSites();
  param.X[3] = GJP.TnodeSites();

  param.anisotropy = GJP.XiBare();
  param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
  param.gauge_order = QUDA_CPS_WILSON_GAUGE_ORDER;
  param.gauge_fix = QUDA_GAUGE_FIXED_NO;
  param.type = QUDA_WILSON_LINKS;  
  param.t_boundary = QUDA_PERIODIC_T;  // or QUDA_ANTI_PERIODIC_T
  param.ga_pad = 0; 

  //////////////////////////////////////////////////////////
  //               QUDA MATRIX INVERSION                  //
  //////////////////////////////////////////////////////////

  inv_param = newQudaInvertParam();

  inv_param.clover_cpu_prec = QUDA_DOUBLE_PRECISION; //QUDA precision type
  inv_param.dslash_type = QUDA_CLOVER_WILSON_DSLASH;
  
  //Options: QUDA_CG_INVERTER, QUDA_GCR_INVERTER, QUDA_BICGSTAB_INVERTER
  switch(cg_arg->Inverter){    
  case CG : 
    inv_param.inv_type       = QUDA_CG_INVERTER;
    printf("QUDA CG \n");
    break;
  case BICGSTAB : 
    inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
    printf("QUDA BICGSTAB \n");
    break;
  case QUDA_GCR : 
    inv_param.inv_type = QUDA_GCR_INVERTER;
    inv_param.gcrNkrylov = 30;
    printf("QUDA GCR \n");
    break;
  default : 
    inv_param.inv_type   = QUDA_GCR_INVERTER;
    inv_param.gcrNkrylov = 30; 
    printf("QUDA DEFAULT = GCR\n");
  }

  inv_param.clover_cuda_prec = QUDA_DOUBLE_PRECISION;
  inv_param.cuda_prec = QUDA_DOUBLE_PRECISION;
  inv_param.cuda_prec_sloppy = QUDA_DOUBLE_PRECISION;
  inv_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
  inv_param.clover_cuda_prec_sloppy = QUDA_DOUBLE_PRECISION;
  inv_param.clover_cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
  inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
  
  inv_param.twist_flavor = QUDA_TWIST_NO;
  inv_param.maxiter = cg_arg->max_num_iter; //dirac_arg.max_num_iter;
  inv_param.reliable_delta = 1e-3;          //QudaParam.reliable_delta;
  inv_param.tol = cg_arg->stop_rsd;
  inv_param.Ls = 0;
  inv_param.cl_pad = 0;
  inv_param.sp_pad = 0;

  double mass = cg_arg->mass;
  inv_param.kappa = 1.0/(2.0*(4.0+mass));
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  inv_param.dagger = QUDA_DAG_NO;
  //changed from QUDA_MATPCDAG_MATPC_SOLUTION to current.
  inv_param.solution_type = QUDA_MAT_SOLUTION;
  // Options: QUDA_DIRECT_SOLVE, QUDA_NORMOP_SOLVE, QUDA_DIRECT_PC_SOLVE, QUDA_NORMOP_PC_SOLVE, QUDA_INVALID_SOLVE
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE; 
  inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS; 
  // Options: QUDA_CPS_WILSON_DIRAC_ORDER
  inv_param.dirac_order = QUDA_DIRAC_ORDER; 
  inv_param.tune = QUDA_TUNE_NO;
  inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;
  // QUDA_VERBOSE, SILENT, DEBUG_VERBOSE, SUMMARIZE
  inv_param.verbosity = QUDA_SUMMARIZE; 
  // QUDA_PACKED_CLOVER_ORDER: even-odd, packed
  // QUDA_LEX_PACKED_CLOVER_ORDER: lexicographical order, packed 
  inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
}
#endif
//End QUDA_CPS
