#include <config.h>
#include <util/vector.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <comms/scu.h>
#include <util/qcdio.h>
#include <alg/qpropw.h>
//#include <alg/alg_ghb.h>
//#include <alg/ghb_arg.h>
#include <alg/alg_hmd.h>                                                                                                             
#include <alg/alg_pbp.h>
#include <alg/alg_w_spect.h>
#include <alg/alg_plaq.h>
#include <alg/do_arg.h>
#include <alg/no_arg.h>
#include <alg/alg_eig.h>

#include <alg/fix_gauge_arg.h>  
#include <alg/alg_fix_gauge.h> 
#include <util/random.h>

#include <util/WriteLatticePar.h>
#include <util/ReadLatticePar.h>

#include <cstdlib>
#include <cstdio>

#define MAX_NUM_TSITES 64

USING_NAMESPACE_CPS

// Local function prototypes.                                                            
void setup_do_arg(DoArg& do_arg, int seed, Float BETA);
void setup_hmd_arg(HmdArg& hmd_arg, Float mass, Float tau);
Float MMDag_re_tr(WilsonMatrix WM);


#define NSITES_3D 4
#define NSITES_T 8
//#define BETA 5.45
//#define MASS -0.775
#define STOP_RSD 1.0e-6
#define MAX_NUM_ITER 50000

#define TAU_INIT 0.5
#define TAU_SIM 0.5
#define TAU_STEPS 10

//if READ==1, read in lattice configs.
//#define READ 
#define SMEAR 1 // if SMEAR==1, use smearing.
#define NDATA 1
#define NSKIP 10
#define NTHERM 10

int gaugecounter = 1;

int main(int argc, char *argv[]) {

  Start(&argc,&argv);
  int seed = atoi(argv[1]);
  cout<<seed<<endl;
  double MASS = atof(argv[2]);
  //MASS *= -1.0/20;
  cout<<MASS<<endl;

  double BETA = 5.60;
  DoArg do_arg;
  setup_do_arg(do_arg,seed,BETA);
  GJP.Initialize(do_arg);

  GwilsonFclover lat;
  CommonArg c_arg;

  //Declare args for Gaussian Smearing
  QPropWGaussArg g_arg;
  g_arg.gauss_link_smear_type=GKLS_STOUT; //Link smearing
  g_arg.gauss_link_smear_coeff=0.15;      //Link smearing
  g_arg.gauss_link_smear_N=3;             //Link smearing hits
  g_arg.gauss_N=8;                               //Source/Sink smearing hits
  g_arg.gauss_W=sqrt(0.125*4*g_arg.gauss_N);     //Smearing parameter
  
  int sweep_counter = 0;
  int total_updates = NTHERM + NSKIP*(NDATA-1);

  double mass = MASS;
  HmdArg hmd_arg;
  char filename[256];
  sprintf(filename, "HMC_data_%.2f_%.3f.dat", BETA, MASS);
  c_arg.results = CAST_AWAY_CONST(filename);
  setup_hmd_arg(hmd_arg, mass, TAU_INIT);
  AlgHmcPhi hmc(lat, &c_arg, &hmd_arg);

  QPropWArg arg0;
  arg0.t=0;
  arg0.x=0;
  arg0.y=0;
  arg0.z=0;
  arg0.cg.mass = MASS;
  arg0.cg.stop_rsd = STOP_RSD;
  arg0.cg.max_num_iter = MAX_NUM_ITER;
  arg0.cg.Inverter = BICGSTAB;
  arg0.cg.bicgstab_n = 1;

  int x2[4];

  WilsonMatrix t4;
  WilsonMatrix t4c;
		
  Float d0_t4t4c_re_tr = 0.0;

  //////////////////////
  // Start simulation //
  ////////////////////// 

  while (sweep_counter < total_updates) {
    for (int n = 1; n <= NSKIP; n++) {
#ifdef READ == 1
      // do nothing
#else
      hmc.run();
#endif
      sweep_counter++;
      if (!UniqueID()) {
        printf("step %d complete\n",sweep_counter);
        fflush(stdout);
      }
    }

    if (sweep_counter == NTHERM) printf("thermalization complete. \n");
    if (sweep_counter >= NTHERM) {

      // Use this code to specify a gauge configuration.
      char lattice[256];
      sprintf(lattice,"lat_HMC_B%.2f_M%.3f_%d-%d_%d.dat", BETA, MASS, NSITES_3D, NSITES_T, sweep_counter);
      //sprintf(lattice,"/home/howarth/latt_configs_BGQ/lat_HMC_%d-%d_%d.dat", NSITES_3D, NSITES_T, sweep_counter);
      //sprintf(lattice,"/home/howarth/latt_configs_BGQ/lat_HMC_B%.2f_M%.3f_%d-%d_%d.dat", BETA, MASS, NSITES_3D, NSITES_T, sweep_counter);
#ifdef READ == 1
      ReadLatticeParallel(lat,lattice);
#else
      //WriteLatticeParallel(lat,lattice);
#endif      
      gaugecounter = 1;

      //int P[3]={0,0,0};
      //Get Momentum Propagator	    
      //QPropWMomSrc qprop_mom(lat, &arg_z, P, &c_arg);
      //QPropWMomSrcSmeared qprop_mom(lat, &arg0, P, &g_arg, &c_arg);
      //qprop_mom.GaussSmearSinkProp(g_arg);
      
      QPropWGaussSrc qprop0(lat, &arg0, &g_arg, &c_arg);
      qprop0.GaussSmearSinkProp(g_arg);
      char file[256];
      //sprintf(file, "B%.2f_M%.3f_N%d_W%.3f_1pion_HMC_smear_%d-%d.dat", BETA, MASS, g_arg.gauss_N, g_arg.gauss_W, NSITES_3D, NSITES_T);
      sprintf(file, "CLOVER_FRM_STANDARD_CPU_B%.2f_M%.3f_N%d_W%.3f_1pion_HMC_smear_%d-%d.dat",BETA,MASS, g_arg.gauss_N, g_arg.gauss_W, NSITES_3D, NSITES_T);
      
      //Sum over x2
      for (x2[3]=0; x2[3]<GJP.TnodeSites(); x2[3]++) {
	//Reinitialise trace
	d0_t4t4c_re_tr *= 0.0;	
	for (x2[2]=0; x2[2]<GJP.ZnodeSites(); x2[2]++)
	  for (x2[1]=0; x2[1]<GJP.YnodeSites(); x2[1]++)
	    for (x2[0]=0; x2[0]<GJP.XnodeSites(); x2[0]++) {
	      int x2_idx = lat.GsiteOffset(x2)/4;
	      
	      t4 = qprop0[x2_idx];
	      //t4c = t4.conj_cp();
	      d0_t4t4c_re_tr += MMDag_re_tr(t4);
	    }
	
	//////////////////////////
	// Write trace to file. //
	//////////////////////////
	
	FILE *t4tr=Fopen(file,"a");
	Fprintf(t4tr,"%d %d %d %.16e\n", sweep_counter, x2[3], 0, d0_t4t4c_re_tr);
	Fclose(t4tr);
	cout<<"time slice = "<<x2[3]<<" complete."<<endl;
	
	//////////////////////////////////////////
	// End trace summation at time slice t. //
	//////////////////////////////////////////
	
      }
    }	
  }
  
  ////////////////////
  // End simulation //
  ////////////////////
  
  End();
  return 0;
}

void setup_do_arg(DoArg& do_arg, int seed, Float BETA)
{

  do_arg.x_node_sites = NSITES_3D/SizeX();
  do_arg.y_node_sites = NSITES_3D/SizeY();
  do_arg.z_node_sites = NSITES_3D/SizeZ();
  do_arg.t_node_sites = NSITES_T/SizeT();
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_INPUT;
  do_arg.start_seed_value = seed;
  do_arg.beta = BETA;

  do_arg.clover_coeff = 1.0;
  do_arg.xi_bare = 1;
  do_arg.xi_v = 1;
  do_arg.xi_v_xi = 1;

}

void setup_hmd_arg(HmdArg& hmd_arg, Float mass, Float tau)
{
  hmd_arg.n_frm_masses = 1;
  hmd_arg.frm_mass[0] = mass;
  hmd_arg.field_type[0] = FERMION;
  hmd_arg.frm_flavors[0] = 2;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.steps_per_traj = TAU_STEPS;
  hmd_arg.step_size = tau/hmd_arg.steps_per_traj;
  hmd_arg.metropolis = METROPOLIS_YES;
  //  hmd_arg.metropolis = METROPOLIS_NO;
  hmd_arg.reunitarize = REUNITARIZE_YES;
  //  hmd_arg.isz = 0; // Location of smallest polar shift
  //  hmd_arg.sw = 2; // Sexton-Weingarten term (gauge contribution per fermion)
  hmd_arg.max_num_iter[0] = 5000;
  hmd_arg.stop_rsd[0] = 1.0e-7;
  //  hmd_arg.stop_rsd_md[0] = 1e-6;
  //  hmd_arg.stop_rsd_mc[0] = 1e-10;
  //  hmd_arg.valid_approx[0] = 0;
  
  hmd_arg.hmd_link_smear_N = 1;
  hmd_arg.hmd_link_smear_coeff = 0.01;
  //hmd_arg.hmd_link_smear_type = HMDLS_STOUT;
}

Float MMDag_re_tr(WilsonMatrix WM) {

  int s1 = 0;
  int c1 = 0;
  int s2 = 0;
  int c2 = 0;
  int sc_idx = 0;

  Float re_tr = 0.0;

  for(s1=0;s1<4;s1++)
    for(c1=0;c1<3;c1++)
      for(s2=0;s2<4;s2++)
	for(c2=0;c2<3;c2++){
	  sc_idx = 2*(c2 + 3*s2 + 12*c1 + 36*s1);

	  re_tr += WM(s1,c1,s2,c2).real()*WM(s1,c1,s2,c2).real() + WM(s1,c1,s2,c2).imag()*WM(s1,c1,s2,c2).imag();
	}
  return re_tr;
}
