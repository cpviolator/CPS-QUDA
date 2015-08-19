#include "args.h"

void setup_do_arg(DoArg& do_arg, int seed)
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

void setup_hmd_arg(HmdArg& hmd_arg)
{

  hmd_arg.n_frm_masses = 1;
  hmd_arg.frm_mass[0] = MASS;
  hmd_arg.field_type[0] = FERMION;
  hmd_arg.frm_flavors[0] = 2;
  hmd_arg.n_bsn_masses = 0;
  hmd_arg.steps_per_traj = TAU_STEPS;
  hmd_arg.step_size = TAU_INIT/hmd_arg.steps_per_traj;
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
}

void setup_qpropwarg_cg(QPropWArg& q_arg) {
  q_arg.cg.mass = MASS;
  q_arg.cg.stop_rsd = STOP_RSD;
  q_arg.cg.max_num_iter = MAX_NUM_ITER;
  q_arg.cg.Inverter = INVERTER_TYPE;
  q_arg.cg.bicgstab_n = BICGSTAB_N;
}

void setup_g_arg(QPropWGaussArg& g_arg) {
  g_arg.gauss_link_smear_type = GAUSS_LS_TYPE;    //Link smearing
  g_arg.gauss_link_smear_coeff = GAUSS_LS_COEFF;  //Link smearing
  g_arg.gauss_link_smear_N = GAUSS_LS_N;          //Link smearing hits
  g_arg.gauss_N = GAUSS_N;                        //Source/Sink smearing hits
  g_arg.gauss_W = sqrt(KAPPA*4*g_arg.gauss_N);    //Smearing parameter.
}

