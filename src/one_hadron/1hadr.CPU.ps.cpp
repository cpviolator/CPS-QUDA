//local include directory
#include <CPS_utils.h>

//here
#include "params.h"
#include "args.h"

USING_NAMESPACE_CPS

int gaugecounter = 1;

int main(int argc, char *argv[]) {

  int seed = atoi(argv[1]);

  Start(&argc,&argv);

  DoArg do_arg;
  setup_do_arg(do_arg, seed);
  GJP.Initialize(do_arg);

  //VRB.DeactivateAll();
  
  GwilsonFclover lat;
  CommonArg c_arg;

  //Declare args for Gaussian Smearing
  QPropWGaussArg g_arg;
  setup_g_arg(g_arg);

  char is_qu[5];
  #ifdef QUENCH
    GhbArg ghb_arg;
    ghb_arg.num_iter = 1;
    AlgGheatBath hb(lat, &c_arg, &ghb_arg);
    strcpy(is_qu,"QUEN");
  #else
    HmdArg hmd_arg;
    setup_hmd_arg(hmd_arg);
    AlgHmcPhi hmc(lat, &c_arg, &hmd_arg);
    strcpy(is_qu,"UNQU");
  #endif

  int sweep_counter = 0;
  int total_updates = NTHERM + NSKIP*(NDATA-1);

  QPropWArg arg0;
  setup_qpropwarg_cg(arg0);
  arg0.t=0;
  arg0.x=0;
  arg0.y=0;
  arg0.z=0;

  int x2[4];
  WilsonMatrix t4;		
  Float d0_t4t4c_re_tr = 0.0;
  int x2_idx = 0;
  int vol3d = pow(NSITES_3D,3);
  char lattice[256]; //lattice config file
  char file[256];  //output file

  //////////////////////
  // Start simulation //
  ////////////////////// 

  while (sweep_counter < total_updates) {
    for (int n = 1; n <= NSKIP; n++) {
#ifdef READ
      //do nothing
#else
      #ifdef QUENCH 
	hb.run();
      #else 
	hmc.run();
      #endif
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
      #ifdef QUENCH
        sprintf(lattice, LATT_PATH"QU/lat_hb_B%.2f_%d-%d_%d.dat", BETA, NSITES_3D, NSITES_T, sweep_counter);
      #else
	sprintf(lattice, LATT_PATH"UNQ/lat_hmc_B%.2f_M%.3f_%d-%d_%d.dat", BETA, NSITES_3D, NSITES_T, sweep_counter);
      #endif
#ifdef READ
      ReadLatticeParallel(lat,lattice);
#else
      WriteLatticeParallel(lat,lattice);
#endif
      gaugecounter = 1;

      // Get Point Source Propagator	    
      // This will place a point source (Kronecker Delta) at the coordinates
      // specified by arg0. It will then be smeared using the parameters specified
      // by g_arg.
      QPropWGaussSrc qprop0(lat, &arg0, &g_arg, &c_arg);
      // Smear the sink with the same g_arg parameters.
      qprop0.GaussSmearSinkProp(g_arg);
      
      //Sum over x2
      for (x2[3]=0; x2[3]<GJP.TnodeSites(); x2[3]++) {
	//Reinitialise trace
	d0_t4t4c_re_tr *= 0.0;	
	for (x2[2]=0; x2[2]<GJP.ZnodeSites(); x2[2]++)
	  for (x2[1]=0; x2[1]<GJP.YnodeSites(); x2[1]++)
	    for (x2[0]=0; x2[0]<GJP.XnodeSites(); x2[0]++) {
	      x2_idx = lat.GsiteOffset(x2)/4;

	      //Get propagator sinked at x2.
	      t4 = qprop0[x2_idx];
	      //Get the real part of the trace.
	      d0_t4t4c_re_tr += MMDag_re_tr(t4);
	    }
	
	//////////////////////////
	// Write trace to file. //
	//////////////////////////
	
	
	//Write data file so that the data can be reproduced from the name of the file.
	sprintf(file, DATAPATH"PS_CPU_%d_B%.2f_M%.3f_N%d_W%.3f_n%d_xi%.2f_1pion_%s_stout_%d-%d.dat",
	seed, BETA, MASS, g_arg.gauss_N, g_arg.gauss_W, g_arg.gauss_link_smear_N, 
	g_arg.gauss_link_smear_coeff, is_qu, NSITES_3D, NSITES_T);
	
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
