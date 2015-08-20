//here
#include "params.h"
#include "args.h"

USING_NAMESPACE_CPS

// We will compute the trace of all graphs
// over all lattice sites
// using a point sources at 0 and Z

int gaugecounter = 1;

int main(int argc, char *argv[]) {

  Start(&argc,&argv);
  int seed = atoi(argv[1]);         //
  int SINPz_Pz   = atof(argv[2]);   // integer percentage of the tolerance of sin(p)/p at Z.
  int SINPxy_Pxy = atof(argv[3]);   // integer percentage of the tolerance of sin(p)/p at XY.
  //int t_in = atoi(argv[5]);         //

  DoArg do_arg;
  setup_do_arg(do_arg, seed); 
  GJP.Initialize(do_arg);  

  GwilsonFclover lat;
  CommonArg c_arg;

  //Declare args for Gaussian Smearing
  QPropWGaussArg g_arg_mom;
  setup_g_arg(g_arg_mom);


  int sweep_counter = 0;
  int total_updates = NTHERM + NSKIP*(NDATA-1);

  #ifdef QUENCH
  GhbArg ghb_arg;
  ghb_arg.num_iter = 1;
  AlgGheatBath hb(lat, &c_arg, &ghb_arg);
  #else
  HmdArg hmd_arg; 
  setup_hmd_arg(hmd_arg);
  AlgHmcPhi hmc(lat, &c_arg, &hmd_arg); 
  #endif

  //Declare args for source at 0.
  QPropWArg arg_0;
  setup_qpropwarg_cg(arg_0);
  arg_0.x = 0;
  arg_0.y = 0;
  arg_0.z = 0;
  arg_0.t = 0;

  //Declare args for source at z.
  QPropWArg arg_z;
  setup_qpropwarg_cg(arg_z);

  // Propagator calculation objects and memory allocation
  //
  // Using x[4] = X(x,y,z,t)
  //       y[4] = Y(x,y,z,t)
  //       z[4] = Z(x,y,z,t)
  int x[4];
  int y[4];
  int z[4];
  int x_idx4d, x_idx3d, y_idx4d, y_idx3d, z_idx4d, z_idx3d;
  int vol4d = GJP.XnodeSites()*GJP.YnodeSites()*GJP.ZnodeSites()*GJP.TnodeSites();
  int vol3d = GJP.XnodeSites()*GJP.YnodeSites()*GJP.ZnodeSites();
  int xnodes = GJP.XnodeSites();
  int ynodes = GJP.YnodeSites();
  int znodes = GJP.ZnodeSites();
  double norm = pow(vol3d, -0.5);
  
  int max_mom = NSITES_3D;
  mom3D mom(max_mom, SINPz_Pz/(1.0*100));

  int s1 = 0;
  int c1 = 0;
  int s2 = 0;
  int c2 = 0;
  int sc_idx = 0;

  //use t to represent the time slice.
  //int t = 0;

  //In these arrays, we will use the index convention [sink_index + vol3d*source_index]
  WilsonMatrix *t3_arr = (WilsonMatrix*)smalloc(vol3d*vol3d*sizeof(WilsonMatrix));
  WilsonMatrix *t2_arr = (WilsonMatrix*)smalloc(vol3d*vol3d*sizeof(WilsonMatrix));
  //Initialise
  for (int i=0; i<vol3d*vol3d; i++) {
    t3_arr[i]    *= 0.0;
    t2_arr[i]    *= 0.0;
  }

  //Arrays to store the trace data
  fftw_complex *FT_t4 = (fftw_complex*)smalloc(vol3d*sizeof(fftw_complex));
  fftw_complex *FT_t2 = (fftw_complex*)smalloc(vol3d*vol3d*sizeof(fftw_complex));
  fftw_complex *FT_t3 = (fftw_complex*)smalloc(vol3d*vol3d*sizeof(fftw_complex));
  
  //Use this array several times for 9d D0, D1, D2.
  fftw_complex *FT_9d  = (fftw_complex*)smalloc(vol3d*vol3d*vol3d*sizeof(fftw_complex));
  
  //Momentum source array.
  fftw_complex *FFTW_mom_arr  = (fftw_complex*)smalloc(vol3d*sizeof(fftw_complex));
  //Initialise
  for (int i=0; i<vol3d*vol3d*vol3d; i++) {
    for(int a=0; a<2; a++){
      FT_9d[i][a]  = 0.0;    
      if(i<vol3d*vol3d) {
	FT_t3[i][a]  = 0.0;
	FT_t2[i][a]  = 0.0;
      }
      if(i<vol3d) {
	FT_t4[i][a]  = 0.0;
	FFTW_mom_arr[i][a]  = 0.0;
      }
    }
  }
 //gaahhbage 
  FFT_F(9, NSITES_3D, FT_9d);
  FFT_B(9, NSITES_3D, FT_9d);

  FFT_F(6, NSITES_3D, FT_t2);
  FFT_B(6, NSITES_3D, FT_t2);
 
  FFT_F(3, NSITES_3D, FFTW_mom_arr);
  FFT_B(3, NSITES_3D, FFTW_mom_arr); 

  WilsonMatrix t1;
  WilsonMatrix t1c;
  WilsonMatrix t4;
  WilsonMatrix t4c;
  WilsonMatrix t4t1c;
  WilsonMatrix t2t3c;
  WilsonMatrix t3;
  WilsonMatrix t3c;
  WilsonMatrix t2;
  WilsonMatrix t2c;
		
  //Rcomplex mom_src;
  //WilsonMatrix temp;

  Rcomplex t1t1c_tr;
  Rcomplex t4t4c_tr;
  Rcomplex d2_tr;
  Rcomplex t2t2c_tr;
  Rcomplex t3t3c_tr;

  //////////////////////
  // Start simulation //
  ////////////////////// 

  Float *time = (Float*)smalloc(10*sizeof(Float));
  for(int a=0; a<10; a++) time[a] = 0.0;

  char lattice[256];
  
  while (sweep_counter < total_updates) {
    for (int n = 1; n <= NSKIP; n++) {
#ifndef READ
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
    
    if (sweep_counter == NTHERM) {
      printf("thermalization complete. \n");
    }
    if (sweep_counter >= NTHERM) {
      // Use this code to specify a gauge configuration.
      #ifdef QUENCH
        sprintf(lattice, LATT_PATH"QU/lat_hb_B%.2f_%d-%d_%d.dat", BETA, NSITES_3D, NSITES_T, sweep_counter);
      #else
        sprintf(lattice, LATT_PATH"UNQ/lat_hmc_B%.2f_M%.3f_%d-%d_%d.dat", BETA, MASS, NSITES_3D, NSITES_T, sweep_counter);
      #endif
#ifdef READ
      ReadLatticeParallel(lat,lattice);
#else
      WriteLatticeParallel(lat,lattice);
#endif
      
      gaugecounter = 1;
      
      // We will compute two arrays of momentum source propagators.
      // One array is of t2 S(x,z)
      // One array is of t3 S(y,z)
      // Each array will be indexed arr[sink_index + vol*source_index].
      
      // The sources for these arrays are calculated using the backaward FT of momentum states.
      // E.G., momemtum state P_0=(0,0,0) is used to calculated the position space state
      // X_0[n] = \frac{1}{sqrt(V)} * \sum_{m} e^{(-2i*pi/N)*n*m} * P_0[m].
      // This source is then used in the inversion to calculate an propagator M_0. M_0 <P_0|  has,
      // strong overlap with the P_0 state. This is repeated for small momenta (e.g. |P| < 1) and the propagators
      // from each inversion are summed and normalised by the number of momenta used k:
      // M = 1/sqrt(k) sum_k M_k <P_k|  The resulting propagator M has strong overlap with the low momentum states.
      // N.B. One can show that using all possible momenta K, the full propagator matrix is recovered.
      
      // The 0-mom source at the origin is calculated outside the time loop.
      int P0[3] = {0,0,0};
      
      arg_0.t = 0;
      QPropWMomSrcSmeared qprop_0(lat, &arg_0, P0, &g_arg_mom, &c_arg);
      qprop_0.GaussSmearSinkProp(g_arg_mom);
      cout<<"Sink Smear 0 complete."<<endl;
      
      //////////////////////////////////
      // Begin loop over time slices. //
      //////////////////////////////////

      for (int t=0; t<GJP.TnodeSites(); t++) {
	//Reinitialise all propagator arrays.
        for (int i=0; i<vol3d*vol3d; i++) {
	  t2_arr[i]    *= 0.0;
	  t3_arr[i]    *= 0.0;
        }
	
	stopwatchStart();
	
	//Generate momentum source
	int n_mom_srcs    = 0;

	for (mom.P[2] = 0; mom.P[2] < max_mom; mom.P[2]++)
	  for (mom.P[1] = 0; mom.P[1] < max_mom; mom.P[1]++)
	    for (mom.P[0] = 0; mom.P[0] < max_mom; mom.P[0]++) {
	      
	      cout<<"MOM = "<<mom.P[0]<<" "<<mom.P[1]<<" "<<mom.P[2]<<endl;
	      
	      cout<<"NORM_MOM_SZE = "<<mom.mod()/M_PI<<endl;
	      //frac = sin(p)/p
	      Float frac = sin(mom.mod())/(mom.mod());
	      cout<<"SIN(Pz)/Pz = "<<frac<<endl;
	      
	      if(frac > mom.sin_cutoff || (mom.P[0] == 0 && mom.P[1] == 0 && mom.P[2] == 0) ){
		
		//Set momentum
		int P[3] = {mom.P[0], mom.P[1], mom.P[2]};
		
		// The CPS momentum source function uses an unnormalised
		// source, so we take the product of both normalisation
		// factors and place them here on the FFTW_mom_arr.
		// A further normalisation to perform comes from the number n_mom_srcs
		// of momentum sources. This is done later in when the trace of
		// of the propagators is caculated.
		
		
		//Get Momentum Propagator
		arg_z.t = t;
		//QPropWMomSrc qprop_mom(lat, &arg_z, P, &c_arg);
		QPropWMomSrcSmeared qprop_mom(lat, &arg_z, P, &g_arg_mom, &c_arg);
		cout<<"Inversion "<<(n_mom_srcs+1)<<" complete."<<endl;
		qprop_mom.GaussSmearSinkProp(g_arg_mom);
		cout<<"Sink Smear "<<(n_mom_srcs+1)<<" complete."<<endl;
		
		int z_idx4d, z_idx3d, x_idx4d, x_idx3d, y_idx4d, y_idx3d;
		//Loop over sources at z.
		z[3] = t;
		for (z[2]=0; z[2]<znodes; z[2]++)
		  for (z[1]=0; z[1]<ynodes; z[1]++)
		    for (z[0]=0; z[0]<xnodes; z[0]++) {
		      z_idx4d = lat.GsiteOffset(z)/4;
		      z_idx3d = z_idx4d - vol3d*z[3];
		      
		      cout<<"mom_src "<<qprop_mom.mom_src(z_idx4d)<<endl;
		      
		      //Loop over sinks at x.
		      x[3] = 0;
		      for (x[2]=0; x[2]<znodes; x[2]++)
			for (x[1]=0; x[1]<ynodes; x[1]++)
			  for (x[0]=0; x[0]<xnodes; x[0]++) {
			    x_idx4d = lat.GsiteOffset(x)/4;
			    x_idx3d = x_idx4d - vol3d*x[3];
			    
			    //Build t2 array.
			    t2_arr[x_idx3d + vol3d*z_idx3d] += qprop_mom[x_idx4d]*conj(qprop_mom.mom_src(z_idx4d));
			  }
		      
		      //Loop over sinks at y.
		      y[3] = t;
		      for (y[2]=0; y[2]<znodes; y[2]++)
			for (y[1]=0; y[1]<ynodes; y[1]++)
			  for (y[0]=0; y[0]<xnodes; y[0]++) {
			    y_idx4d = lat.GsiteOffset(y)/4;
			    y_idx3d = y_idx4d - vol3d*y[3];
			    
			    //Build t3 array.
			    t3_arr[y_idx3d + vol3d*z_idx3d] += qprop_mom[y_idx4d]*conj(qprop_mom.mom_src(z_idx4d));
			  }
		    }
		n_mom_srcs++; 
		cout << "momentum sources: "<<1+mom.P[2]*max_mom*max_mom + mom.P[1]*max_mom + mom.P[0]<<" / "<<pow(max_mom,3)<<" checked"<<endl;
	      }
	    }
	
	cout<<"FLAG 1"<<endl;
	//inversions + fill      
	time[1] = stopwatchReadSeconds();
	stopwatchStart();
	
	//////////////////////////////////////////////////////////////////
	// End momentum source propagator calculation for time slice t. //
	//////////////////////////////////////////////////////////////////
	
	
	///////////////////////////////////////////////
	// Begin summation of trace at time slice t. //
	///////////////////////////////////////////////
	      
	// The t1, t1c, t4, and t4c propagators are calculated 'on the fly'
	// within the trace summation.
      
	//Reinitialise all trace variables
	
	t1  *= 0.0;
	t1c *= 0.0;
	t2  *= 0.0;
	t2c *= 0.0;
	t3  *= 0.0;
	t3c *= 0.0;
	t4  *= 0.0;
	t4c *= 0.0;
	t4t1c *= 0.0;
	t2t3c *= 0.0;
      
	t1t1c_tr *= 0.0;
	t2t2c_tr *= 0.0;
	t3t3c_tr *= 0.0;
	t4t4c_tr *= 0.0;
	d2_tr *= 0.0;
            
	for (int i=0; i<vol3d*vol3d*vol3d; i++) 
	  for(int a=0; a<2; a++) {
	    FT_9d[i][a] = 0.0;
	    if(i<vol3d*vol3d) {
	      FT_t3[i][a] = 0.0;
	      FT_t2[i][a] = 0.0;
	    }
	    if(i<vol3d) {
	      FT_t4[i][a] = 0.0;
	    }
	  }
	//Sum over X
	x[3] = 0;
	for (x[2]=0; x[2]<znodes; x[2]++)
	  for (x[1]=0; x[1]<ynodes; x[1]++)
	    for (x[0]=0; x[0]<xnodes; x[0]++) {
	      x_idx4d = lat.GsiteOffset(x)/4;
	      x_idx3d = x_idx4d - vol3d*x[3];
	      
	      t1 = qprop_0[x_idx4d];
	      t1c = t1.conj_cp();
	      
	      //Sum over Y
	      y[3] = t;
	      for (y[2]=0; y[2]<znodes; y[2]++)
		for (y[1]=0; y[1]<ynodes; y[1]++)
		  for (y[0]=0; y[0]<xnodes; y[0]++) {
		    y_idx4d = lat.GsiteOffset(y)/4;
		    y_idx3d = y_idx4d - vol3d*y[3];
		    
		    t4 = qprop_0[y_idx4d];
		  
		    // Use this condition so that t4t4c is calculated only once
		    // over X per time slice.
		    if (x_idx3d == 0) {
		      //Perform t4t4c trace sum for D0 graph.
		      FT_t4[y_idx3d][0] = MMDag_re_tr(t4);
		      FT_t4[y_idx3d][1] = 0.0;
		    }
		    
		    //Declare new Wilson Matrix t4*t1c for D2 and compute
		    t4t1c = t4;
		    t4t1c *= t1c;
		    
		    //Sum over Z.
		    z[3] = t;
		    for (z[2]=0; z[2]<znodes; z[2]++)
		      for (z[1]=0; z[1]<ynodes; z[1]++)
			for (z[0]=0; z[0]<xnodes; z[0]++) {
			  z_idx4d = lat.GsiteOffset(z)/4;
			  z_idx3d = z_idx4d - vol3d*z[3];
			  
			  //Declare new Wilson Matrix t2*t3c and compute it.			  
			  t2t3c = t2_arr[x_idx3d + vol3d*z_idx3d];
			  t3c   = t3_arr[y_idx3d + vol3d*z_idx3d].conj_cp();
			  t2t3c *= t3c;
			  
			  //Perform t4t1c * t2t3c trace sum for D2 graph.
			  d2_tr = Trace(t4t1c, t2t3c);
			  
			  //Create 9d array for D2.			  
			  FT_9d[x_idx3d + vol3d*(y_idx3d + vol3d*z_idx3d)][0] = d2_tr.real();
			  FT_9d[x_idx3d + vol3d*(y_idx3d + vol3d*z_idx3d)][1] = d2_tr.imag();
			  
			  ///////////////////////////////////////////////////////////////////
			  // Use this condition so that t2t2c is calculated only over
			  // x1 and x3 loops per time slice. 
			  if (y_idx3d == 0) {
			    //Retrieve propagators for t2t2c trace sum.
			    FT_t2[x_idx3d + vol3d*z_idx3d][0] = MMDag_re_tr(t2_arr[x_idx3d + vol3d*z_idx3d]);
			    FT_t2[x_idx3d + vol3d*z_idx3d][1] = 0.0;
			  }
			  // Use this condition so that t3t3c is calculated only over
			  // x2 and x3 loops per time slice. 
			  if (x_idx3d == 0) {
			    
			    //Retrieve propagators for t3t3c trace sum.
			    FT_t3[y_idx3d + vol3d*z_idx3d][0] = MMDag_re_tr(t3_arr[y_idx3d + vol3d*z_idx3d]);
			    FT_t3[y_idx3d + vol3d*z_idx3d][1] = 0.0;
			  }
			  ///////////////////////////////////////////////////////////////////
			}
		  }
	    }
	
	//Fill the trace arrays
	time[2] = stopwatchReadSeconds();
	
	cout<<"FLAG 3"<<endl;

	///////////////////////////////////////////////
	// Write traces to file for post-processing. //
	///////////////////////////////////////////////
      
	char file[256];
  	FFT_F(6, NSITES_3D, FT_t2);
	FFT_F(6, NSITES_3D, FT_t3);
	FFT_F(3, NSITES_3D, FT_t4);

	// if(t==0) {    
	// sprintf(file, "%d-%d_3-0.1_msmsFT_6d_data/t1t1c_TR_%d_%d-%d_%d_%d.dat",  NSITES_3D, NSITES_T, n_mom_srcs, NSITES_3D, NSITES_T, sweep_counter, t);	  
	// FILE *qt1tr   = Fopen(file, "a");
	// for(int snk =0; snk<vol3d; snk++) {
	// Fprintf(qt1tr,  "%d %d %d %.16e %.16e\n", sweep_counter, t, snk, FT_t4[snk][0], FT_t4[snk][1]);
	// }
	// Fclose(qt1tr);
	// }
	
	sprintf(file, DATAPATH"t4t4c_TR_%d_%d-%d_%d_%d.dat",  n_mom_srcs, NSITES_3D, NSITES_T, sweep_counter, t);
	FILE *qt4tr   = Fopen(file, "a");
	for(int snk =0; snk<vol3d; snk++) {
	  Fprintf(qt4tr,  "%d %d %d %.16e %.16e\n", sweep_counter, t, snk, FT_t4[snk][0], FT_t4[snk][1]);
	}
	
	sprintf(file, DATAPATH"t2t2c_TR_%d_%d-%d_%d_%d.dat",  n_mom_srcs, NSITES_3D, NSITES_T, sweep_counter, t);
	FILE *qt2tr   = Fopen(file, "a");
	sprintf(file, DATAPATH"t3t3c_TR_%d_%d-%d_%d_%d.dat",  n_mom_srcs, NSITES_3D, NSITES_T, sweep_counter, t);
	FILE *qt3tr   = Fopen(file, "a");
	
	for(int src =0; src<vol3d; src++) {
	  for(int snk =0; snk<vol3d; snk++) {
	    Fprintf(qt2tr,"%d %d %d %d %.16e %.16e\n", sweep_counter, t, src, snk, FT_t2[snk + vol3d*src][0], FT_t2[snk + vol3d*src][1]);
	    Fprintf(qt3tr,"%d %d %d %d %.16e %.16e\n", sweep_counter, t, src, snk, FT_t3[snk + vol3d*src][0], FT_t3[snk + vol3d*src][1]);
	  }
	}
	

	Fclose(qt2tr);
	Fclose(qt3tr);
	Fclose(qt4tr);

	//////////////////////////
	// FFT the 9D D2 array. //
	//////////////////////////
      
	stopwatchStart();
	
	FFT_F(9, NSITES_3D, FT_9d);
	//time for D2 6d FFT
	time[4] = stopwatchReadSeconds();      
	//wtf == 'write to file', include/FFTW_functions.cpp
	FFT_wtf_ZYX(FT_9d, 2, SINPz_Pz, SINPxy_Pxy, n_mom_srcs, NSITES_3D, NSITES_T, sweep_counter, t);
	
	//sprintf(file, "T_data/times_%d-%d_%d_%d.dat", NSITES_3D, NSITES_T, sweep_counter, t);
	//FILE *time_fp = Fopen(file, "a");
	//Fprintf(time_fp, "%.4f %.4f %.4f %.4f\n", time[1], time[2], time[3], time[4]);
	//Fclose(time_fp); 
	
	//////////////////////////////////////////
	// End trace summation at time slice t. //
	//////////////////////////////////////////
      }
    }
  }
  ////////////////////
  // End simulation //
  ////////////////////
  
  sfree(t2_arr);
  sfree(t3_arr);

  //sfree(FT_t1);
  sfree(FT_t4);
  sfree(FT_t2);
  sfree(FT_t3);
  sfree(FT_9d);  


  sfree(time);

  //End();
  return 0;
}
