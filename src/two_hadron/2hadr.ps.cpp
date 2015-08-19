//here
#include "params.h"
#include "args.h"

USING_NAMESPACE_CPS


// We will compute the trace of all graphs
// over all lattice sites
// using a point sources at 0 and Z

#define NMOM 1

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
  QPropWGaussArg g_arg;
  setup_g_arg(g_arg);

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
  int vol4d = GJP.XnodeSites()*GJP.YnodeSites()*GJP.ZnodeSites()*GJP.TnodeSites();
  int vol3d = GJP.XnodeSites()*GJP.YnodeSites()*GJP.ZnodeSites();
  int x[4];
  int y[4];
  int z[4];
  
  int x_idx4d = 0;
  int x_idx3d = 0;
  int y_idx4d = 0;
  int y_idx3d = 0;
  int z_idx4d = 0;
  int z_idx3d = 0;

  int s1 = 0;
  int c1 = 0;
  int s2 = 0;
  int c2 = 0;
  int sc_idx = 0;

  //use t to represent the time slice.
  int t = 0;

  //Arrays to store the trace data
  fftw_complex *FT_t1  = (fftw_complex*)smalloc(vol3d*sizeof(fftw_complex));
  fftw_complex *FT_t4  = (fftw_complex*)smalloc(vol3d*sizeof(fftw_complex));
  fftw_complex *FT_t2  = (fftw_complex*)smalloc(vol3d*vol3d*sizeof(fftw_complex));
  fftw_complex *FT_t3  = (fftw_complex*)smalloc(vol3d*vol3d*sizeof(fftw_complex));
  
  fftw_complex *FT_d2_9d  = (fftw_complex*)smalloc(vol3d*vol3d*vol3d*sizeof(fftw_complex));
  
  WilsonMatrix *t3_arr = (WilsonMatrix*)smalloc(vol3d*vol3d*sizeof(WilsonMatrix));
  WilsonMatrix *t2_arr = (WilsonMatrix*)smalloc(vol3d*vol3d*sizeof(WilsonMatrix));
  //Initialise
  for (int i=0; i<vol3d*vol3d; i++) {
    t3_arr[i]    *= 0.0;
    t2_arr[i]    *= 0.0;
  }


  //Initialise
  for (int i=0; i<vol3d*vol3d*vol3d; i++) {
    for(int a=0; a<2; a++){
      FT_d2_9d[i][a]  = 0.0;
      if(i<vol3d*vol3d) {
	FT_t2[i][a]  = 0.0;
	FT_t3[i][a]  = 0.0;
      }
      if(i<vol3d) {
	FT_t1[i][a]  = 0.0;
	FT_t4[i][a]  = 0.0;
      }
    }
  }
  
  FFT_F(3, NSITES_3D, FT_t4);
  FFT_F(6, NSITES_3D, FT_t3);
  FFT_F(9, NSITES_3D, FT_d2_9d);
//end initialising fftw-- first transform is always garbage, so get it out of the way with useless data
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
		
  Rcomplex d2_tr;

  double time[5]; //hold stopwatch values
  char lattice[256]; //lattice config file
  char file[256]; //output file

  //////////////////////
  // Start simulation //
  ////////////////////// 

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

    if (sweep_counter == NTHERM) printf("thermalization complete. \n");
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
      // Use this code to write a gauge configuration.
      //char lattice[256];
      //sprintf(lattice,"latGPUWIL4-8.dta");
      //WriteLatticeParallel(lat,lattice);
      
      gaugecounter = 1;
      
      //Reinitialise corr[t]
      //for (int a=0; a < GJP.TnodeSites(); a++) corr[a] = 0.0;

      // We will compute two arrays of point source propagators.
      // One array is of t2
      // One array is of t3
      // Each array will be indexed arr[sink_index + vol*source_index]
      // We will the sum over all spatial points and calculate the 4 traces of the 4 graphs,
      // and use \gamma_5 conjugation to obtain down type quark propagators.
      // N.B. for the fermion loops, we must multiply by the trace by -1 in accordance with the minus
      // sign each picks up from Wick contractions of the quark operators.


      // The point source is calculated outside the time loop.
      //QPropWPointSrc qprop0(lat, &arg0, &c_arg);

      QPropWGaussSrc qprop_0(lat, &arg_0, &g_arg, &c_arg);
      qprop_0.GaussSmearSinkProp(g_arg);

      //////////////////////////////////
      // Begin loop over time slices. //
      //////////////////////////////////

      for (t=0; t<GJP.TnodeSites(); t++) {
	//for (t=t_in; t<t_in+1; t++) {
	//for (t=t_in; t<GJP.TnodeSites(); t++) {
	//for (t=t_in; t<t_in+2; t++) {
	
	//Reinitialise all propagator arrays.
	for (int i=0; i<vol3d*vol3d; i++) {
	  t2_arr[i]    *= 0.0;
	  t3_arr[i]    *= 0.0;
	}

	//Loop over sources at z.
	z[3] = t;
	for (z[2]=0; z[2]<GJP.ZnodeSites(); z[2]++)
	  for (z[1]=0; z[1]<GJP.YnodeSites(); z[1]++)
	    for (z[0]=0; z[0]<GJP.XnodeSites(); z[0]++) {
	      z_idx4d = lat.GsiteOffset(z)/4;
	      z_idx3d = z_idx4d - vol3d*z[3];

	      arg_z.t=z[3];
	      arg_z.z=z[2];
	      arg_z.y=z[1];
	      arg_z.x=z[0];

	      QPropWGaussSrc qprop_z(lat, &arg_z, &g_arg, &c_arg);
	      qprop_z.GaussSmearSinkProp(g_arg);
	      //QPropWPointSrc qpropx3(lat, &argz, &c_arg);

	      //Loop over sinks at x.
	      x[3] = 0;
              for (x[2]=0; x[2]<GJP.ZnodeSites(); x[2]++)
                for (x[1]=0; x[1]<GJP.YnodeSites(); x[1]++)
                  for (x[0]=0; x[0]<GJP.XnodeSites(); x[0]++) {
                    x_idx4d = lat.GsiteOffset(x)/4;
		    x_idx3d = x_idx4d - vol3d*x[3];

		    //Build t2 array.
		    t2_arr[x_idx3d + vol3d*z_idx3d] = qprop_z[x_idx4d];
		  }

	      //Loop over sinks at y.
	      y[3] = t;
	      for (y[2]=0; y[2]<GJP.ZnodeSites(); y[2]++)
		for (y[1]=0; y[1]<GJP.YnodeSites(); y[1]++)
		  for (y[0]=0; y[0]<GJP.XnodeSites(); y[0]++) {
		    y_idx4d = lat.GsiteOffset(y)/4;
		    y_idx3d = y_idx4d - vol3d*y[3];
		    
		    //Build t3 array
		    t3_arr[y_idx3d + vol3d*z_idx3d] = qprop_z[y_idx4d];
		  }
	    }
	
	
	///////////////////////////////////////////////////////////////
	// End point source propagator calculation for time slice t. //
	///////////////////////////////////////////////////////////////


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
      
	d2_tr *= 0.0;
	
	for (int i=0; i<vol3d*vol3d*vol3d; i++) 
	  for(int a=0; a<2; a++) {
	    FT_d2_9d[i][a] = 0.0;
	    if(i<vol3d*vol3d) {
	      FT_t3[i][a] = 0.0;
	      FT_t2[i][a] = 0.0;
	      if(i<vol3d) {
		FT_t1[i][a] = 0.0;
		FT_t4[i][a] = 0.0;
	      }
	    }
	  }
	
	//Sum over X
	x[3] = 0;
	for (x[2]=0; x[2]<GJP.ZnodeSites(); x[2]++)
	  for (x[1]=0; x[1]<GJP.YnodeSites(); x[1]++)
	    for (x[0]=0; x[0]<GJP.XnodeSites(); x[0]++) {
	      x_idx4d = lat.GsiteOffset(x)/4;
	      x_idx3d = x_idx4d - vol3d*x[3];
	      
	      t1 = qprop_0[x_idx4d];
	      t1c = t1.conj_cp();
	      
	      //Perform t1t1c trace sum for D1 graph.
	      //t1t1c_tr = Trace(t1, t1c);
	      //d1_t1t1c_re_tr += t1t1c_tr.real();
	      
	      //Sum over Y
	      y[3] = t;
	      for (y[2]=0; y[2]<GJP.ZnodeSites(); y[2]++)
		for (y[1]=0; y[1]<GJP.YnodeSites(); y[1]++)
		  for (y[0]=0; y[0]<GJP.XnodeSites(); y[0]++) {
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
		    
		    //Sum over z
		    z[3] = t;
		    for (z[2]=0; z[2]<GJP.ZnodeSites(); z[2]++)
		      for (z[1]=0; z[1]<GJP.YnodeSites(); z[1]++)
			for (z[0]=0; z[0]<GJP.XnodeSites(); z[0]++) {
			  z_idx4d = lat.GsiteOffset(z)/4;
			  z_idx3d = z_idx4d - vol3d*z[3];
			  
			  //Declare new Wilson Matrix t2*t3c and compute it.			  
			  t2t3c = t2_arr[x_idx3d + vol3d*z_idx3d];
			  t3c   = t3_arr[y_idx3d + vol3d*z_idx3d].conj_cp();
			  t2t3c *= t3c;
			  
			  //Perform t4t1c * t2t3c trace sum for D2 graph.
			//instead of summing Re Tr
			  d2_tr = Trace(t4t1c, t2t3c);
			  FT_d2_9d[vol3d*vol3d*z_idx3d + vol3d*y_idx3d + x_idx3d][0] = d2_tr.real();
			  FT_d2_9d[vol3d*vol3d*z_idx3d + vol3d*y_idx3d + x_idx3d][1] = d2_tr.imag();
			  
			  ////////////////////////////////////////////////////////////
			  
			  // Use this condition so that t2t2c is calculated only over
			  // x1 and x3 loops per time slice. 
			  if (y_idx3d == 0) {
			    //Retrieve propagators for t2t2c trace sum.
			    FT_t2[vol3d*z_idx3d + x_idx3d][0] = MMDag_re_tr(t2_arr[vol3d*z_idx3d + x_idx3d]);
			    FT_t2[vol3d*z_idx3d + x_idx3d][1] = 0.0;
			  }
			   
			   ////////////////////////////////////////////////////////////
			   
			   // Use this condition so that t3t3c is calculated only over
			   // x2 and x3 loops per time slice. 
			   if (x_idx3d == 0) {
			     
			     //Retrieve propagators for t3t3c trace sum.
			     FT_t3[vol3d*z_idx3d + y_idx3d][0] = MMDag_re_tr(t3_arr[vol3d*z_idx3d+ y_idx3d]);
			     FT_t3[vol3d*z_idx3d + y_idx3d][1] = 0.0;
			   }
			}
		  }
	    }
	
	//corr[t] += ((d0_t4t4c_re_tr * d0_t2t2c_re_tr) + (d1_t1t1c_re_tr * d1_t3t3c_re_tr));
	//corr[t] -= ((d2_re_tr) + (d3_re_tr));

	////////////////////////////////////
	// Write traces of loops to file. //
	////////////////////////////////////

//fft everything, not just 0-mom mode
  	FFT_F(6, NSITES_3D, FT_t2);
	FFT_F(6, NSITES_3D, FT_t3);
	FFT_F(3, NSITES_3D, FT_t4);

	if(t==0) {
	  sprintf(file, DATAPATH"t1t1c_FT_%d_%d-%d_%d_%d.dat", NMOM, NSITES_3D, NSITES_T, sweep_counter, t);	  
	  FILE *qt1tr   = Fopen(file, "a");
	  for(int snk =0; snk<vol3d; snk++) {
	    Fprintf(qt1tr,  "%d %d %d %.16e %.16e\n", sweep_counter, t, snk, FT_t4[snk][0], FT_t4[snk][1]);
	  }
	  Fclose(qt1tr);
	}
	
	sprintf(file, DATAPATH"t4t4c_FT_%d_%d-%d_%d_%d.dat", NMOM, NSITES_3D, NSITES_T, sweep_counter, t);
	FILE *qt4tr   = Fopen(file, "a");
	for(int snk =0; snk<vol3d; snk++) {
	  Fprintf(qt4tr,  "%d %d %d %.16e %.16e\n", sweep_counter, t, snk, FT_t4[snk][0], FT_t4[snk][1]);
	}
	
	sprintf(file, DATAPATH"t2t2c_FT_%d_%d-%d_%d_%d.dat", NMOM, NSITES_3D, NSITES_T, sweep_counter, t);
	FILE *qt2tr   = Fopen(file, "a");
	sprintf(file, DATAPATH"t3t3c_FT_%d_%d-%d_%d_%d.dat", NMOM, NSITES_3D, NSITES_T, sweep_counter, t);
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

	////////////////////////
	// FFT the D2 arrays. //
	////////////////////////
      
	stopwatchStart();
	
	FFT_F(9, NSITES_3D, FT_d2_9d);
	//time for D2 6d FFT
	time[4] = stopwatchReadSeconds();	
	FFT_wtf_ZYX(FT_d2_9d, 2, SINPz_Pz, SINPxy_Pxy, NMOM, NSITES_3D, NSITES_T, sweep_counter, t);	
      }
      
      /////////////////////////
      // End loop over time. //
      /////////////////////////
      
    }
  }
  
  ////////////////////
  // End simulation //
  ////////////////////

  sfree(t2_arr);
  sfree(t3_arr);

  End();
  return 0;
}
