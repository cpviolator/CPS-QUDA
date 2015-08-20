//here
#include "params.h"
#include "args.h"

USING_NAMESPACE_CPS

// We will compute the trace of all graphs
// over all lattice sites
// using a point sources at 0 and Z


int main(int argc, char *argv[]) {
  
  int gaugecounter = 1;
  
  Start(&argc,&argv);
  int seed = atoi(argv[1]);
  int t = atoi(argv[2]); //still called argv[3] in comments..
   
  DoArg do_arg;
  setup_do_arg(do_arg, seed); 
  GJP.Initialize(do_arg);  
  
  GwilsonFclover lat;
  CommonArg c_arg;

  //Declare args for Gaussian Smearing
  QPropWGaussArg g_arg;
  setup_g_arg(g_arg);

  QPropWRandArg r_arg;
  r_arg.rng = SRC_TYPE; 
  
  int sweep_counter = 0;
  int total_updates = NTHERM + NSKIP*(NDATA-1);

  string is_qu;
  #ifdef QUENCH
  GhbArg ghb_arg;
  ghb_arg.num_iter = 1;
  AlgGheatBath hb(lat, &c_arg, &ghb_arg);
  is_qu = "QU";
  #else
  HmdArg hmd_arg; 
  setup_hmd_arg(hmd_arg);
  AlgHmcPhi hmc(lat, &c_arg, &hmd_arg); 
  is_qu = "UNQ";
  #endif
  
  QPropWArg arg_ps;
  setup_qpropwarg_cg(arg_ps);
  arg_ps.t=0;
  arg_ps.x=0;
  arg_ps.y=0;
  arg_ps.z=0;
  
  //Declare args for source at x3.
  QPropWArg arg_sto;
  setup_qpropwarg_cg(arg_sto);
  
  // Propagator calculation objects and memory allocation
  //
  // Using x1[4] = X(x,y,z,t)
  //       x2[4] = Y(x,y,z,t)
  //       x3[4] = Z(x,y,z,t)

  int vol4d = GJP.XnodeSites()*GJP.YnodeSites()*GJP.ZnodeSites()*GJP.TnodeSites();
  int vol3d = GJP.XnodeSites()*GJP.YnodeSites()*GJP.ZnodeSites();
  int xnodes = GJP.XnodeSites();
  int ynodes = GJP.YnodeSites();
  int znodes = GJP.ZnodeSites();
  int x[4];
  int y[4];
  int z[4];
  int x_idx3d = 0;
  int x_idx4d = 0;
  int y_idx3d = 0;
  int y_idx4d = 0;
  int z_idx3d = 0;
  int z_idx4d = 0;
  int z_arr_idx = 0;

  //Dilution variables.
  int L[3]={2, 2, 2};
  int index=0;
  int subvol=L[0]*L[1]*L[2];

  //All >8-byte variables used inside simulation are declared here
  //and grouped into vectors for easy re-initialization
  //In these arrays, we will use the index convention [sink_index + vol3d*source_index]
  
  WilsonMatrix* t2A_full_arr = new WilsonMatrix[vol3d*vol3d];
  WilsonMatrix* t3A_full_arr = new WilsonMatrix[vol3d*vol3d];
  WilsonMatrix* t2B_full_arr = new WilsonMatrix[vol3d*vol3d];
  WilsonMatrix* t3B_full_arr = new WilsonMatrix[vol3d*vol3d];
  vector<WilsonMatrix*> allWMarrs = {t2A_full_arr, t3A_full_arr, t2B_full_arr, t3B_full_arr};
  //Initialise
  for (int r=0; r<allWMarrs.size(); r++)
    for (int v=0; v<vol3d*vol3d; v++) 
      allWMarrs[r][v] *= 0.0;
  
  WilsonMatrix t1, t1c, t4, t4c, t4t1c, t1t4c, t2t3c, t3t2c, t2, t2c, t2A, t2Ac, t2B, t2Bc, t3, t3c, t3A, t3Ac, t3B, t3Bc;
  vector<WilsonMatrix*> allWMs = {&t1,&t1c,&t4,&t4c,&t4t1c,&t1t4c,&t2t3c,&t3t2c,&t2,&t2c,&t2A,&t2Ac,&t2B,&t2Bc,&t3,&t3c,&t3A,&t3Ac,&t3B,&t3Bc};
  for (int r=0; r<allWMs.size(); r++)
    *(allWMs[r]) *= 0.0;

  Rcomplex tr;

  Float t1t1c_re_tr, t4t4c_re_tr, /*t3t3c_1_re_tr, t3t3c_1_re_tr, d2_1_re_tr, d3_1_re_tr, */t2t2c_2_re_tr, t3t3c_2_re_tr, d2_2_re_tr, d3_2_re_tr;
  vector<Float*> allReTrs = {&t1t1c_re_tr,&t4t4c_re_tr,/*&t3t3c_1_re_tr,&t3t3c_1_re_tr,&d2_1_re_tr,&d3_1_re_tr,*/&t2t2c_2_re_tr,&t3t3c_2_re_tr,&d2_2_re_tr,&d3_2_re_tr};
  for (int r=0; r<allReTrs.size(); r++)
    *(allReTrs[r]) = 0.0;
  //for writing to file
  vector<string> allReTrNames = {"t1t1c_re_tr","t4t4c_re_tr",/*&t3t3c_1_re_tr",",&t3t3c_1_re_tr",",&d2_1_re_tr",",&d3_1_re_tr",*/"t2t2c_2_re_tr","t3t3c_2_re_tr","d2_2_re_tr","d3_2_re_tr"};

  double time[5]; //hold stopwatch values
  char lattice[256]; //lattice config file
  char file[256]; //output file
  
  //////////////////////
  // Start simulation //
  ////////////////////// 

  stopwatchStart();

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
      time[0] = stopwatchReadSeconds();
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
      
      // The point source is calculated outside the time loop.
      QPropWGaussSrc qprop_point(lat, &arg_ps, &g_arg, &c_arg);
      qprop_point.GaussSmearSinkProp(g_arg);
      
      //////////////////////////////////
      // Begin loop over time slices. //
      //////////////////////////////////
      
      //for (t=0; t<GJP.TnodeSites(); t++) {
	
	//Reinitialise all propagator arrays.
	for (int r=0; r<allWMarrs.size(); r++)
	  for (int v=0; v<vol3d*vol3d; v++)
	    allWMarrs[r][v] *= 0.0;
	
	/////////////////////////
	// Compute t2A and t3A //
	/////////////////////////
	time[1] = stopwatchReadSeconds();
	
	for(int j=0; j<NRAND; j++) {
//THIS IS ARGV[2]
	  arg_sto.t = t;	  
	  for (index=0;index<subvol;index++) {
	    QPropWRandWallSrcDilutedSmeared qprop_stoch(lat, &arg_sto, &r_arg, &c_arg, &g_arg, L, index);
	    qprop_stoch.GaussSmearSinkProp(g_arg);
	    //Loop over sources at z.
//THIS IS ARGV[2]
	    z[3] = t;
	    for (z[2]=0; z[2]<znodes; z[2]++)
	      for (z[1]=0; z[1]<ynodes; z[1]++)
		for (z[0]=0; z[0]<xnodes; z[0]++) {
		  z_idx4d = lat.GsiteOffset(z)/4;
		  z_idx3d = z_idx4d - vol3d*z[3];
		  z_arr_idx = vol3d*z_idx3d + vol3d*vol3d*(((j+1)/(NRAND))-1);
		  
		  //Loop over sinks at x.
		  x[3] = 0;
		  for (x[2]=0; x[2]<znodes; x[2]++)
		    for (x[1]=0; x[1]<ynodes; x[1]++)
		      for (x[0]=0; x[0]<xnodes; x[0]++) {
			x_idx4d = lat.GsiteOffset(x)/4;
			x_idx3d = x_idx4d - vol3d*x[3];
			
			t2A_full_arr[x_idx3d + z_arr_idx] += qprop_stoch[x_idx4d]*conj(qprop_stoch.rand_src(z_idx4d));
			if((j+1)%(NRAND)==0) {
			  t2A_full_arr[x_idx3d + z_arr_idx] *= 1.0/((j+1));
			}
		      }
		  
		  //Loop over sinks at y.
		  y[3] = t;
		  for (y[2]=0; y[2]<znodes; y[2]++)
		    for (y[1]=0; y[1]<ynodes; y[1]++)
		      for (y[0]=0; y[0]<xnodes; y[0]++) {
			y_idx4d = lat.GsiteOffset(y)/4;
			y_idx3d = y_idx4d - vol3d*y[3];
			
			t3A_full_arr[y_idx3d + z_arr_idx] += qprop_stoch[y_idx4d]*conj(qprop_stoch.rand_src(z_idx4d));
			if((j+1)%(NRAND)==0) {
			  t3A_full_arr[y_idx3d + z_arr_idx] *= 1.0/((j+1));
			}
		      }
		}
	  }
	}
	
	/////////////////////////
	// Compute t2B and t3B //
	/////////////////////////
	
	for(int j=0; j<NRAND; j++) {
//THIS TOO IS ARGV[2]
	  arg_sto.t = t;
	  for (index=0;index<subvol;index++) {
	    QPropWRandWallSrcDilutedSmeared qprop_stoch(lat, &arg_sto, &r_arg, &c_arg, &g_arg, L, index);
	    qprop_stoch.GaussSmearSinkProp(g_arg);	  
	    
	    
	    //Loop over sources at z.
//THIS TOO IS ARGV[2]
	    z[3] = t;
	    for (z[2]=0; z[2]<znodes; z[2]++)
	      for (z[1]=0; z[1]<ynodes; z[1]++)
		for (z[0]=0; z[0]<xnodes; z[0]++) {
		  z_idx4d = lat.GsiteOffset(z)/4;
		  z_idx3d = z_idx4d - vol3d*z[3];
		  z_arr_idx = vol3d*z_idx3d + vol3d*vol3d*(((j+1)/(NRAND))-1); 
		  //Loop over sinks at x.
		  x[3] = 0;
		  for (x[2]=0; x[2]<znodes; x[2]++)
		    for (x[1]=0; x[1]<ynodes; x[1]++)
		      for (x[0]=0; x[0]<xnodes; x[0]++) {
			x_idx4d = lat.GsiteOffset(x)/4;
			x_idx3d = x_idx4d - vol3d*x[3];
			
			t2B_full_arr[x_idx3d + z_arr_idx] += qprop_stoch[x_idx4d]*conj(qprop_stoch.rand_src(z_idx4d));
			if((j+1)%(NRAND)==0) {
			  t2B_full_arr[x_idx3d + z_arr_idx] *= 1.0/((j+1));
			}
		      }
		  
		  //Loop over sinks at y.
//THIS TOO IS ARGV[2]
		  y[3] = t;
		  for (y[2]=0; y[2]<znodes; y[2]++)
		    for (y[1]=0; y[1]<ynodes; y[1]++)
		      for (y[0]=0; y[0]<xnodes; y[0]++) {
			y_idx4d = lat.GsiteOffset(y)/4;
			y_idx3d = y_idx4d - vol3d*y[3];
			
			t3B_full_arr[y_idx3d + z_arr_idx] += qprop_stoch[y_idx4d]*conj(qprop_stoch.rand_src(z_idx4d));
			if((j+1)%(NRAND)==0) {
			  t3B_full_arr[y_idx3d + z_arr_idx] *= 1.0/((j+1));
			}
		      }
		}
	  }
	}
	time[2] = stopwatchReadSeconds();
	///////////////////////////////////////////////
	// Begin summation of trace at time slice t. //
	///////////////////////////////////////////////
	
	// The t1, t1c, t4, and t4c propagators are calculated 'on the fly'
	// within the trace summation.
	
	//Reinitialise Wilson matrix and trace variables	
	for (int r=0; r<allWMs.size(); r++)
	  *(allWMs[r]) *= 0.0;
	for (int r=0; r<allReTrs.size(); r++)
	  *(allReTrs[r]) = 0.0;

	//Sum over x
	x[3] = 0;
	for (x[2]=0; x[2]<znodes; x[2]++)
	  for (x[1]=0; x[1]<ynodes; x[1]++)
	    for (x[0]=0; x[0]<xnodes; x[0]++) {
	      x_idx4d = lat.GsiteOffset(x)/4;
	      x_idx3d = x_idx4d - vol3d*x[3];
	      
	      t1 = qprop_point[x_idx4d];
	      t1c = t1.conj_cp();
	      
	      //Perform t1t1c trace sum for D1 graph.
	      tr = Trace(t1, t1c);
	      t1t1c_re_tr += tr.real();
	      
	      //Sum over y
	      y[3] = t;
	      for (y[2]=0; y[2]<znodes; y[2]++)
		for (y[1]=0; y[1]<ynodes; y[1]++)
		  for (y[0]=0; y[0]<xnodes; y[0]++) {
		    y_idx4d = lat.GsiteOffset(y)/4;
		    y_idx3d = y_idx4d - vol3d*y[3];
		    
		    t4 = qprop_point[y_idx4d];
		    t4c = t4.conj_cp();
		    
		    // Use this condition so that t4t4c is calculated only once
		    // over y per time slice.
		    if (x_idx3d == 0) {
		      //Perform t4t4c trace sum.
		      tr = Trace(t4, t4c);
		      t4t4c_re_tr += tr.real();
		    }
		    
		    //Declare new Wilson Matrix t4*t1c for D2 and compute
		    t4t1c = t4;
		    t4t1c *= t1c;
		    
		    //Declare new Wilson Matrix t1c*t4 for D3 and compute
		    t1t4c = t1;
		    t1t4c *= t4c;
		    
		    //Sum over z
		    z[3] = t;
		    for (z[2]=0; z[2]<znodes; z[2]++)
		      for (z[1]=0; z[1]<ynodes; z[1]++)
			for (z[0]=0; z[0]<xnodes; z[0]++) {
			  z_idx4d = lat.GsiteOffset(z)/4;
			  z_idx3d = z_idx4d - vol3d*z[3];
			  
			  /*
			    ///////////////////////
			    //ONE SOURCE SECTION //
			    ///////////////////////
			    
			    //Perform t4t1c * t2t3c trace sum.		  
			    t2t3c = t2A_full_arr[x_idx3d + vol3d*z_idx3d + a*vol3d*vol3d];
			    t3c   = t3A_full_arr[y_idx3d + vol3d*z_idx3d + a*vol3d*vol3d].conj_cp();
			    t2t3c *= t3c;
			    
			    tr = Trace(t4t1c, t2t3c);
			    d2_1_re_tr[a] += tr.real();
			    
			    //Perform t1ct4c * t3t2c trace sum.
			    t3t2c = t3A_full_arr[y_idx3d + vol3d*z_idx3d + a*vol3d*vol3d];
			    t2c   = t2A_full_arr[x_idx3d + vol3d*z_idx3d + a*vol3d*vol3d].conj_cp();;
			    t3t2c *= t2c;
			    
			    tr = Trace(t1t4c, t3t2c);
			    d3_1_re_tr[a] += tr.real();
			    
			    // Use this condition so that t2t2c is calculated only over
			    // x and z loops per time slice. 
			    if (y_idx3d == 0) {
			      t2  = t2A_full_arr[x_idx3d + vol3d*z_idx3d + a*vol3d*vol3d];
			      t2c = t2A_full_arr[x_idx3d + vol3d*z_idx3d + a*vol3d*vol3d].conj_cp();
			      
			      tr = Trace(t2, t2c);
			      t2t2c_1_re_tr[a] += tr.real();
			    }
			    
			    // Use this condition so that t3t3c is calculated only over
			    // y and z loops per time slice. 
			    if (x_idx3d == 0) {
			      t3  = t3A_full_arr[y_idx3d + vol3d*z_idx3d + a*vol3d*vol3d];
			      t3c = t3A_full_arr[y_idx3d + vol3d*z_idx3d + a*vol3d*vol3d].conj_cp();
			      
			      tr = Trace(t3, t3c);
			      t3t3c_1_re_tr[a] += tr.real();
			    }
			    */
			    ///////////////////////
			    //TWO SOURCE SECTION //
			    ///////////////////////
			    
			    //Perform t4t1c * t2t3c trace sum for D2
			    t2t3c = t2A_full_arr[x_idx3d + vol3d*z_idx3d + vol3d*vol3d];
			    t3c   = t3B_full_arr[y_idx3d + vol3d*z_idx3d + vol3d*vol3d].conj_cp();
			    t2t3c *= t3c;
			    
			    tr = Trace(t4t1c, t2t3c);
			    d2_2_re_tr += tr.real();
			    
			    //Perform t1t4c * t3t2c trace sum for D3
			    //t3t2c = t3A_full_arr[y_idx3d + vol3d*z_idx3d + vol3d*vol3d];
			    //t2c   = t2B_full_arr[x_idx3d + vol3d*z_idx3d + vol3d*vol3d].conj_cp();;
			    //t3t2c *= t2c;
			    
			    //tr = Trace(t1t4c, t3t2c);
			    //d3_2_re_tr += tr.real();
			    
			    // Use this condition so that t2t2c is calculated only over
			    // x and z loops per time slice. 
			    if (y_idx3d == 0) {
			      t2  = t2A_full_arr[x_idx3d + vol3d*z_idx3d + vol3d*vol3d];
			      t2c = t2B_full_arr[x_idx3d + vol3d*z_idx3d + vol3d*vol3d].conj_cp();
			      
			      tr = Trace(t2, t2c);
			      t2t2c_2_re_tr += tr.real();
			    }
			    
			    // Use this condition so that t3t3c is calculated only over
			    // y and z loops per time slice. 
			    if (x_idx3d == 0) {
			      t3  = t3A_full_arr[y_idx3d + vol3d*z_idx3d + vol3d*vol3d];
			      t3c = t3B_full_arr[y_idx3d + vol3d*z_idx3d + vol3d*vol3d].conj_cp();
			      
			      //Perform t3*t3c trace sum for D1 graph.
			      tr = Trace(t3, t3c);
			      t3t3c_2_re_tr += tr.real();
			    }
			}
		  }
	    }

	time[3] = stopwatchReadSeconds();
	
	///////////////////////////
	// Write traces to file. //
	///////////////////////////
	if (allReTrs.size() != allReTrNames.size())
	  cerr << "writing to file: data & name mismatch!" << endl;
	for (int r=0; r<allReTrs.size(); r++) {
	  sprintf(file, DATAPATH"%s_%s_sto_%d-%d_%d_%d_%d_%d%d%d.dat", is_qu.c_str(), allReTrNames[r].c_str(), 
	          NSITES_3D, NSITES_T, sweep_counter, t, NRAND, L[0], L[1], L[2]);
	  FILE* f = Fopen(file, "a");
	  Fprintf(f,"%d %d %d %.16e\n", sweep_counter, t, 0, allReTrs[r]); 
	  Fclose(f);
	}
	
	//}
      /////////////////////////
      // End loop over time. //
      /////////////////////////
    }
  }
  ////////////////////
  // End simulation //
  ////////////////////
  
  //free heap variables
  for (int r=0; r<allWMarrs.size(); r++)
    delete[] allWMarrs[r];

  return 0;
}
