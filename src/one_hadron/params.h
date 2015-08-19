//params.h

///////////////////
//common includes//
///////////////////

//C++ STL
#include <cstdlib>
#include <sys/time.h>
#include <cmath>
#include <string>
#include <vector>

//CPS
#include <config.h>        //Autoconf thing
#include <util/vector.h>   //Use utilities for Vector data structure.
#include <util/lattice.h>  //Use utilities for Lattice data structure.
#include <util/gjp.h>      //Global Job Paramters
#include <comms/scu.h>     //communications
#include <util/qcdio.h>
#include <alg/qpropw.h>    //The quark propagator class
#include <util/random.h>   //Use CPS's RNG
#include <util/WriteLatticePar.h>
#include <util/ReadLatticePar.h>
#include <alg/alg_ghb.h>   //Heat bath == quenched
#include <alg/ghb_arg.h>   // "  "  "  "  
#include <alg/alg_hmd.h>   //Hybrid molecular dynamics == unquenched

//local include directory
#include <CPS_utils.h>
#include <mom3D.h>

//////////////
//parameters//
//////////////

//simulation parameters
#define STOP_RSD 1.0e-6 //fractional error for the BiCGSTAB solver to stop at
#define MAX_NUM_ITER 50000 //fail if a solver takes this many steps
#define NDATA 100 //number of measurements to collect
#define NSKIP 200 //number of lattice steps between measurements. Recommended: quenched ~200, unquenched ~20
#define NTHERM 500 //number of thermalizations. Recommended: quenched ~500, unquenched ~5000

//physical parameters
#define NSITES_3D 16        //lattice size, space dimensions
#define NSITES_T 32         //lattice size, time dimension
#define BETA 5.45
#define MASS -0.775
#define QUENCH              //quenched or unquenched simulation

//Gaussian kernel link smearing parameters
#define GAUSS_LS_TYPE GKLS_STOUT     //Link smearing
#define GAUSS_LS_COEFF 0.1           //Link smearing
#define GAUSS_LS_N 3                 //Link smearing hits
#define GAUSS_N 32                   //Source/Sink smearing hits
#define KAPPA 0.125

//Inverter parameters
#define INVERTER_TYPE BICGSTAB //BICGSTAB or CG
#define BICGSTAB_N 1 //1-10, see CPS doc for details

//Hybrid Monte Carlo (HMC) parameters (unquenched only)
#define TAU_INIT 0.5
//#define TAU_SIM 0.5 //unused
#define TAU_STEPS 70

//phase (momentum source only)
//P = { P1, P2, P3 }
#define P1 0
#define P2 0
#define P3 0

//path to look for lattice objects
#define LATT_PATH "../../latt_configs/"
//path to send output files to
#define DATAPATH "../../data"

/*
If you have run a simulation before with these parameters unchanged:
	QUENCH=true, BETA, NSITES_3D, NSITES_T
		or
	QUENCH=false, BETA, MASS, NSITES_3D, NSITES_T
uncomment READ. The lattice objects will be read from LATT_PATH
instead of re-generated, which is significantly faster. If you are unsure,
grep in /latt_configs for a set of .dat files with the relevant parameters.
*/

//#define READ
