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
#include <config.h>
#include <util/vector.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <comms/scu.h>
#include <util/qcdio.h>
#include <alg/qpropw.h>
//#include <alg/fix_gauge_arg.h>
//#include <alg/alg_fix_gauge.h>
#include <util/random.h>
#include <util/WriteLatticePar.h>
#include <util/ReadLatticePar.h>
#include <alg/alg_ghb.h>   //Heat bath == quenched
#include <alg/ghb_arg.h>   // "  "  "  "  
#include <alg/alg_hmd.h>   //Hybrid molecular dynamics == unquenched

//FFTW
#include <fftw3.h>

//local include directory
#include <FFTW_functions.h>
#include <stopwatch.h>
#include <CPS_utils.h>
#include <mom1D.h>

//////////////
//parameters//
//////////////

//simulation parameters
#define STOP_RSD 1.0e-6 //fractional error for the inverter to stop at
#define MAX_NUM_ITER 50000 //fail if an inversion takes this many steps
#define NDATA 36 //number of measurements to collect
#define NSKIP 200 //number of lattice steps between measurements. Recommended: quenched ~200, unquenched ~20
#define NTHERM 500 //number of thermalizations. Recommended: quenched ~500, unquenched ~5000

//physical parameters
#define NSITES_3D 10        //lattice size, space dimensions
#define NSITES_T 20         //lattice size, time dimension
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

//stochastic source parameters
#define SRC_TYPE ZTWOXZTWO  //GAUSS, UONE, ZTWO, ZTWOXTWO
#define NRAND 10 //number of sources per timestep

//path to look for lattice objects
#define LATT_PATH "../../latt_configs/"
//path to send output files to
#define DATAPATH "../../data"

/*
If you have run a simulation before with these parameters unchanged:
        Quenched, BETA, NSITES_3D, NSITES_T
                or
        Unquenched, BETA, MASS, NSITES_3D, NSITES_T
uncomment READ. The lattice objects will be read from LATT_PATH
instead of re-generated, which is significantly faster. If you are unsure,
grep in LATT_PATH for a set of .dat files with the relevant parameters.
*/
#define READ

//for file i/o and debugging purposes
//#define NAME(a) #a
