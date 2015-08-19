//args.h
//functions to initialize CPS args
#include <alg/alg_hmd.h>   //Hybrid molecular dynamics
#include "params.h"

void setup_do_arg(DoArg& do_arg, int seed);

void setup_hmd_arg(HmdArg& hmd_arg);

void setup_qpropwarg_cg(QPropWArg& q_arg);

void setup_g_arg(QPropWGaussArg& g_arg);

