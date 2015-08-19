Welcome to the QUDA enabled Columbia Physics System

	GPU SUPPORT

This CPS library has been modified to offer limited GPU support. Currently 
supported fermion inverters are listed below:

	  1. Clover Fermions
	  2. Wilson Fermions*

(*One may set the clover coefficient `do_arg.clover_coeff = 0.0` to simulate
Wilson fermions. A less inelegant solution is forthcoming.)

One may utilise the CG or BICGSTAB QUDA inverter type simply by setting the
CPS inverter type as required. An new CPS enumerator QUDA_GCR_INVERTER 
has been included to uitise the QUDA GCR inverter.

The functions added to CPS that enable QUDA support are contained in 
`src/util/quda_functions/quda_functions.C` and the header file is in 
include/util/quda_functions.h. The code changes to CPS are delineated by 
the compiler flag `#ifdef USEQUDA`

Examples on how to use the GPU enabled library are given in this release
for both single hadron and two hadron channels. The modifications 
required to allow existing code to use this library should be minimal.
For details, lease refer to the differences between GPU enabled and CPU 
only files, highlighted in our production code package. Instructions for 
Configuring and making this package are given below.

	CONFIGURING AND COMPILATION

To configure this library on scalar, BGQ, or other CPU architecture, 
please refer to the README file in CPS. To configure for GPU, run:
```sh
$ ./configure CXXFLAGS=" -DUSEQUDA -I/path/to/quda/include -I/path/to/cuda/include" 
```
from the `QUDA-CPS_v5_0_26/` directory. You must already have working 
installations of QUDA and CUDA.

CPS can be configured for multiple architectures from the same source tree.
In this case, we recommend keeping the source tree and build trees separate.
An example build for scalar and GPU follows:
```sh
$ mkdir cpscpu cpsgpu 
$ cd cpscpu 
$ ../QUDA-CPS_v5_0_26/configure CXXFLAGS="-Wno-write-strings" 
$ make 
$ cd ../cpsgpu 
$ ../QUDA-CPS_v5_0_26/configure CXXFLAGS="-Wno-write-strings 
-DUSEQUDA -I/path/to/quda/include -I/path/to/cuda/include" 
$ make 
```
(The `-Wno-write-strings` flag is optional, but recommended to silence
compiler warnings.)

	BUILD INSTRUCTIONS

The CPS library should be built separately from the simulation code that
composes the rest of this package. To build the rest of the package 
after building CPS, edit the top-level `makefile.defs` and then run `make`.
A top-level `make` will fail if CPS is unbuilt or if `makefile.defs` is not
properly configured. An example `makefile.defs` corresponding to the above
configuration example is:
```make
CPSSRC = /path/to/QUDA-CPS/QUDA-CPS_v5_0_26#CPS source code
CPSGPU = /path/to/QUDA-CPS/cpsgpu#built GPU version
CPSCPU = /path/to/QUDA-CPS/cpscpu#built scalar version
QUDA = /path/to/quda
CUDA = /path/to/cuda
FFTW = /path/to/fftw

BUILD_GPU = yes
BUILD_SCALAR = yes
```

Dependencies for all (GPU and scalar) simulation code are: QUDA-CPS 5.0.26 ; FFTW 3

Dependencies for GPU simulation code are: QUDA 0.5.0+ ; CUDA 4.1+


	EXTRA FEATURES

Several extra features not found in the current CPS release are included 
in this package:

   1. Stout Gaussian Kernel Link Smearing.
      This routine applies the link smearing procedure outlined in
      http://arxiv.org/pdf/hep-lat/0311018v1.pdf
      For the time being only GKLS is supported. Incorporating stout
      smearing into the HMC routines is ongoing.

   2. Z_2 X Z_2 stochastic sources.
      We have included an option to use Z_2 X Z_2 stochastic sources.
      The enum type is ZTWOXZTWO.   

   3. Diluted/Smeared Stochastic Sources
      The QPropWRandWall class has been extended to include interlace 
      dilution and/or Jacobi smearing. Please refer to
      src/alg/alg_qprop/QPropW.C for the necessary arguments one needs
      to pass to these classes.

   4. Diluted/Smeared Momentum Shell sources
      The QPropWMom class now inherits directly from QPropW. It includes
      a new member function 'mom_src[]' that is analogous to the 
      the QPropWRand member function 'rand_src[]'. This new member function 
      is necessary for the calculation of momentum shell sources. The
      structure of the new QPropWMom and derived classes is the same as
      the new The QPropWRandWall classes in that one may employ dilution, 
      Jacobi smearing, or both.

      
The current release of CPS contains a bug in the function GFWallSource(). We 
have fixed the bug in this release. All changes to the current vanilla 
version of CPS are delineated by

	Begin QUDA_CPS
	      ... code changes ...
	End QUDA_CPS

to allow for easy identification. We request that any changes or additions 
in pushes to the public repository be delineated in the same fashion to ease 
the merging procedure.

    	CONTACT

We welcome all questions, comments, suggestions, complaints... but most 
of all, we like bug reports. Thank you.

   [Dean Howarth](https://github.com/cpviolator) howard3 at rpi dot edu
   [Matthew Bernstein](https://github.com/bernsm3) bernsm3 at rpi dot edu
