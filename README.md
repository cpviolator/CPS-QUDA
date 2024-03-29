Welcome to the QUDA enabled Columbia Physics System!

	GPU SUPPORT Mx=y solve

This CPS library has been modified to offer limited GPU support for 
fermion inverters. Currently supported fermions are:

	  1. Clover Fermions
	  2. Wilson Fermions
	  3. Twisted Mass (Wilson)

One may utilise the `CG` or `BICGSTAB` QUDA inverter type simply by setting the
CPS inverter type as required. A new CPS enumerator `QUDA_GCR_INVERTER` 
has been included to utilise the QUDA `GCR` inverter.

	GPU SUPPORT M^{\dagger}Mx=y solve

We have implemented partial GPU support in the HMC routines. The fermion
matrix inversion is done on the device, but the gauge force and 
pseudo-heatbath routines are presently still done on the host.
Currently supported fermions are:

	  1. Wilson Fermions
	  2. Twisted Mass (Wilson)

N.B. When using the twisted mass fermion action, the epsilon parameter
must be set as a Global Job Parameter in the DoArg structure. Simply
set

do_arg.epsilonTm = Your Value;

at compile time.

All extra files and functions needed to enable QUDA support are contained in the
directory `src/util/quda_interface` and the header file 
`include/util/quda_interface.h`. Please review the file 
`src/util/quda_interface/quda_interface.cpp` to adjust the QUDA settings as you need.
All modified files in CPS pertaining to GPU support are delineated by the compiler flag 

`#ifdef USEQUDA`

and all changes, including the extra features, are delineated by

	Begin QUDA_CPS
	      ... code changes ...
	End QUDA_CPS

Examples on how to use the GPU enabled library are given in 
this release for both single hadron and two hadron channels. The modifications 
required to allow existing code to use this library should be minimal.
For details, please refer to the differences between GPU enabled and CPU 
only files, highlighted in our production code package. Instructions for 
configuring and making this package are given below.

	CONFIGURING AND COMPILATION

To configure this library on scalar, BG/Q, or other CPU architecture, 
please refer to the `README` file in CPS. To configure for GPU, run:
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

Dependencies for GPU simulation code are: QUDA 0.5.0x 0.6.x ; CUDA 4.1+


	EXTRA FEATURES

Several extra features not found in the current CPS release are included 
in this package:

   1. Stout Gaussian Kernel Link Smearing

      This routine applies the link smearing procedure outlined in
      http://arxiv.org/pdf/hep-lat/0311018v1.pdf
      For the time being only `GKLS_STOUT` is supported. Incorporating stout
      smearing into the HMC routines is ongoing.

   2. Z<sub>2</sub> X Z<sub>2</sub> stochastic sources

      We have included an option to use Z<sub>2</sub> X Z<sub>2</sub> stochastic sources.
      The enum type is `ZTWOXZTWO`.   

   3. Diluted/Smeared Stochastic Sources

      The `QPropWRandWall` class has been extended to include interlace 
      dilution and/or Jacobi smearing. Please refer to
      `src/alg/alg_qprop/QPropW.C` for the necessary arguments one needs
      to pass to these classes.

   4. Diluted/Smeared Momentum Shell sources

      The `QPropWMom` class now inherits directly from `QPropW`. It includes
      a new member function `mom_src[src_phase]` that is analogous to the 
      the `QPropWRand` member function `rand_src[src_value]`. This new member function 
      is necessary for the calculation of momentum shell sources. The
      structure of the new `QPropWMom` and derived classes is the same as
      the new `QPropWRandWall` classes in that one may employ dilution, 
      Jacobi smearing, or both.

   5. Smearing and Gauge Fixing

      In order to properly apply the gauge fixing matrices to the source
      and sink fermion fields, we needed to change the order in which 
      gauge fixing and smearing are performed. For the source this is a 
      simple case of reodering some function calls, but for the sink we
      must move the gauge fixing routine out of `QPropW.Run()` and into
      `GaussSinkSmearProp()`. As a result, when using gauge fixing and 
      smeared operators, one must set the new `QPropWArg` flag 
      `smeared_snk = 1` so that the internal routines can perform gauge
      fixing in the correct order.   		       
      
      N.B. The current release of CPS contains a bug in the function 
      `GFWallSource()`. We have fixed the bug in this release. 

        DEVELOPMENT      

All changes to the current vanilla version of CPS are delineated by

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
