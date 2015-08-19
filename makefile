TOPLVL = .
include makefile.defs

ifeq ($(BUILD_GPU),yes)
ifeq (`which nvcc`,)
$(error nvcc is not in your path. Add it before running make.)
endif
endif

all: 
	$(MAKE) -C include
	$(MAKE) -C src

cps: 
	cd $(CPS_SRC)
ifeq ($(BUILD_SCALAR),yes)
	cd $(CPSCPU); $(CPS_SRC)/configure CXXFLAGS='-Wno-write-strings'; make
endif
ifeq ($(BUILD_GPU),yes)
	cd $(CPSGPU); $(CPS_SRC)/configure \
CXXFLAGS='-Wno-write-strings -DUSEQUDA -I$(QUDA)/include -I$(CUDA)/include'; make
endif

clean:
	$(MAKE) clean -C include
	$(MAKE) clean -C src

cleanlib:
	$(MAKE) clean -c QUDA-CPS_v5_0_26
