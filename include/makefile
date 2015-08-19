include ../makefile.defs

ifeq ($(BUILD_SCALAR),yes)
ALL += stopwatch.o FFTW_functions.o
WHICHFLAGS = $(CXXFLAGS_CPU)
endif
ifeq ($(BUILD_GPU),yes)
ALL += tr_ABdag.o
endif
ALL ?= null
ifneq ($(ALL),null)
WHICHFLAGS ?= $(CXXFLAGS_GPU)
ALL += #CPS_utils.o#CPS_utils_1pion.o CPS_utils_sum.o
endif

.PHONY: all clean $(ALL)

all: $(ALL)

tr_ABdag.o:
	nvcc -c -I$(CUDA)/include tr_ABdag.cu

stopwatch.o: stopwatch.h stopwatch.cpp ; $(CXX) -c -I. stopwatch.cpp

FFTW_functions.o: FFTW_functions.h FFTW_functions.cpp mom1D.h mom1D.cpp ; $(CXX) -c $(CXXFLAGS_FFTW) FFTW_functions.cpp mom1D.cpp
#FFTW_functions.o: test_FFT_wtf.cpp mom1D.cpp ; $(CXX) -c $(CXXFLAGS_FFTW) test_FFT_wtf.cpp mom1D.cpp

CPS_utils_1pion.o: CPS_utils_1pion.h CPS_utils.cpp ; $(CXX) -c $(CXXFLAGS) $(WHICHFLAGS) CPS_utils_1pion.h CPS_utils.cpp
CPS_utils_sum.o: CPS_utils_sum.h CPS_utils.cpp ; $(CXX) -c $(CXXFLAGS) $(WHICHFLAGS) CPS_utils_sum.h CPS_utils.cpp
CPS_utils.o: CPS_utils.h CPS_utils.cpp ; $(CXX) -c $(CXXFLAGS) $(WHICHFLAGS) CPS_utils.cpp

null:
	@echo include: Nothing built.

clean:
	rm -f *.o
