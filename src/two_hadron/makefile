TOPLVL = ../..
include ../../makefile.defs

.PHONY: all clean

#edit this list to change targets
CPP = \
2hadr.stoch.dil_smear.cpp \
2hadr.ps.cpp \
2hadr.mom.6d.cpp \
2hadr.mom.9d.cpp \

ifeq ($(BUILD_SCALAR),yes)
ALL = exe
endif
ALL ?= null
CXXFLAGS += $(CXXFLAGS_CPU) $(CXXFLAGS_FFTW)

all: $(ALL) ;

exe: obj ; $(foreach nm, $(CPP), $(CXX) -o $(nm:.cpp=.exe) $(nm:.cpp=.o) \
$(TOPLVL)/include/stopwatch.o $(TOPLVL)/include/FFTW_functions.o $(TOPLVL)/include/mom3D.o\
 args.o $(CPSCPU)/cps.a $(FFTW)/lib/libfftw3.a;)

obj: $(CPP) ; $(foreach nm, $(CPP) args.cpp, $(CXX) -c $(CXXFLAGS) $(nm);)

null:
	@echo two_hadron: Nothing built.

clean:
	rm -f *.o *.exe
