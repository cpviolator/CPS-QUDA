TOPLVL = ../..
include ../../makefile.defs

.PHONY: clean

#edit this list to change targets
CPU_CPP = \
1hadr.CPU.mom.cpp \
1hadr.CPU.ps.cpp \

#edit this list to change targets
GPU_CPP = \
1hadr.GPU.mom.cpp \
1hadr.GPU.ps.cpp \

ifeq ($(BUILD_SCALAR),yes)
ALL += cpu_exe
endif
ifeq ($(BUILD_GPU),yes)
ALL += gpu_exe
endif
ALL ?= null

all: $(ALL)

cpu_exe: cpu_o $(CPSCPU)/cps.a ; $(foreach nm, $(CPU_CPP), $(CXX) -o $(nm:.cpp=.exe) $(nm:.cpp=.o) \
args.o $(CPSCPU)/cps.a;)

cpu_o: $(CPU_CPP) ; $(foreach nm, $(CPU_CPP) args.cpp, $(CXX) -c $(CXXFLAGS) $(CXXFLAGS_CPU) $(nm);)

gpu_exe: gpu_o $(CPSGPU)/cps.a ; $(foreach nm, $(GPU_CPP), $(CXX) -o $(nm:.cpp=.exe) $(nm:.cpp=.o) \
args.o $(CPSGPU)/cps.a $(QUDA)/lib/libquda.a -L$(CUDA)/lib64 -lcudart;)

gpu_o: $(GPU_CPP) ; $(foreach nm, $(GPU_CPP) args.cpp, nvcc -c $(CXXFLAGS) $(CXXFLAGS_GPU) $(nm) args.cpp;)

null:
	@echo one_hadron: Nothing built.

clean: 
	rm -f *.o *.exe
