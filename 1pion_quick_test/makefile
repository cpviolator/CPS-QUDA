CPSCPU=/path/to/CPU CPS build
CPSGPU=/path/to/GPU CPS build
QUDA=/path/to/quda
CUDA=/path/to/cuda

NVCC = nvcc -arch=sm_20
CXX = g++
CPUCXXFLAGS = -I. -I$(CPSCPU) -I$(CPSCPU)/include
GPUCXXFLAGS = -I. -I$(CPSGPU) -I$(CPSGPU)/include -I$(QUDA)/include

all : 1pion.exe 1pion.GPU.exe

1pion.o: 1pion.cpp ; ${CXX} -c ${CPUCXXFLAGS} 1pion.cpp

1pion.exe: 1pion.o $(CPSCPU)/cps.a ; ${CXX} -o 1pion.exe \
   1pion.o $(CPSCPU)/cps.a\
   -lgmpxx -lgmp -lmpfr

1pion.GPU.o: 1pion.GPU.cpp ; ${CXX} -c ${GPUCXXFLAGS} 1pion.GPU.cpp

1pion.GPU.exe: 1pion.GPU.o $(CPSGPU)/cps.a ; ${CXX} -o 1pion.GPU.exe \
   1pion.GPU.o $(CPSGPU)/cps.a $(QUDA)/lib/libquda.a -L$(CUDA)/lib64 -lcudart\
   -lgmpxx -lgmp -lmpfr

clean: 
	rm -f *.o
	rm -f *.exe

