CXX=g++
MPICXX=mpic++
PARAM=-fopenmp 
ARMAPARAMS=-DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -DARMA_USE_HDF5 -larmadillo -lopenblas -lhdf5 -lsuperlu
OPTFLAGS= -O3 -march=native -Ofast#-frename-registers  -funroll-loops -ffast-math

#PARAM=-larmadillo  -lfftw3_omp -lfftw3 -lm -O3 -fopenmp 
ALLFILES=main.cpp IncBeam.cpp SetConfig.cpp PlasProfile.cpp
NUMBERS=Numbers.cpp SetConfig.cpp
CXXFLAGS= -I/usr/include/hdf5/serial/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial -std=c++14 
 
Out=../main
OutTri=../numbers

all:
	$(CXX) -o $(Out) $(ALLFILES) $(PARAM) $(ARMAPARAMS) $(OPTFLAGS) $(CXXFLAGS)
	../main

mpi: 
	$(MPICXX) -o $(Out) $(ALLFILES) $(PARAM) $(ARMAPARAMS) $(OPTFLAGS) $(CXXFLAGS)
	../main
	
debug: 
	$(MPICXX) -g -o $(Out) $(ALLFILES) $(PARAM) $(ARMAPARAMS) $(OPTFLAGS) $(CXXFLAGS)
	../main	

generate:
	$(CXX) -o $(Out) $(ALLFILES) $(PARAM) $(ARMAPARAMS) $(OPTFLAGS) $(CXXFLAGS) -fprofile-generate
	../main
	
use:
	$(CXX) -o $(Out) $(ALLFILES) $(PARAM) $(ARMAPARAMS) $(OPTFLAGS) $(CXXFLAGS) -fprofile-use -fprofile-correction
	../main

numbers:
	$(CXX) -o $(OutTri) $(NUMBERS) $(PARAM) $(ARMAPARAMS) $(OPTFLAGS) $(CXXFLAGS)
	../numbers