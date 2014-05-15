
########################################################################
# CUDA Makefile                                                        #
# Author       : Daniel Savage                                         #
# Date         : 6/30/2013                                             #
# Description  : Compiles Spectral code with CUDA implementation       #
# Information  : Makefile http://mrbook.org/tutorials/make/            #
#              : Compiler options: http://www.pgroup.com/doc/pgiug.pdf # 
########################################################################

# FLAGS sets the following compilation options: compiler version, -fast: combines 
# several compilation optimization features, Mconcur=innermost: auto paralization of loops 
FLAGS = -V13.4 -fast -Mconcur=innermost 
# FLAGS_CUDA sets the following compilation options for CUDA: Toolkit version, Compute version,
# target processor type 
FLAGS_CUDA =-Mcuda=cuda5.0,cc35 -tp:x64 -ta=nvidia,time
#Sets compiler for extension type
F90=pgfortran

# Variables 
SOURCES = Variables.f90 CUDA_Kernels.f90 CUDA_Spectral.f90
OBJECTS = $(SOURCES: .f90=.o)

all: build

build: Spectral_GPU

Spectral_GPU: $(OBJECTS)
	$(F90) $(FLAGS) $(FLAGS_CUDA) $(OBJECTS) -o $@ 
.f90:
	$(F90) $(FLAGS) $(FLAGS_CUDA) $< -o $@ 
run: Spectral_GPU
	./Spectral_GPU
	make tidy
tidy:
	@find * | egrep "#" | xargs rm -f
	@find * | egrep "\~" | xargs rm -f	

# Cleans 
.PHONY: clean
clean:
	rm *.mod *.o Spectral_GPU
