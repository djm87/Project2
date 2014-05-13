
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
FLAGS = -V13.4  -fast -mp #-Minfo=ccff  #-Mvect=levels:4 -lpthread#-Mipa=fast,inline -fast
# FLAGS_CUDA sets the following compilation options for CUDA: Toolkit version, Compute version,
# target processor type 
FLAGS_CUDA = #-Mcuda=cuda5.0,cc35 -tp:x64

#Sets compiler for extension type
F90=pgfortran

#Variables 
SOURCES = Variables.f90 CPSOLV.f90 cpty.f90  #CUDA_Kernels.f90
OBJECTS = $(SOURCES: .f90=.o)
#EXECUTABLE = CUDA_Parent

all: build

build: OMP_Parent


OMP_Parent: $(OBJECTS)
	$(F90) $(FLAGS) $(OBJECTS) -o $@ #$(FLAGS_CUDA) 
.f90:
	$(F90) $(FLAGS)  $< -o $@ 
run: OMP_Parent
	./OMP_Parent
	

# Cleans 
.PHONY: clean
clean:
	rm *.mod *.o OMP_Parent
