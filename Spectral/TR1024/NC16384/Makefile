
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
FLAGS = -V14.1 -fast -r8 
# FLAGS_CUDA sets the following compilation options for CUDA: Toolkit version, Compute version,
# target processor type 
FLAGS_CUDA =#-Mcuda=cuda5.0,cc35 -tp:x64 -ta=nvidia,time
#Sets compiler for extension type
F90=pgf90

# Variables 
SOURCES = mk_32.f90
OBJECTS = $(SOURCES: .f90=.o)

all: build

build: Spectral

Spectral: $(OBJECTS)

	$(F90) $(FLAGS) $(FLAGS_CUDA) $(OBJECTS) -o $@ 
.f90:
	$(F90) $(FLAGS) $(FLAGS_CUDA) $< -o $@ 
run: Spectral
	./Spectral
tidy:
	@find * | egrep "#" | xargs rm -f
	@find * | egrep "\~" | xargs rm -f	

# Cleans 
.PHONY: clean
clean:
	rm *.mod Spectral_OMP
