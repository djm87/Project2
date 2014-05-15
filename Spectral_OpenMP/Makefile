
########################################################################
# CUDA Makefile                                                        #
# Author       : Daniel Savage                                         #
# Date         : 6/30/2013                                             #
# Description  : Compiles and builds results for test cases OpenMP     #
# Information  : Makefile http://mrbook.org/tutorials/make/            #
#              : Compiler options: http://www.pgroup.com/doc/pgiug.pdf # 
########################################################################
<<<<<<< HEAD

export OMP_NUM_THREADS=10	

# Project folders that contain test cases
PROJECTS = $(shell find TR1024/NC256 TR1024/NC512 TR1024/NC1024 TR1024/NC2048 TR1024/NC4096 TR1024/NC8192 TR1024/NC16384 -name Makefile)
	

%.ph_build :
	+@$(MAKE) -C $(dir $*) $(MAKECMDGOALS)

%.ph_clean : 
	+@$(MAKE) -C $(dir $*) clean $(USE_DEVICE)

all:  $(addsuffix .ph_build,$(PROJECTS))
	@echo "Finished building Spectral_OMP Cases"

build: $(addsuffix .ph_build,$(PROJECTS))
       

run: $(addsuffix .ph_build,$(PROJECTS))
	make tidy
	make plot 
	@echo "Finished building Spectral_OMP Cases"

plot: 
	./gnuplot.bash
	eog Spectral_OMP_Scaling.png &

=======
export OMP_NUM_THREADS=10
# FLAGS sets the following compilation options: compiler version, -fast: combines 
# several compilation optimization features, Mconcur=innermost: auto paralization of loops 
FLAGS = -V14.1 -fast -r8 -mp 
# FLAGS_CUDA sets the following compilation options for CUDA: Toolkit version, Compute version,
# target processor type 
FLAGS_CUDA =#-Mcuda=cuda5.0,cc35 -tp:x64 -ta=nvidia,time
#Sets compiler for extension type
F90=pgf90

# Variables 
SOURCES = Spectral_OMP.f90
OBJECTS = $(SOURCES: .f90=.o)

all: build

build: Spectral_OMP

Spectral_OMP: $(OBJECTS)
	$(F90) $(FLAGS) $(FLAGS_CUDA) $(OBJECTS) -o $@ 
.f90:
	$(F90) $(FLAGS) $(FLAGS_CUDA) $< -o $@ 
run: Spectral_OMP
	./Spectral_OMP
>>>>>>> Develop_GPU_CPSOLV
tidy:
	@find * | egrep "#" | xargs rm -f
	@find * | egrep "\~" | xargs rm -f

clean: tidy $(addsuffix .ph_clean,$(PROJECTS)) 
