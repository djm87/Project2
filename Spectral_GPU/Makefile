
########################################################################
# CUDA Makefile                                                        #
# Author       : Daniel Savage                                         #
# Date         : 6/30/2013                                             #
# Description  : Compiles and builds results for test cases OpenMP     #
# Information  : Makefile http://mrbook.org/tutorials/make/            #
#              : Compiler options: http://www.pgroup.com/doc/pgiug.pdf # 
########################################################################

export OMP_NUM_THREADS=10	

# Project folders that contain test cases
PROJECTS ?= $(shell find PROJECTS ?= $(shell find TR1024/NC256 TR1024/NC512 TR1024/NC1024 TR1024/NC2048 TR1024/NC4096 TR1024/NC8192 TR1024/NC16384 -name Makefile) -name Makefile)
	

%.ph_build :
	+@$(MAKE) -C $(dir $*) $(MAKECMDGOALS)

%.ph_clean : 
	+@$(MAKE) -C $(dir $*) clean $(USE_DEVICE)

all:  $(addsuffix .ph_build,$(PROJECTS))
	@echo "Finished building Spectral_GPU Cases"

build: $(addsuffix .ph_build,$(PROJECTS))
       

run: $(addsuffix .ph_build,$(PROJECTS))
	make tidy
	make plot 
	@echo "Finished building Spectral_GPU Cases"

plot: 
	./gnuplot.bash
	eog Spectral_GPU_Scaling.png &

tidy:
	@find * | egrep "#" | xargs rm -f
	@find * | egrep "\~" | xargs rm -f

clean: tidy $(addsuffix .ph_clean,$(PROJECTS)) 
