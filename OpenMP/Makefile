
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
PROJECTS ?= $(shell find NC256 NC512 NC1024 NC2048 NC4096 NC8192 NC16384 -name Makefile)

%.ph_build :
	+@$(MAKE) -C $(dir $*) $(MAKECMDGOALS)

%.ph_clean : 
	+@$(MAKE) -C $(dir $*) clean $(USE_DEVICE)

all:  $(addsuffix .ph_build,$(PROJECTS))
	@echo "Finished building OpenMP Cases"

build: $(addsuffix .ph_build,$(PROJECTS))
       

run: $(addsuffix .ph_build,$(PROJECTS))
	./gnuplot.bash
	make tidy

	@echo "Finished Running OpenMP Cases"
tidy:
	@find * | egrep "#" | xargs rm -f
	@find * | egrep "\~" | xargs rm -f

clean: tidy $(addsuffix .ph_clean,$(PROJECTS))


