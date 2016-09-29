
BOUT_TOP = ../../../..

SOURCEC         = cartesianmesh.cxx
SOURCEH         = cartesianmesh.hxx
TARGET          = lib

include $(BOUT_TOP)/make.config

#lib: makefile $(BOUT_TOP)/make.config $(BOUT_TOP)/include $(BOUT_TOP)/lib $(OBJ) generated_stencils.cxx

cartesianmesh.o:generated_stencils.cxx

generated_stencils.cxx:stencils_cleaned.cxx gen_stencils.py
	python gen_stencils.py > $@

stencils_cleaned.cxx:../../stencils.cxx
	gcc -E $< > $@
