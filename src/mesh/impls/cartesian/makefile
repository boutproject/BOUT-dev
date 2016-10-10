
BOUT_TOP = ../../../..

SOURCEC         = cartesianmesh.cxx
SOURCEH         = cartesianmesh.hxx
TARGET          = lib

include $(BOUT_TOP)/make.config

#lib: makefile $(BOUT_TOP)/make.config $(BOUT_TOP)/include $(BOUT_TOP)/lib $(OBJ) generated_stencils.cxx

cartesianmesh.o:generated_stencils.cxx generated_derivs.cxx

tables_cleaned.cxx:../../diff_tables.cxx
	gcc -E $< > $@

#generated_stencils.cxx:gen_stencils.py stencils_cleaned.cxx
#	python $< > $@

stencils_cleaned.cxx:../../stencils.cxx
	gcc -E $< > $@

generated_derivs.cxx:gen_derivs.py tables_cleaned.cxx gen_stencils.py stencils_cleaned.cxx
	python $< > $@
