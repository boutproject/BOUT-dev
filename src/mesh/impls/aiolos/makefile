
BOUT_TOP = ../../../..

SOURCEC         = aiolosmesh.cxx
SOURCEH         = aiolosmesh.hxx
TARGET          = lib

include $(BOUT_TOP)/make.config

generated=generated_header.hxx generated_derivs.cxx generated_stencils.cxx generated_init.cxx
#lib: makefile $(BOUT_TOP)/make.config $(BOUT_TOP)/include $(BOUT_TOP)/lib $(OBJ) generated_stencils.cxx

aiolosmesh.o:$(generated)

tables_cleaned.cxx:../../diff_tables.cxx
	gcc -E $< > $@

:generated_derivs.cxx

stencils_cleaned.cxx:../../stencils.cxx
	gcc -E $< > $@

$(generated):gen_derivs.py \
		tables_cleaned.cxx gen_stencils.py stencils_cleaned.cxx \
		common.py
	python $< > generated_derivs.cxx || ( rm $@ ; exit 1)
	@(astyle -n $(generated) ) || exit 0 # ignore failure

clean:: clean-generated

clean-generated:
	@rm $(generated)
