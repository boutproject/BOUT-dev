
BOUT_TOP = ../../../..

SOURCEC         = aiolosmesh.cxx
SOURCEH         = aiolosmesh.hxx
TARGET          = lib

include $(BOUT_TOP)/make.config

generated=generated_header.hxx generated_derivs.cxx generated_stencils.cxx generated_init.cxx
#lib: makefile $(BOUT_TOP)/make.config $(BOUT_TOP)/include $(BOUT_TOP)/lib $(OBJ) generated_stencils.cxx

aiolosmesh.o:$(generated)

tables_cleaned.cxx:../../index_derivs.cxx
	@ grep ' * Lookup tables of functions. Map between names, codes and functions' -B 1 -A 10000 $< |\
	grep ' * Routines to use the above tables to map between function codes, names' -A 2 -B 10000 |\
	gcc -E - > $@

:generated_derivs.cxx

stencils_cleaned.cxx:../../index_derivs.cxx
	@ grep ' * All expect to have an input g' -B 2 -A 10000 $< | \
	grep ' * Lookup tables of functions. Map between names, codes and functions' -A 1 -B 10000 |\
	cat  |\
	gcc -E - > $@

$(generated):gen_derivs.py \
		tables_cleaned.cxx gen_stencils.py stencils_cleaned.cxx \
		common.py
	python $< > generated_derivs.cxx || ( rm $@ ; exit 1)
	@(clang-format -i $(generated) || astyle -n $(generated) ) || exit 0 # ignore failure

clean:: clean-generated

clean-generated:
	@rm -f $(generated)
