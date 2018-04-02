
BOUT_TOP = ../../../..

SOURCEC         = aiolosmesh.cxx aiolos_interp_to.cxx aiolos_stencils.cxx \
		  aiolos_init.cxx aiolos_derivs.cxx
SOURCEH         = aiolosmesh.hxx aiolos_header.hxx aiolos_init.hxx aiolos_stencils.hxx

TARGET          = lib

include $(BOUT_TOP)/make.config



%xx: %xx.in.py derivs.py \
               tables_cleaned.cxx stencils.py stencils_cleaned.cxx \
               common.py
	@python3 common_test.py &> test.log || ( cat test.log ; exit 1 )
	@python3 test_interp.py &> test.log || ( cat test.log ; exit 1 )
	@echo "  Generating $@"
	@python3 $< > $@.tmp
	@clang-format -i $@.tmp || echo "  Formatting failed - do not git commit this file"
	@mv $@.tmp $@

# Explicitly list cxx files, to prevent make from removing
# intermediate targets.

aiolos_init.o: aiolos_init.hxx aiolos_init.cxx

aiolos_derivs.o: aiolos_init.hxx aiolos_stencils.hxx aiolos_derivs.cxx

aiolosmesh.o: aiolos_header.hxx aiolosmesh.hxx

aiolos_interp_to.o: aiolos_interp_to.cxx

aiolos_stencils.o: aiolos_stencils.cxx

tables_cleaned.cxx:../../index_derivs.cxx
	@grep ' * Lookup tables of functions. Map between names, codes and functions' -B 1 -A 10000 $< |\
	grep ' * Routines to use the above tables to map between function codes, names' -A 2 -B 10000 |\
	gcc -E - > $@
	@echo "  Creating $@"


stencils_cleaned.cxx:../../index_derivs.cxx
	@grep ' * All expect to have an input g' -B 2 -A 10000 $< | \
	grep ' * Lookup tables of functions. Map between names, codes and functions' -A 1 -B 10000 |\
	cat  |\
	gcc -E - > $@
	@echo "  Creating $@"


clean:: clean-generated

clean-generated:
	@rm -f stencils_cleaned.cxx tables_cleaned.cxx *.pyc aiolos_*xx
	@rm -rf __pycache__
