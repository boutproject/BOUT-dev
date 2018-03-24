
BOUT_TOP = ../../../..

SOURCEC         = aiolosmesh.cxx
SOURCEH         = aiolosmesh.hxx
TARGET          = lib

include $(BOUT_TOP)/make.config

generated=generated_header.hxx generated_derivs.cxx generated_stencils.cxx generated_init.cxx
#lib: makefile $(BOUT_TOP)/make.config $(BOUT_TOP)/include $(BOUT_TOP)/lib $(OBJ) generated_stencils.cxx

aiolosmesh.o: $(generated)

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

.generated.stamp:gen_derivs.py \
		tables_cleaned.cxx gen_stencils.py stencils_cleaned.cxx \
		common.py
	@echo "  Generating code for aiolos mesh"
# Testing scripts
	@python3 common_test.py &> test.log || ( cat test.log ; exit 1 )
	@python3 test_interp.py &> test.log || ( cat test.log ; exit 1 )
# Generating
	@python3 $< > generated_derivs.cxx || ( rm -f $@ ; exit 1)
	@(clang-format -i $(generated) || astyle -n $(generated) ) || exit 0 # ignore failure
	@touch .generated.stamp

# Solution from automake
# https://www.gnu.org/savannah-checkouts/gnu/automake/manual/html_node/Multiple-Outputs.html
$(generated): .generated.stamp
## Recover from the removal of $@
	@ if test -f $@; then :; else \
	  trap 'rm -rf generated.lock .generated.stamp' 1 2 13 15; \
	  if mkdir generated.lock 2>/dev/null; then \
	    rm -f .generated.stamp; \
	    $(MAKE) .generated.stamp; \
	    result=$$?; rm -rf generated.lock; exit $$result; \
	  else \
	    while test -d generated.lock; do sleep 1; done; \
	    test -f .generated.stamp; \
	  fi; \
	fi

clean:: clean-generated

clean-generated:
	@rm -f $(generated) generated.stamp stencils_cleaned.cxx tables_cleaned.cxx *.pyc
	@rm -rf __pycache__ generated.lock
