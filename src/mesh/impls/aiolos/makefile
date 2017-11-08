
BOUT_TOP = ../../../..

SOURCEC         = aiolosmesh.cxx
SOURCEH         = aiolosmesh.hxx
TARGET          = lib

include $(BOUT_TOP)/make.config

generated=generated_header.hxx generated_derivs.cxx generated_stencils.cxx generated_init.cxx
#lib: makefile $(BOUT_TOP)/make.config $(BOUT_TOP)/include $(BOUT_TOP)/lib $(OBJ) generated_stencils.cxx

aiolosmesh.cxx: generated.stamp

tables_cleaned.cxx:../../index_derivs.cxx
	@ grep ' * Lookup tables of functions. Map between names, codes and functions' -B 1 -A 10000 $< |\
	grep ' * Routines to use the above tables to map between function codes, names' -A 2 -B 10000 |\
	gcc -E - > $@


stencils_cleaned.cxx:../../index_derivs.cxx
	@ grep ' * All expect to have an input g' -B 2 -A 10000 $< | \
	grep ' * Lookup tables of functions. Map between names, codes and functions' -A 1 -B 10000 |\
	cat  |\
	gcc -E - > $@

generated.stamp:gen_derivs.py \
		tables_cleaned.cxx gen_stencils.py stencils_cleaned.cxx \
		common.py
	python3 $< > generated_derivs.cxx || ( rm $@ ; exit 1)
	@(clang-format -i $(generated) || astyle -n $(generated) ) || exit 0 # ignore failure

# Solution from automake
# https://www.gnu.org/savannah-checkouts/gnu/automake/manual/html_node/Multiple-Outputs.html
$(generated): generated.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  trap 'rm -rf generated.lock generated.stamp' 1 2 13 15; \
## mkdir is a portable test-and-set
	  if mkdir generated.lock 2>/dev/null; then \
## This code is being executed by the first process.
	    rm -f generated.stamp; \
	    $(MAKE) generated.stamp; \
	    result=$$?; rm -rf generated.lock; exit $$result; \
	  else \
## This code is being executed by the follower processes.
## Wait until the first process is done.
	    while test -d generated.lock; do sleep 1; done; \
## Succeed if and only if the first process succeeded.
	    test -f generated.stamp; \
	  fi; \
	fi

clean:: clean-generated

clean-generated:
	@rm -f $(generated)
