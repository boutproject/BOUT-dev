
BOUT_TOP = ../..


DIRS            = impls parallel data interpolation
SOURCEC		= difops.cxx interpolation.cxx mesh.cxx boundary_standard.cxx \
		  boundary_factory.cxx boundary_region.cxx \
		  surfaceiter.cxx coordinates.cxx index_derivs.cxx \
	  	  parallel_boundary_region.cxx parallel_boundary_op.cxx fv_ops.cxx \
		  boundary_nonuniform.cxx
SOURCEH		= $(SOURCEC:%.cxx=%.hxx)
TARGET		= lib

boundary_nonuniform.%: boundary_nonuniform.%.in.py
	@echo "  Generating $@"
	@./$< > $@.tmp || (fail=$$?; echo "touch $@ to ignore failed generation" ; exit $$fail)
	@clang-format -i $@.tmp || echo "Formatting failed"
	@mv $@.tmp $@

boundary_factory.o: boundary_nonuniform.hxx boundary_nonuniform.cxx

include $(BOUT_TOP)/make.config
