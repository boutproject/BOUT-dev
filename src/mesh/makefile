
BOUT_TOP = ../..


DIRS            = impls parallel data interpolation
SOURCEC		= difops.cxx interpolation.cxx mesh.cxx boundary_standard.cxx \
		  boundary_factory.cxx boundary_region.cxx meshfactory.cxx \
		  surfaceiter.cxx coordinates.cxx index_derivs.cxx \
	  	  parallel_boundary_region.cxx parallel_boundary_op.cxx fv_ops.cxx \
		  generator-context.cxx
SOURCEH		= $(SOURCEC:%.cxx=%.hxx)
TARGET		= lib

include $(BOUT_TOP)/make.config
