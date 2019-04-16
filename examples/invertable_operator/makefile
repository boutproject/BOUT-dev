BOUT_TOP	= ../..

SOURCEC		= invertable_operator.cxx

target: requires_petsc invertable_operator

requires_petsc:
	@test $$($(BOUT_TOP)/bin/bout-config --has-petsc ) = yes || \
	 (echo "ERROR: This example requires PETSc. Please compile BOUT++ with PETSc." ; exit 1)

include $(BOUT_TOP)/make.config
