test-compile-examples
=====================

This test ensures that all models in examples can be compiled.

This currently works by finding any makefile located beneath
examples. It skips any directories matching doc in order to skip
documentation.  Building documentation requires often additional
tools, and is thus skipped.  The test is split in two parts, one for
examples without PETSc, and one for examples that require PETSc.

PETSc detection works by checking the makefile whether any target
requires `requires_petsc`.
