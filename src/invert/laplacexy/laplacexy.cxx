// NOLINTBEGIN(misc-include-cleaner, unused-includes)
#include "impls/hypre/laplacexy-hypre.hxx"
#include "impls/petsc/laplacexy-petsc.hxx"
#include "impls/petsc2/laplacexy-petsc2.hxx"
// NOLINTEND(misc-include-cleaner, unused-includes)

#include <bout/invert/laplacexy.hxx>

// DO NOT REMOVE: ensures linker keeps all symbols in this TU
void LaplaceXYFactory::ensureRegistered() {}
