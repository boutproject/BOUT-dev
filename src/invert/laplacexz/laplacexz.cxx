#include "impls/cyclic/laplacexz-cyclic.hxx"
#include "impls/petsc/laplacexz-petsc.hxx"

#include <bout/invert/laplacexz.hxx>

// DO NOT REMOVE: ensures linker keeps all symbols in this TU
void LaplaceXZFactory::ensureRegistered() {}

constexpr decltype(LaplaceXZFactory::type_name) LaplaceXZFactory::type_name;
constexpr decltype(LaplaceXZFactory::section_name) LaplaceXZFactory::section_name;
constexpr decltype(LaplaceXZFactory::option_name) LaplaceXZFactory::option_name;
constexpr decltype(LaplaceXZFactory::default_type) LaplaceXZFactory::default_type;
