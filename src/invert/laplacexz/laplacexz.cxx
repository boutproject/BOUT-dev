#include "impls/cyclic/laplacexz-cyclic.hxx"
#include "impls/petsc/laplacexz-petsc.hxx"

#include <bout/invert/laplacexz.hxx>

std::unique_ptr<LaplaceXZ> LaplaceXZ::create(Mesh *m, Options *options, CELL_LOC loc) {
  return LaplaceXZFactory::getInstance().create(m, options, loc);
}
