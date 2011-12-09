
#include "petsc_laplace.hxx"

const FieldPerp LaplacePetsc::solve(const FieldPerp &b) {
  return solve(b,b);
}

const FieldPerp LaplacePetsc::solve(const FieldPerp &b, const FieldPerp &x0) {
  
}
