
#ifdef BOUT_HAS_PETSC

#include "petsc_laplace3d.hxx"

Laplace3DPetsc::Laplace3DPetsc(Options *opt) : Laplace3D(opt) {
  // Get options
}

Laplace3DPetsc::~Laplace3DPetsc() {

}

const Field3D Laplace3DPetsc::solve(const Field3D &b) {
  // Set matrix coefficients
  
  // Solve
  
}

#endif // BOUT_HAS_PETSC
