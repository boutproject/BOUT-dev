
#include <globals.hxx>
#include "petsc_laplace.hxx"

LaplacePetsc::LaplacePetsc(Options *opt = NULL) : Laplacian(opt) {
  
  // Create a vector using X communicator
  VecCreate(mesh->getXcomm(), &x);
  
  int localN, globalN;
  

  // Set size of vector
  VecSetSizes(x, localN, globalN);
  
  // Set vector options
  VecSetFromOptions(x);
  
  // Duplicate vector
  VecDuplicate(x,&b);
  
  
}

const FieldPerp LaplacePetsc::solve(const FieldPerp &b) {
  return solve(b,b);
}

const FieldPerp LaplacePetsc::solve(const FieldPerp &b, const FieldPerp &x0) {
  
}
