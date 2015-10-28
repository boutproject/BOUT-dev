#include "petsc.hxx"

LaplaceXZpetsc::LaplaceXZpetsc(Mesh *m, Options *options)
  : mesh(m) {

  if(opt == NULL) {
    // If no options supplied, use default
    opt = Options::getRoot()->getSection("laplacexz");
  }

  // Get MPI communicator
  MPI_Comm comm = communicator();

  // Local size
  int localN = localSize();

  // Create Vectors 
  VecCreate( comm, &xs );
  VecSetSizes( xs, localN, PETSC_DETERMINE );
  VecSetFromOptions( xs );
  VecDuplicate( xs , &bs );

  // Set size of Matrix on each processor to localN x localN
  MatCreate( comm, &MatA );                                
  MatSetSizes( MatA, localN, localN, PETSC_DETERMINE, PETSC_DETERMINE );
  MatSetFromOptions(MatA);

  
}

LaplaceXZpetsc::~LaplaceXZpetsc() {
  
}

void LaplaceXZpetsc::setCoefs(const Field2D &A, const Field2D &B) {

}

Field3D LaplaceXZpetsc::solve(const Field3D &b, const Field3D &x0) {
  
}
