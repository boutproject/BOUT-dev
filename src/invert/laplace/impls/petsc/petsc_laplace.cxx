
#include <globals.hxx>
#include <boutexception.hxx>
#include "petsc_laplace.hxx"

LaplacePetsc::LaplacePetsc(Options *opt = NULL) : Laplacian(opt) {
  
  // Get communicator for group of processors in X
  comm = mesh->getXcomm();

  // Create a vector using X communicator
  VecCreate(comm, &x);
  
  int localN, globalN;
  
  localN = (mesh->xend - mesh->xstart + 1) * (mesh->ngz-1);
  if(mesh->firstX())
    localN += mesh->xstart * (mesh->ngz-1);
  if(mesh->lastX())
    localN += mesh->xstart * (mesh->ngz-1);
  
  // All reduce -> get globalN
  if(MPI_Allreduce(&localN, &globalN, 1, MPI_DOUBLE, 
                   MPI_SUM, comm) != MPI_SUCCESS)
    throw BoutException("Error in MPI_Allreduce during LaplacePetsc initialisation");
  
  
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
  int y = b.getIndex(); // Get the Y index
  
  // Starting row number
  int row = 0; // From PETSC?
  
  
  // Lower boundary
  if(mesh->firstX()) {
    // X = 0 to mesh->xstart-1
    for(int x=0; x < mesh->xstart-1; x++)
      for(int z=0;z < mesh->ngz-1; z++) {
        // Do something with mesh->g11[x][y], coefficients A[x][y][z]
        
        row++; // Increment row in Petsc matrix
      }
  }

  // Main domain with Laplacian operator
  for(int x=mesh->xstart; x <= mesh->xend; x++)
    for(int z=0;z < mesh->ngz-1; z++) {
      // Do something with mesh->g11[x][y], coefficients A[x][y][z]
      
      row++;
    }
  
  // Upper boundary
  if(mesh->lastX()) {
    // X = mesh->xend+1 to mesh->ngx-1;
    for(int x=mesh->xend+1; x < mesh->ngx; x++)
      for(int z=0;z < mesh->ngz-1; z++) {
        // Do something with mesh->g11[x][y], coefficients A[x][y][z]
        
        row++; // Increment row in Petsc matrix
      }
  }
  
  
}
