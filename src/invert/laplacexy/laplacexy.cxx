
#ifdef BOUT_HAS_PETSC

#include <petscksp.h>

#include <bout/invert/laplacexy.hxx>

#include <bout/assert.hxx>

LaplaceXY::LaplaceXY(Mesh *m, Options *opt) : mesh(m) {
  
  if(opt == NULL) {
    // If no options supplied, use default
    opt = Options::getRoot()->getSection("laplacexy");
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

  // Determine which row/columns of the matrix are locally owned
  int Istart, Iend;
  MatGetOwnershipRange( MatA, &Istart, &Iend );
  
  //////////////////////////////////////////////////
  // Specify global indices. This creates a mapping
  // from local indices to global index, using a Field2D object

  indexXY = -1;  // Set all points to -1, indicating out of domain

  int ind = Istart;
  
  // Y boundaries
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    
    indexXY(it.ind, mesh->ystart-1) = ind++;
  }
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    indexXY(it.ind, mesh->yend+1) = ind++;
  }
  
  int xstart = mesh->xstart;
  if(mesh->firstX())
    xstart -= 1; // Include X guard cells
  int xend = mesh->xend;
  if(mesh->lastX())
    xend += 1;
  for(int x=xstart;x<=xend;x++)
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      indexXY(x,y) = ind++;
    }
  
  ASSERT0(ind == Iend+1); // Reached end of range
  
  // Now communicate to fill guard cells
  mesh->communicate(indexXY);
  
  //////////////////////////////////////////////////
  // Pre-allocate storage
  
  PetscInt *d_nnz, *o_nnz;
  PetscMalloc( (localN)*sizeof(PetscInt), &d_nnz );
  PetscMalloc( (localN)*sizeof(PetscInt), &o_nnz );
  
  for(int i=0;i<localN;i++) {
    // Non-zero elements on this processor
    d_nnz[i] = 5; // Star pattern in 2D
    // Non-zero elements on neighboring processor
    o_nnz[i] = 0;
  }
  
  // X boundaries
  if(mesh->firstX()) {
    // Lower X boundary
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = globalIndex(mesh->xstart-1,y) - Istart;
    
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = globalIndex(mesh->xstart,y) - Istart;
      d_nnz[localIndex] -= 1;
      o_nnz[localIndex] += 1;
    }
  }
  if(mesh->lastX()) {
    // Upper X boundary
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = globalIndex(mesh->xend+1,y) - Istart;
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = globalIndex(mesh->xend,y) - Istart;
      d_nnz[localIndex] -= 1;
      o_nnz[localIndex] += 1;
    }
  }
  // Y boundaries
  for(int x=mesh->xstart; x <=mesh->xend; x++) {
    // Default to no boundary
    int localIndex = globalIndex(x, mesh->ystart) - Istart;
    d_nnz[localIndex] -= 1;
    o_nnz[localIndex] += 1;
    
    localIndex = globalIndex(x, mesh->yend) - Istart;
    d_nnz[localIndex] -= 1;
    o_nnz[localIndex] += 1;
  }
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    int localIndex = globalIndex(it.ind, mesh->ystart-1) - Istart;
    d_nnz[localIndex] = 2; // Diagonal sub-matrix
    o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    
    localIndex = globalIndex(it.ind, mesh->ystart) - Istart;
    d_nnz[localIndex] += 1;
    o_nnz[localIndex] -= 1;
  }
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    int localIndex = globalIndex(it.ind, mesh->yend+1) - Istart;
    d_nnz[localIndex] = 2; // Diagonal sub-matrix
    o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    
    localIndex = globalIndex(it.ind, mesh->yend) - Istart;
    d_nnz[localIndex] += 1;
    o_nnz[localIndex] -= 1;
  }
  
  // Pre-allocate
  MatMPIAIJSetPreallocation( MatA, 0, d_nnz, 0, o_nnz );
  MatSetUp(MatA); 
  
  PetscFree( d_nnz );
  PetscFree( o_nnz );
  
  //////////////////////////////////////////////////
  // Set Matrix elements
  // 
  // (1/J) d/dx ( J * g11 d/dx ) + (1/J) d/dy ( J * g22 d/dy )
  
  for(int x=mesh->xstart; x <= mesh->xend; x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      // stencil entries
      PetscScalar c, xm, xp, ym, yp;
      
      // XX component
      
      // Metrics on x+1/2 boundary
      BoutReal J = 0.5*(mesh->J(x,y) + mesh->J(x+1,y));
      BoutReal g11 = 0.5*(mesh->g11(x,y) + mesh->g11(x+1,y));
      BoutReal dx = 0.5*(mesh->dx(x,y) + mesh->dx(x+1,y));
      
      BoutReal val = J * g11 / (mesh->J(x,y) * dx * mesh->dx(x,y));
      xp = val;
      c  = -val;
      
      // Metrics on x-1/2 boundary
      J = 0.5*(mesh->J(x,y) + mesh->J(x-1,y));
      g11 = 0.5*(mesh->g11(x,y) + mesh->g11(x-1,y));
      dx = 0.5*(mesh->dx(x,y) + mesh->dx(x-1,y));
      
      val = J * g11 / (mesh->J(x,y) * dx * mesh->dx(x,y));
      xm = val;
      c  -= val;
      
      // YY component
      
      // Metrics at y+1/2
      J = 0.5*(mesh->J(x,y) + mesh->J(x,y+1));
      BoutReal g22 = 0.5*(mesh->g22(x,y) + mesh->g22(x,y+1));
      BoutReal dy = 0.5*(mesh->dy(x,y) + mesh->dy(x,y+1));
      
      val = J * g22 / (mesh->J(x,y) * dy * mesh->dy(x,y));
      yp = val;
      c -= val;
      
      // Metrics at y-1/2
      J = 0.5*(mesh->J(x,y) + mesh->J(x,y-1));
      g22 = 0.5*(mesh->g22(x,y) + mesh->g22(x,y-1));
      dy = 0.5*(mesh->dy(x,y) + mesh->dy(x,y-1));
      
      val = J * g22 / (mesh->J(x,y) * dy * mesh->dy(x,y));
      ym = val;
      c -= val;
      
      /////////////////////////////////////////////////
      // Now have a 5-point stencil for the Laplacian
      
      int row = globalIndex(x,y);
      
      // Set the centre (diagonal)
      MatSetValues(MatA,1,&row,1,&row,&c,INSERT_VALUES);
      
      // X + 1
      int col = globalIndex(x+1, y);
      MatSetValues(MatA,1,&row,1,&col,&xp,INSERT_VALUES);
      
      // X - 1
      col = globalIndex(x-1, y);
      MatSetValues(MatA,1,&row,1,&col,&xm,INSERT_VALUES);
      
      // Y + 1
      col = globalIndex(x, y+1);
      MatSetValues(MatA,1,&row,1,&col,&yp,INSERT_VALUES);
      
      // Y - 1
      col = globalIndex(x, y-1);
      MatSetValues(MatA,1,&row,1,&col,&ym,INSERT_VALUES);
    }
  }
  
  // X boundaries
  if(mesh->firstX()) {
    // Neumann on inner X boundary
    
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int row = globalIndex(mesh->xstart-1,y);
      PetscScalar val = 1.0;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      int col = globalIndex(mesh->xstart,y);
      val = -1.0;
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
  }
  if(mesh->lastX()) {
    // Dirichlet on outer X boundary
    
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int row = globalIndex(mesh->xend+1,y);
      PetscScalar val = 0.5;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      int col = globalIndex(mesh->xstart,y);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
  }
  
  // Dirichlet on Y boundaries
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    int row = globalIndex(it.ind, mesh->ystart-1);
    PetscScalar val = 0.5;
    MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
    
    int col = globalIndex(it.ind, mesh->ystart);
    MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
  }
  
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    int row = globalIndex(it.ind, mesh->yend+1);
    PetscScalar val = 0.5;
    MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
    
    int col = globalIndex(it.ind, mesh->yend);
    MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
  }
  
  // Assemble Matrix
  MatAssemblyBegin( MatA, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( MatA, MAT_FINAL_ASSEMBLY );

  //////////////////////////////////////////////////
  // Set up KSP
  
  // Declare KSP Context 
  KSPCreate( comm, &ksp ); 
  
  // Configure Linear Solver               
  KSPSetOperators( ksp,MatA,MatA,DIFFERENT_NONZERO_PATTERN );
  
  bool direct;
  OPTION(opt, direct, false);
  
  if(direct) {
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
    PCFactorSetMatSolverPackage(pc,"mumps");
  }else {
    
    // Convergence Parameters. Solution is considered converged if |r_k| < max( rtol * |b| , atol )
    // where r_k = b - Ax_k. The solution is considered diverged if |r_k| > dtol * |b|.
    BoutReal rtol, atol, dtol;
    int maxits; ///< Maximum iterations
    
    OPTION(opt, rtol, 1e-5);     // Relative tolerance 
    OPTION(opt, atol, 1e-10);    // Absolute tolerance
    OPTION(opt, dtol, 1e3);      // Diverged threshold
    OPTION(opt, maxits, 100000); // Maximum iterations
    
    // Get KSP Solver Type
    string ksptype;
    opt->get("ksptype", ksptype, "gmres");
    
    // Get PC type
    string pctype;
    opt->get("pctype", pctype, "none", true);

    KSPSetType( ksp, ksptype.c_str() );
    KSPSetTolerances( ksp, rtol, atol, dtol, maxits );
    
    KSPSetInitialGuessNonzero( ksp, (PetscBool) true );
    
    KSPGetPC(ksp,&pc);
    PCSetType(pc, pctype.c_str());
  }
  
  KSPSetFromOptions( ksp );
}

LaplaceXY::~LaplaceXY() {
  KSPDestroy( &ksp );
  VecDestroy( &xs );
  VecDestroy( &bs );
  MatDestroy( &MatA );
}

const Field2D LaplaceXY::solve(const Field2D &rhs, const Field2D &x0) {
  
  // Load initial guess x0 into xs and rhs into bs
  
  for(int x=mesh->xstart;x<= mesh->xend;x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(x,y);
      
      PetscScalar val = x0(x,y);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      
      val = rhs(x,y);
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  }
  
  // Inner X boundary (Neumann)
  if(mesh->firstX()) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(mesh->xstart-1,y);
      
      PetscScalar val = x0(mesh->xstart-1,y);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      
      val = x0(mesh->xstart-1,y) - x0(mesh->xstart,y);
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  }
  
  // Outer X boundary (Dirichlet)
  if(mesh->lastX()) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(mesh->xend+1,y);
      
      PetscScalar val = x0(mesh->xend+1,y);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      
      val = 0.5*(x0(mesh->xend,y) + x0(mesh->xend+1,y));
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  }
  
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, mesh->ystart-1);
    
    PetscScalar val = x0(it.ind,mesh->ystart-1);
    VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
    
    val = 0.5*(x0(it.ind, mesh->ystart-1) + x0(it.ind, mesh->ystart));
    VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
  }
  
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, mesh->yend+1);
    
    PetscScalar val = x0(it.ind,mesh->yend+1);
    VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
    
    val = 0.5*(x0(it.ind, mesh->yend+1) + x0(it.ind, mesh->yend));
    VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
  }
  
  // Assemble RHS Vector
  VecAssemblyBegin(bs);
  VecAssemblyEnd(bs);

  // Assemble Trial Solution Vector
  VecAssemblyBegin(xs);
  VecAssemblyEnd(xs);
  
  // Solve the system
  KSPSolve( ksp, bs, xs );
  
  KSPConvergedReason reason;
  KSPGetConvergedReason( ksp, &reason );
  
  if(reason <= 0) {
    throw BoutException("LaplaceXY failed to converge. Reason %d", reason);
  }
  
  //////////////////////////
  // Copy data into result
  
  Field2D result;
  result.allocate();
  
  for(int x=mesh->xstart;x<= mesh->xend;x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(x,y);
      
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      result(x,y) = val;
    }
  }
  
  // Inner X boundary
  if(mesh->firstX()) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(mesh->xstart-1,y);
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      result(mesh->xstart-1,y) = val;
    }
  }

  // Outer X boundary
  if(mesh->lastX()) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(mesh->xend+1,y);
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      result(mesh->xend+1,y) = val;
    }
  }
  
  // Lower Y boundary
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, mesh->ystart-1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val );
    result(it.ind, mesh->ystart-1) = val;
  }
  
  // Upper Y boundary
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, mesh->yend+1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val );
    result(it.ind, mesh->yend+1) = val;
  }
  
  return result;
}

///////////////////////////////////////////////////////////////

int LaplaceXY::localSize() {
  
  // Bulk of points
  int nx = mesh->xend - mesh->xstart + 1;
  int ny = mesh->yend - mesh->ystart + 1;
  
  int n = nx * ny;  
  
  // X boundaries
  if(mesh->firstX())
    n += ny;
  if(mesh->lastX())
    n += ny;
  
  // Y boundaries
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    n++;
  }
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    n++;
  }
  
  return n;
}

int LaplaceXY::globalIndex(int x, int y) {
  if( (x < 0) || (x >= mesh->ngx) ||
      (y < 0) || (y >= mesh->ngy) )
    return -1; // Out of range
 
  // Get the index from a Field2D, round to integer
  return roundInt(indexXY(x,y));
}

int LaplaceXY::roundInt(BoutReal f) {
  if(f > 0.0) {
    return (int) (f + 0.5);
  }
  return (int) (f - 0.5);
}

#endif // BOUT_HAS_PETSC
