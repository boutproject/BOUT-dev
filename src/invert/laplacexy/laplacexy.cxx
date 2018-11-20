
#ifdef BOUT_HAS_PETSC

#include <petscksp.h>

#include <bout/invert/laplacexy.hxx>

#include <bout/assert.hxx>

#include <bout/boutcomm.hxx>
#include <bout/utils.hxx>
#include <bout/sys/timer.hxx>

#include <bout/output.hxx>

#undef __FUNCT__
#define __FUNCT__ "laplacePCapply"
static PetscErrorCode laplacePCapply(PC pc,Vec x,Vec y) {
  int ierr;
  
  // Get the context
  LaplaceXY *s;
  ierr = PCShellGetContext(pc,(void**)&s);CHKERRQ(ierr);
  
  PetscFunctionReturn(s->precon(x, y));
}

LaplaceXY::LaplaceXY(Mesh *m, Options *opt, const CELL_LOC loc) : mesh(m), location(loc) {
  Timer timer("invert");

  if (opt == nullptr) {
    // If no options supplied, use default
    opt = Options::getRoot()->getSection("laplacexy");
  }
  
  // Get MPI communicator
  MPI_Comm comm = BoutComm::get();
  
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
  
  //////////////////////////////////////////////////
  // Specify local indices. This creates a mapping
  // from local indices to index, using a Field2D object

  indexXY = -1;  // Set all points to -1, indicating out of domain

  int ind = 0;
  
  // Y boundaries
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    indexXY(it.ind, mesh->ystart-1) = ind++;
  }
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    indexXY(it.ind, mesh->yend+1) = ind++;
  }
  
  xstart = mesh->xstart;
  if(mesh->firstX())
    xstart -= 1; // Include X guard cells
  xend = mesh->xend;
  if(mesh->lastX())
    xend += 1;
  for(int x=xstart;x<=xend;x++)
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      indexXY(x,y) = ind++;
    }
  
  ASSERT1(ind == localN); // Reached end of range

  //////////////////////////////////////////////////
  // Allocate storage for preconditioner
  
  nloc    = xend - xstart + 1; // Number of X points on this processor
  nsys = mesh->yend - mesh->ystart + 1; // Number of separate Y slices

  acoef = Matrix<BoutReal>(nsys, nloc);
  bcoef = Matrix<BoutReal>(nsys, nloc);
  ccoef = Matrix<BoutReal>(nsys, nloc);
  xvals = Matrix<BoutReal>(nsys, nloc);
  bvals = Matrix<BoutReal>(nsys, nloc);

  // Create a cyclic reduction object
  // FIXME: replace with make_unique when we upgrade to C++14 or add our own version
  cr = std::unique_ptr<CyclicReduce<BoutReal>>(
      new CyclicReduce<BoutReal>(mesh->getXcomm(), nloc));

  //////////////////////////////////////////////////
  // Pre-allocate PETSc storage
  
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
      int localIndex = indexXY(mesh->xstart-1,y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = indexXY(mesh->xstart,y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] -= 1;
      o_nnz[localIndex] += 1;
    }
  }
  if(mesh->lastX()) {
    // Upper X boundary
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = indexXY(mesh->xend+1,y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int localIndex = indexXY(mesh->xend,y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] -= 1;
      o_nnz[localIndex] += 1;
    }
  }
  // Y boundaries
  
  for(int x=mesh->xstart; x <=mesh->xend; x++) {
    // Default to no boundary
    // NOTE: This assumes that communications in Y are to other
    //   processors. If Y is communicated with this processor (e.g. NYPE=1)
    //   then this will result in PETSc warnings about out of range allocations
    
    int localIndex = indexXY(x, mesh->ystart);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    //d_nnz[localIndex] -= 1;  // Note: Slightly inefficient
    o_nnz[localIndex] += 1;
    
    localIndex = indexXY(x, mesh->yend);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    //d_nnz[localIndex] -= 1; // Note: Slightly inefficient
    o_nnz[localIndex] += 1;
  }
  
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    int localIndex = indexXY(it.ind, mesh->ystart-1);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    d_nnz[localIndex] = 2; // Diagonal sub-matrix
    o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    
    localIndex = indexXY(it.ind, mesh->ystart);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    d_nnz[localIndex] += 1;
    o_nnz[localIndex] -= 1;
  }
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    int localIndex = indexXY(it.ind, mesh->yend+1);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    d_nnz[localIndex] = 2; // Diagonal sub-matrix
    o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    
    localIndex = indexXY(it.ind, mesh->yend);
    ASSERT1( (localIndex >= 0) && (localIndex < localN) );
    d_nnz[localIndex] += 1;
    o_nnz[localIndex] -= 1;
  }
  // Pre-allocate
  MatMPIAIJSetPreallocation( MatA, 0, d_nnz, 0, o_nnz );
  MatSetUp(MatA); 
  
  PetscFree( d_nnz );
  PetscFree( o_nnz );
  
  // Determine which row/columns of the matrix are locally owned
  int Istart, Iend;
  MatGetOwnershipRange( MatA, &Istart, &Iend );
  
  // Convert indexXY from local index to global index
  indexXY += Istart;
  
  // Now communicate to fill guard cells
  mesh->communicate(indexXY);

  //////////////////////////////////////////////////
  // Set up KSP
  
  // Declare KSP Context 
  KSPCreate( comm, &ksp ); 
  
  // Configure Linear Solver
  
  bool direct;
  OPTION(opt, direct, false);
  
  if(direct) {
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
#if PETSC_VERSION_GE(3,9,0)
    PCFactorSetMatSolverType(pc,"mumps");
#else
    PCFactorSetMatSolverPackage(pc,"mumps");
#endif
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
    std::string ksptype;
    opt->get("ksptype", ksptype, "gmres");
    
    // Get PC type
    std::string pctype;
    opt->get("pctype", pctype, "none", true);

    KSPSetType( ksp, ksptype.c_str() );
    KSPSetTolerances( ksp, rtol, atol, dtol, maxits );
    
    KSPSetInitialGuessNonzero( ksp, (PetscBool) true );
    
    KSPGetPC(ksp,&pc);
    PCSetType(pc, pctype.c_str());

    if(pctype == "shell") {
      // Using tridiagonal solver as preconditioner
      PCShellSetApply(pc,laplacePCapply);
      PCShellSetContext(pc,this);
      
      bool rightprec;
      OPTION(opt, rightprec, true);
      if(rightprec) {
        KSPSetPCSide(ksp, PC_RIGHT); // Right preconditioning
      }else
        KSPSetPCSide(ksp, PC_LEFT);  // Left preconditioning
    }
  }
  
  KSPSetFromOptions( ksp );

  ///////////////////////////////////////////////////
  // Decide boundary condititions
  if(mesh->periodicY(mesh->xstart)) {
    // Periodic in Y, so in the core
    opt->get("core_bndry_dirichlet", x_inner_dirichlet, false);
  }else {
    // Non-periodic, so in the PF region
    opt->get("pf_bndry_dirichlet", x_inner_dirichlet, true);
  }
  opt->get("y_bndry_dirichlet", y_bndry_dirichlet, false);

  ///////////////////////////////////////////////////
  // Including Y derivatives?

  OPTION(opt, include_y_derivs, true);
  
  ///////////////////////////////////////////////////
  // Set the default coefficients
  setCoefs(1.0, 0.0);
}

void LaplaceXY::setCoefs(const Field2D &A, const Field2D &B) {
  Timer timer("invert");

  Coordinates *coords = mesh->getCoordinates(location);
  
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
      BoutReal J = 0.5*(coords->J(x,y) + coords->J(x+1,y));
      BoutReal g11 = 0.5*(coords->g11(x,y) + coords->g11(x+1,y));
      BoutReal dx = 0.5*(coords->dx(x,y) + coords->dx(x+1,y));
      BoutReal Acoef = 0.5*(A(x,y) + A(x+1,y));
      
      BoutReal val = Acoef * J * g11 / (coords->J(x,y) * dx * coords->dx(x,y));
      xp = val;
      c  = -val;
      
      // Metrics on x-1/2 boundary
      J = 0.5*(coords->J(x,y) + coords->J(x-1,y));
      g11 = 0.5*(coords->g11(x,y) + coords->g11(x-1,y));
      dx = 0.5*(coords->dx(x,y) + coords->dx(x-1,y));
      Acoef = 0.5*(A(x,y) + A(x-1,y));
      
      val = Acoef * J * g11 / (coords->J(x,y) * dx * coords->dx(x,y));
      xm = val;
      c  -= val;

      c += B(x,y);
      
      // Put values into the preconditioner, X derivatives only
      acoef(y - mesh->ystart, x - xstart) = xm;
      bcoef(y - mesh->ystart, x - xstart) = c;
      ccoef(y - mesh->ystart, x - xstart) = xp;

      if( include_y_derivs ) {
        // YY component
        // Metrics at y+1/2
        J = 0.5*(coords->J(x,y) + coords->J(x,y+1));
        BoutReal g_22 = 0.5*(coords->g_22(x,y) + coords->g_22(x,y+1));
        BoutReal g23  = 0.5*(coords->g23(x,y) + coords->g23(x,y+1));
        BoutReal g_23 = 0.5*(coords->g_23(x,y) + coords->g_23(x,y+1));
        BoutReal dy   = 0.5*(coords->dy(x,y) + coords->dy(x,y+1));
        Acoef = 0.5*(A(x,y+1) + A(x,y));
        
        val = -Acoef * J * g23 * g_23 / (g_22 * coords->J(x,y) * dy * coords->dy(x,y));
        yp = val;
        c -= val;
        
        // Metrics at y-1/2
        J    = 0.5*(coords->J(x,y)    + coords->J(x,y-1));
        g_22 = 0.5*(coords->g_22(x,y) + coords->g_22(x,y-1));
        g23  = 0.5*(coords->g23(x,y)  + coords->g23(x,y-1));
        g_23 = 0.5*(coords->g_23(x,y) + coords->g_23(x,y-1));
        dy   = 0.5*(coords->dy(x,y)   + coords->dy(x,y-1));
        Acoef = 0.5*(A(x,y-1) + A(x,y));
        
        val = -Acoef * J * g23 * g_23 / (g_22 * coords->J(x,y) * dy * coords->dy(x,y));
        ym = val;
        c -= val;
      }
      
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
      
      if( include_y_derivs ) {
        // Y + 1
        col = globalIndex(x, y+1);
        MatSetValues(MatA,1,&row,1,&col,&yp,INSERT_VALUES);
        
        // Y - 1
        col = globalIndex(x, y-1);
        MatSetValues(MatA,1,&row,1,&col,&ym,INSERT_VALUES);
      }
      
    }
  }
  
  // X boundaries
  if(mesh->firstX()) {
    if(x_inner_dirichlet) {

      // Dirichlet on inner X boundary
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int row = globalIndex(mesh->xstart-1,y);
        PetscScalar val = 0.5;
        MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
        
        int col = globalIndex(mesh->xstart,y);
        MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        
        // Preconditioner
        bcoef(y - mesh->ystart, 0) = 0.5;
        ccoef(y - mesh->ystart, 0) = 0.5;
      }
      
    }else {
      
      // Neumann on inner X boundary
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int row = globalIndex(mesh->xstart-1,y);
        PetscScalar val = 1.0;
        MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
        
        int col = globalIndex(mesh->xstart,y);
        val = -1.0;
        MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        
        // Preconditioner
        bcoef(y - mesh->ystart, 0) = 1.0;
        ccoef(y - mesh->ystart, 0) = -1.0;
      }
    }
  }
  if(mesh->lastX()) {
    // Dirichlet on outer X boundary
    
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int row = globalIndex(mesh->xend+1,y);
      PetscScalar val = 0.5;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      int col = globalIndex(mesh->xend,y);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
      
      // Preconditioner
      acoef(y - mesh->ystart, mesh->xend + 1 - xstart) = 0.5;
      bcoef(y - mesh->ystart, mesh->xend + 1 - xstart) = 0.5;
    }
  }

  if(y_bndry_dirichlet) {
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
  }else {
    // Neumann on Y boundaries
    for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
      int row = globalIndex(it.ind, mesh->ystart-1);
      PetscScalar val = 1.0;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
    
      val = -1.0;
      int col = globalIndex(it.ind, mesh->ystart);
      
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
    
    for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
      int row = globalIndex(it.ind, mesh->yend+1);
      PetscScalar val = 1.0;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      val = -1.0;
      int col = globalIndex(it.ind, mesh->yend);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
  }
  
  // Assemble Matrix
  MatAssemblyBegin( MatA, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( MatA, MAT_FINAL_ASSEMBLY );

  // Set the operator
#if PETSC_VERSION_GE(3,5,0)
  KSPSetOperators( ksp,MatA,MatA );
#else
  KSPSetOperators( ksp,MatA,MatA,DIFFERENT_NONZERO_PATTERN );
#endif
  
  // Set coefficients for preconditioner
  cr->setCoefs(acoef, bcoef, ccoef);
}

LaplaceXY::~LaplaceXY() {
  PetscBool is_finalised;
  PetscFinalized(&is_finalised);

  if (!is_finalised) {
    // PetscFinalize may already have destroyed this object
    KSPDestroy(&ksp);
  }

  VecDestroy(&xs);
  VecDestroy(&bs);
  MatDestroy(&MatA);
}

const Field2D LaplaceXY::solve(const Field2D &rhs, const Field2D &x0) {
  Timer timer("invert");
  
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

  if(mesh->firstX()) {
    if(x_inner_dirichlet) {
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int ind = globalIndex(mesh->xstart-1,y);
      
        PetscScalar val = x0(mesh->xstart-1,y);
        VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
        
        val = 0.5*(x0(mesh->xstart-1,y) + x0(mesh->xstart,y));
        VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
      }
    }else {
      // Inner X boundary (Neumann)
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        int ind = globalIndex(mesh->xstart-1,y);
        
        PetscScalar val = x0(mesh->xstart-1,y);
        VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
        
        val = 0.0; //x0(mesh->xstart-1,y) - x0(mesh->xstart,y);
        VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
      }
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

  if(y_bndry_dirichlet) {
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
  } else {
    // Y boundaries Neumann
    for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, mesh->ystart-1);
    
      PetscScalar val = x0(it.ind,mesh->ystart-1);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      
      val = 0.0;
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  
    for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, mesh->yend+1);
      
      PetscScalar val = x0(it.ind,mesh->yend+1);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
    
      val = 0.0;
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
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
  result.setLocation(rhs.getLocation());
  
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
      for(int x=mesh->xstart-1; x >= 0; x--)
        result(x,y) = val;
    }
  }

  // Outer X boundary
  if(mesh->lastX()) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      int ind = globalIndex(mesh->xend+1,y);
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      for(int x=mesh->xend+1;x < mesh->LocalNx;x++)
        result(x,y) = val;
    }
  }
  
  // Lower Y boundary
  for(RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, mesh->ystart-1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val );
    for(int y=mesh->ystart-1;y>=0;y--)
      result(it.ind, y) = val;
  }
  
  // Upper Y boundary
  for(RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, mesh->yend+1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val );
    for(int y=mesh->yend+1;y<mesh->LocalNy;y++)
      result(it.ind, y) = val;
  }
  
  return result;
}

/*! Preconditioner
 * NOTE: For efficiency, this routine does not use globalIndex() 
 * in the inner loop. Instead, the indexing must be ordered in
 * exactly the same way as in the construction of indexXY
 */
int LaplaceXY::precon(Vec input, Vec result) {
  
  // Starting index
  int ind = -1;
  
  RangeIterator itdwn=mesh->iterateBndryLowerY();
  if(!itdwn.isDone()) {
    ind = globalIndex(itdwn.ind, mesh->ystart-1);
    
    for(; !itdwn.isDone(); itdwn++) {
      PetscScalar val;
      VecGetValues(input, 1, &ind, &val ); 
      VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
      ind++;
    }
  }
  RangeIterator itup=mesh->iterateBndryUpperY();
  if(!itup.isDone()) {
    if(ind == -1) {
      // No lower boundary
      ind = globalIndex(itup.ind, mesh->yend+1);
    }
    for(; !itup.isDone(); itup++) {
      PetscScalar val;
      VecGetValues(input, 1, &ind, &val ); 
      VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
      ind++;
    }
  }
  if(ind == -1) {
    // No Y boundaries
    ind = globalIndex(xstart, mesh->ystart);
  }
    
  int ind0 = ind;
  // Load vector x into bvals array
  for(int x=xstart;x<=xend;x++) {
    for(int y=mesh->ystart; y<=mesh->yend;y++) {
      PetscScalar val;
      VecGetValues(input, 1, &ind, &val );
      bvals(y - mesh->ystart, x - xstart) = val;
      ind++;
    }
  }
  
  // Solve tridiagonal systems using CR solver
  cr->solve(bvals, xvals);

  // Save result xvals into y array
  ind = ind0;
  for(int x=xstart;x<=xend;x++) {
    for(int y=mesh->ystart; y<=mesh->yend;y++) {
      PetscScalar val = xvals(y - mesh->ystart, x - xstart);
      VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
      ind++;
    }
  }
  VecAssemblyBegin(result);
  VecAssemblyEnd(result);
  return 0;
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
  if( (x < 0) || (x >= mesh->LocalNx) ||
      (y < 0) || (y >= mesh->LocalNy) )
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
