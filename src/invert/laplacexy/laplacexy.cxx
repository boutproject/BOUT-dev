
#ifdef BOUT_HAS_PETSC

#include <petscksp.h>

#include <bout/invert/laplacexy.hxx>

#include <bout/assert.hxx>

#include <boutcomm.hxx>
#include <derivs.hxx>
#include <globals.hxx>
#include <utils.hxx>
#include <bout/sys/timer.hxx>

#include <output.hxx>

#include <cmath>

#undef __FUNCT__
#define __FUNCT__ "laplacePCapply"
static PetscErrorCode laplacePCapply(PC pc,Vec x,Vec y) {
  int ierr;
  
  // Get the context
  LaplaceXY *s;
  ierr = PCShellGetContext(pc, reinterpret_cast<void**>(&s)); CHKERRQ(ierr);
  
  PetscFunctionReturn(s->precon(x, y));
}

LaplaceXY::LaplaceXY(Mesh *m, Options *opt, const CELL_LOC loc)
    : lib(opt==nullptr ? &(Options::root()["laplacexy"]) : opt),
      localmesh(m==nullptr ? bout::globals::mesh : m), location(loc), monitor(*this) {
  Timer timer("invert");

  if (opt == nullptr) {
    // If no options supplied, use default
    opt = &(Options::root()["laplacexy"]);
  }
  
  finite_volume = (*opt)["finite_volume"].doc(
      "Use finite volume rather than finite difference discretisation."
      ).withDefault(true);

  ///////////////////////////////////////////////////
  // Boundary condititions options
  if (localmesh->periodicY(localmesh->xstart)) {
    // Periodic in Y, so in the core
    x_inner_dirichlet = (*opt)["core_bndry_dirichlet"].withDefault(false);
  } else {
    // Non-periodic, so in the PF region
    x_inner_dirichlet = (*opt)["pf_bndry_dirichlet"].withDefault(true);
  }
  if ((*opt)["y_bndry_dirichlet"].isSet()) {
    throw BoutException("y_bndry_dirichlet has been deprecated. Use y_bndry=dirichlet "
                        "instead.");
  } else {
    y_bndry = (*opt)["y_bndry"].withDefault("neumann");
  }

  // Check value of y_bndry is a supported option
  if (not(
    y_bndry == "dirichlet"
    or y_bndry == "neumann"
    or y_bndry == "free_o3")) {

    throw BoutException("Unrecognized option '{}' for laplacexy:ybndry", y_bndry);
  }

  if (not finite_volume) {
    // Check we can use corner cells
    if (not localmesh->include_corner_cells) {
      throw BoutException(
          "Finite difference form of LaplaceXY allows non-orthogonal x- and "
          "y-directions, so requires mesh:include_corner_cells=true.");
    }
  }

  // Use name of options section as the default prefix for performance logging variables
  default_prefix = opt->name();

  // Get MPI communicator
  auto comm = BoutComm::get();
  
  // Local size
  const int localN = localSize();

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
  for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    // Should not go into corner cells, LaplaceXY stencil does not include them
    if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
      continue;
    }
    indexXY(it.ind, localmesh->ystart-1) = ind++;
  }
  if ((not finite_volume) and localmesh->hasBndryLowerY()) {
    // Corner boundary cells
    if (localmesh->firstX()) {
      indexXY(localmesh->xstart-1, localmesh->ystart-1) = ind++;
    }
    if (localmesh->lastX()) {
      indexXY(localmesh->xend+1, localmesh->ystart-1) = ind++;
    }
  }
  for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    // Should not go into corner cells, LaplaceXY stencil does not include them
    if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
      continue;
    }
    indexXY(it.ind, localmesh->yend+1) = ind++;
  }
  if ((not finite_volume) and localmesh->hasBndryUpperY()) {
    // Corner boundary cells
    if (localmesh->firstX()) {
      indexXY(localmesh->xstart-1, localmesh->yend+1) = ind++;
    }
    if (localmesh->lastX()) {
      indexXY(localmesh->xend+1, localmesh->yend+1) = ind++;
    }
  }
  
  xstart = localmesh->xstart;
  if(localmesh->firstX())
    xstart -= 1; // Include X guard cells
  xend = localmesh->xend;
  if(localmesh->lastX())
    xend += 1;
  for(int x=xstart;x<=xend;x++)
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      indexXY(x,y) = ind++;
    }
  
  ASSERT1(ind == localN); // Reached end of range

  //////////////////////////////////////////////////
  // Allocate storage for preconditioner

  nloc = xend - xstart + 1;                       // Number of X points on this processor
  nsys = localmesh->yend - localmesh->ystart + 1; // Number of separate Y slices

  acoef.reallocate(nsys, nloc);
  bcoef.reallocate(nsys, nloc);
  ccoef.reallocate(nsys, nloc);
  xvals.reallocate(nsys, nloc);
  bvals.reallocate(nsys, nloc);

  // Create a cyclic reduction object
  cr = bout::utils::make_unique<CyclicReduce<BoutReal>>(localmesh->getXcomm(), nloc);

  //////////////////////////////////////////////////
  // Pre-allocate PETSc storage
  
  PetscInt *d_nnz, *o_nnz;
  PetscMalloc( (localN)*sizeof(PetscInt), &d_nnz );
  PetscMalloc( (localN)*sizeof(PetscInt), &o_nnz );
  
  if (finite_volume) {
    setPreallocationFiniteVolume(d_nnz, o_nnz);
  } else {
    setPreallocationFiniteDifference(d_nnz, o_nnz);
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
  // Note, this includes corner cells if necessary
  localmesh->communicate(indexXY);

  //////////////////////////////////////////////////
  // Set up KSP
  
  // Declare KSP Context 
  KSPCreate( comm, &ksp ); 
  
  // Configure Linear Solver
  
  const bool direct = (*opt)["direct"].doc("Use a direct LU solver").withDefault(false);
  
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

    const BoutReal rtol = (*opt)["rtol"].doc("Relative tolerance").withDefault(1e-5);
    const BoutReal atol = (*opt)["atol"]
            .doc("Absolute tolerance. The solution is considered converged if |Ax-b| "
                 "< max( rtol * |b| , atol )")
            .withDefault(1e-10);
    const BoutReal dtol = (*opt)["dtol"]
                        .doc("The solution is considered diverged if |Ax-b| > dtol * |b|")
                        .withDefault(1e3);
    const int maxits = (*opt)["maxits"].doc("Maximum iterations").withDefault(100000);

    // Get KSP Solver Type
    const std::string ksptype = (*opt)["ksptype"].doc("KSP solver type").withDefault("gmres");
    
    // Get PC type
    const std::string pctype = (*opt)["pctype"].doc("Preconditioner type").withDefault("none");

    KSPSetType( ksp, ksptype.c_str() );
    KSPSetTolerances( ksp, rtol, atol, dtol, maxits );

    KSPSetInitialGuessNonzero(ksp, static_cast<PetscBool>(true));

    KSPGetPC(ksp,&pc);
    PCSetType(pc, pctype.c_str());

    if (pctype == "shell") {
      // Using tridiagonal solver as preconditioner
      PCShellSetApply(pc,laplacePCapply);
      PCShellSetContext(pc,this);
      
      const bool rightprec = (*opt)["rightprec"].doc("Use right preconditioning?").withDefault(true);
      if (rightprec) {
        KSPSetPCSide(ksp, PC_RIGHT); // Right preconditioning
      } else {
        KSPSetPCSide(ksp, PC_LEFT);  // Left preconditioning
      }
    }
  }
  
  lib.setOptionsFromInputFile(ksp);

  ///////////////////////////////////////////////////
  // Including Y derivatives?

  include_y_derivs = (*opt)["include_y_derivs"]
                         .doc("Include Y derivatives in operator to invert?")
                         .withDefault(true);

  ///////////////////////////////////////////////////
  // Set the default coefficients
  Field2D one(1., localmesh);
  Field2D zero(0., localmesh);
  one.setLocation(location);
  zero.setLocation(location);
  setCoefs(one, zero);
}

void LaplaceXY::setPreallocationFiniteVolume(PetscInt* d_nnz, PetscInt* o_nnz) {
  const int localN = localSize();

  // This discretisation uses a 5-point stencil
  for(int i=0;i<localN;i++) {
    // Non-zero elements on this processor
    d_nnz[i] = 5; // Star pattern in 2D
    // Non-zero elements on neighboring processor
    o_nnz[i] = 0;
  }

  // X boundaries
  if(localmesh->firstX()) {
    // Lower X boundary
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      const int localIndex = globalIndex(localmesh->xstart - 1, y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );

      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      const int localIndex = globalIndex(localmesh->xstart, y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] -= 1;
      o_nnz[localIndex] += 1;
    }
  }
  if(localmesh->lastX()) {
    // Upper X boundary
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      const int localIndex = globalIndex(localmesh->xend + 1, y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      const int localIndex = globalIndex(localmesh->xend, y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] -= 1;
      o_nnz[localIndex] += 1;
    }
  }
  // Y boundaries

  for(int x=localmesh->xstart; x <=localmesh->xend; x++) {
    // Default to no boundary
    // NOTE: This assumes that communications in Y are to other
    //   processors. If Y is communicated with this processor (e.g. NYPE=1)
    //   then this will result in PETSc warnings about out of range allocations
    {
      const int localIndex = globalIndex(x, localmesh->ystart);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      // d_nnz[localIndex] -= 1;  // Note: Slightly inefficient
      o_nnz[localIndex] += 1;
    }
    {
      const int localIndex = globalIndex(x, localmesh->yend);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      // d_nnz[localIndex] -= 1; // Note: Slightly inefficient
      o_nnz[localIndex] += 1;
    }
  }

  for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    // Should not go into corner cells, LaplaceXY stencil does not include them
    if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
      continue;
    }
    {
      const int localIndex = globalIndex(it.ind, localmesh->ystart - 1);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
    {
      const int localIndex = globalIndex(it.ind, localmesh->ystart);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      d_nnz[localIndex] += 1;
      o_nnz[localIndex] -= 1;
    }
  }
  for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    // Should not go into corner cells, LaplaceXY stencil does not include them
    if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
      continue;
    }
    {
      const int localIndex = globalIndex(it.ind, localmesh->yend + 1);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
    {
      const int localIndex = globalIndex(it.ind, localmesh->yend);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      d_nnz[localIndex] += 1;
      o_nnz[localIndex] -= 1;
    }
  }
}

void LaplaceXY::setPreallocationFiniteDifference(PetscInt* d_nnz, PetscInt* o_nnz) {
  const int localN = localSize();

  // This discretisation uses a 9-point stencil
  for(int i=0;i<localN;i++) {
    // Non-zero elements on this processor
    d_nnz[i] = 9; // Square pattern in 2D
    // Non-zero elements on neighboring processor
    o_nnz[i] = 0;
  }

  // X boundaries
  if(localmesh->firstX()) {
    // Lower X boundary
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      const int localIndex = globalIndex(localmesh->xstart - 1, y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );

      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      const int localIndex = globalIndex(localmesh->xstart, y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] -= 3;
      o_nnz[localIndex] += 3;
    }
  }
  if(localmesh->lastX()) {
    // Upper X boundary
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      const int localIndex = globalIndex(localmesh->xend + 1, y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] = 2; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
  }else {
    // On another processor
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      const int localIndex = globalIndex(localmesh->xend, y);
      ASSERT1( (localIndex >= 0) && (localIndex < localN) );
      d_nnz[localIndex] -= 3;
      o_nnz[localIndex] += 3;
    }
  }
  // Y boundaries
  const int y_bndry_stencil_size = (y_bndry == "free_o3") ? 4 : 2;

  for(int x=localmesh->xstart; x <=localmesh->xend; x++) {
    // Default to no boundary
    // NOTE: This assumes that communications in Y are to other
    //   processors. If Y is communicated with this processor (e.g. NYPE=1)
    //   then this will result in PETSc warnings about out of range allocations
    {
      const int localIndex = globalIndex(x, localmesh->ystart);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      // d_nnz[localIndex] -= 3;  // Note: Slightly inefficient
      o_nnz[localIndex] += 3;
    }
    {
      const int localIndex = globalIndex(x, localmesh->yend);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      // d_nnz[localIndex] -= 3; // Note: Slightly inefficient
      o_nnz[localIndex] += 3;
    }
  }

  for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    // Should not go into corner cells, they are handled specially below
    if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
      continue;
    }
    {
      const int localIndex = globalIndex(it.ind, localmesh->ystart - 1);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      d_nnz[localIndex] = y_bndry_stencil_size; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
    {
      const int localIndex = globalIndex(it.ind, localmesh->ystart);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      //d_nnz[localIndex] += 3;
      o_nnz[localIndex] -= 3;
    }
  }
  if (localmesh->hasBndryLowerY()) {
    if (y_bndry == "dirichlet") {
      // special handling for the corners, since we use a free_o3 y-boundary
      // condition just in the corners when y_bndry=="dirichlet"
      if (localmesh->firstX()) {
        const int localIndex = globalIndex(localmesh->xstart-1, localmesh->ystart-1);
        ASSERT1((localIndex >= 0) && (localIndex < localN));
        d_nnz[localIndex] = 4;
      }
      if (localmesh->lastX()) {
        const int localIndex = globalIndex(localmesh->xend+1, localmesh->ystart-1);
        ASSERT1((localIndex >= 0) && (localIndex < localN));
        d_nnz[localIndex] = 4;
      }
    } else {
      if (localmesh->firstX()) {
        const int localIndex = globalIndex(localmesh->xstart-1, localmesh->ystart-1);
        ASSERT1((localIndex >= 0) && (localIndex < localN));
        d_nnz[localIndex] = y_bndry_stencil_size;
      }
      if (localmesh->lastX()) {
        const int localIndex = globalIndex(localmesh->xend+1, localmesh->ystart-1);
        ASSERT1((localIndex >= 0) && (localIndex < localN));
        d_nnz[localIndex] = y_bndry_stencil_size;
      }
    }
  }
  for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    // Should not go into corner cells, they are handled specially below
    if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
      continue;
    }
    {
      const int localIndex = globalIndex(it.ind, localmesh->yend + 1);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      d_nnz[localIndex] = y_bndry_stencil_size; // Diagonal sub-matrix
      o_nnz[localIndex] = 0; // Off-diagonal sub-matrix
    }
    {
      const int localIndex = globalIndex(it.ind, localmesh->yend);
      ASSERT1((localIndex >= 0) && (localIndex < localN));
      //d_nnz[localIndex] += 3;
      o_nnz[localIndex] -= 3;
    }
  }
  if (localmesh->hasBndryUpperY()) {
    if (y_bndry == "dirichlet") {
      // special handling for the corners, since we use a free_o3 y-boundary
      // condition just in the corners when y_bndry=="dirichlet"
      if (localmesh->firstX()) {
        const int localIndex = globalIndex(localmesh->xstart-1, localmesh->yend+1);
        ASSERT1((localIndex >= 0) && (localIndex < localN));
        d_nnz[localIndex] = 4;
      }
      if (localmesh->lastX()) {
        const int localIndex = globalIndex(localmesh->xend+1, localmesh->yend+1);
        ASSERT1((localIndex >= 0) && (localIndex < localN));
        d_nnz[localIndex] = 4;
      }
    } else {
      if (localmesh->firstX()) {
        const int localIndex = globalIndex(localmesh->xstart-1, localmesh->yend+1);
        ASSERT1((localIndex >= 0) && (localIndex < localN));
        d_nnz[localIndex] = y_bndry_stencil_size;
      }
      if (localmesh->lastX()) {
        const int localIndex = globalIndex(localmesh->xend+1, localmesh->yend+1);
        ASSERT1((localIndex >= 0) && (localIndex < localN));
        d_nnz[localIndex] = y_bndry_stencil_size;
      }
    }
  }
}

void LaplaceXY::setCoefs(const Field2D &A, const Field2D &B) {
  Timer timer("invert");

  ASSERT1(A.getMesh() == localmesh);
  ASSERT1(B.getMesh() == localmesh);
  ASSERT1(A.getLocation() == location);
  ASSERT1(B.getLocation() == location);

  if (finite_volume) {
    setMatrixElementsFiniteVolume(A, B);
  } else {
    setMatrixElementsFiniteDifference(A, B);
  }
  
  // X boundaries
  if(localmesh->firstX()) {
    if(x_inner_dirichlet) {

      // Dirichlet on inner X boundary
      for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
        int row = globalIndex(localmesh->xstart-1,y);
        PetscScalar val = 0.5;
        MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
        
        int col = globalIndex(localmesh->xstart,y);
        MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        
        // Preconditioner
        bcoef(y - localmesh->ystart, 0) = 0.5;
        ccoef(y - localmesh->ystart, 0) = 0.5;
      }
      
    }else {
      
      // Neumann on inner X boundary
      for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
        int row = globalIndex(localmesh->xstart-1,y);
        PetscScalar val = 1.0;
        MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
        
        int col = globalIndex(localmesh->xstart,y);
        val = -1.0;
        MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        
        // Preconditioner
        bcoef(y - localmesh->ystart, 0) = 1.0;
        ccoef(y - localmesh->ystart, 0) = -1.0;
      }
    }
  }
  if(localmesh->lastX()) {
    // Dirichlet on outer X boundary
    
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      int row = globalIndex(localmesh->xend+1,y);
      PetscScalar val = 0.5;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      int col = globalIndex(localmesh->xend,y);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
      
      // Preconditioner
      acoef(y - localmesh->ystart, localmesh->xend + 1 - xstart) = 0.5;
      bcoef(y - localmesh->ystart, localmesh->xend + 1 - xstart) = 0.5;
    }
  }

  if (y_bndry == "dirichlet") {
    // Dirichlet on Y boundaries
    for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int row = globalIndex(it.ind, localmesh->ystart-1);
      PetscScalar val = 0.5;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      int col = globalIndex(it.ind, localmesh->ystart);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
    
    for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int row = globalIndex(it.ind, localmesh->yend+1);
      PetscScalar val = 0.5;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);
      
      int col = globalIndex(it.ind, localmesh->yend);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
  } else if (y_bndry == "neumann") {
    // Neumann on Y boundaries
    for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int row = globalIndex(it.ind, localmesh->ystart-1);
      PetscScalar val = 1.0;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

      val = -1.0;
      int col = globalIndex(it.ind, localmesh->ystart);

      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
    
    for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int row = globalIndex(it.ind, localmesh->yend+1);
      PetscScalar val = 1.0;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

      val = -1.0;
      int col = globalIndex(it.ind, localmesh->yend);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
  } else if (y_bndry == "free_o3") {
    // 'free_o3' extrapolating boundary condition on Y boundaries
    for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int row = globalIndex(it.ind, localmesh->ystart-1);
      PetscScalar val = 1.0;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

      val = -3.0;
      int col = globalIndex(it.ind, localmesh->ystart);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

      val = 3.0;
      col = globalIndex(it.ind, localmesh->ystart+1);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

      val = -1.0;
      col = globalIndex(it.ind, localmesh->ystart+2);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }

    for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int row = globalIndex(it.ind, localmesh->yend+1);
      PetscScalar val = 1.0;
      MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

      val = -3.0;
      int col = globalIndex(it.ind, localmesh->yend);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

      val = 3.0;
      col = globalIndex(it.ind, localmesh->yend-1);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

      val = -1.0;
      col = globalIndex(it.ind, localmesh->yend-2);
      MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
    }
  } else {
    throw BoutException("Unsupported option for y_bndry");
  }

  if (not finite_volume) {
    // Handle corner boundary cells in case we need to include D2DXDY
    // Apply the y-boundary-condition to the cells in the x-boundary - this is an
    // arbitrary choice, cf. connections around the X-point

    if (localmesh->hasBndryLowerY()) {
      if (localmesh->firstX()) {
        if (y_bndry == "neumann") {
          // Neumann y-bc
          // f(xs-1,ys-1) = f(xs-1,ys)
          PetscScalar val = 1.0;
          int row = globalIndex(localmesh->xstart-1, localmesh->ystart-1);
          MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

          val = -1.0;
          int col = globalIndex(localmesh->xstart-1, localmesh->ystart);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        } else if (y_bndry == "free_o3" or y_bndry == "dirichlet") {
          // 'free_o3' extrapolating boundary condition on Y boundaries
          // f(xs-1,ys-1) = 3*f(xs-1,ys) - 3*f(xs-1,ys+1) + f(xs-1,ys+2)
          //
          // Use free_o3 at the corners for Dirichlet y-boundaries because we don't know
          // what value to pass for the corner
          PetscScalar val = 1.0;
          int row = globalIndex(localmesh->xstart-1, localmesh->ystart-1);
          MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

          val = -3.0;
          int col = globalIndex(localmesh->xstart-1, localmesh->ystart);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

          val = 3.0;
          col = globalIndex(localmesh->xstart-1, localmesh->ystart+1);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

          val = -1.0;
          col = globalIndex(localmesh->xstart-1, localmesh->ystart+2);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        } else {
          throw BoutException("Unsupported option for y_bndry");
        }
      }
      if (localmesh->lastX()) {
        if (y_bndry == "neumann") {
          // Neumann y-bc
          // f(xe+1,ys-1) = f(xe+1,ys)
          PetscScalar val = 1.0;
          int row = globalIndex(localmesh->xend+1, localmesh->ystart-1);
          MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

          val = -1.0;
          int col = globalIndex(localmesh->xend+1, localmesh->ystart);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        } else if (y_bndry == "free_o3" or y_bndry == "dirichlet") {
          // 'free_o3' extrapolating boundary condition on Y boundaries
          // f(xe+1,ys-1) = 3*f(xe+1,ys) - 3*f(xe+1,ys+1) + f(xe+1,ys+2)
          //
          // Use free_o3 at the corners for Dirichlet y-boundaries because we don't know
          // what value to pass for the corner
          PetscScalar val = 1.0;
          int row = globalIndex(localmesh->xend+1, localmesh->ystart-1);
          MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

          val = -3.0;
          int col = globalIndex(localmesh->xend+1, localmesh->ystart);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

          val = 3.0;
          col = globalIndex(localmesh->xend+1, localmesh->ystart+1);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

          val = -1.0;
          col = globalIndex(localmesh->xend+1, localmesh->ystart+2);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        } else {
          throw BoutException("Unsupported option for y_bndry");
        }

      }
    }
    if (localmesh->hasBndryUpperY()) {
      if (localmesh->firstX()) {
        if (y_bndry == "neumann") {
          // Neumann y-bc
          // f(xs-1,ys-1) = f(xs-1,ys)
          PetscScalar val = 1.0;
          int row = globalIndex(localmesh->xstart-1, localmesh->yend+1);
          MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

          val = -1.0;
          int col = globalIndex(localmesh->xstart-1, localmesh->yend);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        } else if (y_bndry == "free_o3" or y_bndry == "dirichlet") {
          // 'free_o3' extrapolating boundary condition on Y boundaries
          // f(xs-1,ys-1) = 3*f(xs-1,ys) - 3*f(xs-1,ys+1) + f(xs-1,ys+2)
          //
          // Use free_o3 at the corners for Dirichlet y-boundaries because we don't know
          // what value to pass for the corner
          PetscScalar val = 1.0;
          int row = globalIndex(localmesh->xstart-1, localmesh->yend+1);
          MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

          val = -3.0;
          int col = globalIndex(localmesh->xstart-1, localmesh->yend);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

          val = 3.0;
          col = globalIndex(localmesh->xstart-1, localmesh->yend-1);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

          val = -1.0;
          col = globalIndex(localmesh->xstart-1, localmesh->yend-2);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        } else {
          throw BoutException("Unsupported option for y_bndry");
        }
      }
      if (localmesh->lastX()) {
        if (y_bndry == "neumann") {
          // Neumann y-bc
          // f(xe+1,ys-1) = f(xe+1,ys)
          PetscScalar val = 1.0;
          int row = globalIndex(localmesh->xend+1, localmesh->yend+1);
          MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

          val = -1.0;
          int col = globalIndex(localmesh->xend+1, localmesh->yend);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        } else if (y_bndry == "free_o3" or y_bndry == "dirichlet") {
          // 'free_o3' extrapolating boundary condition on Y boundaries
          // f(xe+1,ys-1) = 3*f(xe+1,ys) - 3*f(xe+1,ys+1) + f(xe+1,ys+2)
          //
          // Use free_o3 at the corners for Dirichlet y-boundaries because we don't know
          // what value to pass for the corner
          PetscScalar val = 1.0;
          int row = globalIndex(localmesh->xend+1, localmesh->yend+1);
          MatSetValues(MatA,1,&row,1,&row,&val,INSERT_VALUES);

          val = -3.0;
          int col = globalIndex(localmesh->xend+1, localmesh->yend);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

          val = 3.0;
          col = globalIndex(localmesh->xend+1, localmesh->yend-1);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);

          val = -1.0;
          col = globalIndex(localmesh->xend+1, localmesh->yend-2);
          MatSetValues(MatA,1,&row,1,&col,&val,INSERT_VALUES);
        } else {
          throw BoutException("Unsupported option for y_bndry");
        }

      }
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

void LaplaceXY::setMatrixElementsFiniteVolume(const Field2D &A, const Field2D &B) {
  //////////////////////////////////////////////////
  // Set Matrix elements
  //
  // (1/J) d/dx ( J * g11 d/dx ) + (1/J) d/dy ( J * g22 d/dy )

  auto coords = localmesh->getCoordinates(location);

  for(int x=localmesh->xstart; x <= localmesh->xend; x++) {
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
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
      acoef(y - localmesh->ystart, x - xstart) = xm;
      bcoef(y - localmesh->ystart, x - xstart) = c;
      ccoef(y - localmesh->ystart, x - xstart) = xp;

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
}

void LaplaceXY::setMatrixElementsFiniteDifference(const Field2D &A, const Field2D &B) {
  //////////////////////////////////////////////////
  // Set Matrix elements
  //
  // Div(A Grad(f)) + B f
  // = A Laplace_perp(f) + Grad_perp(A).Grad_perp(f) + B f
  // = A*(G1*dfdx + (G2-1/J*d/dy(J/g_22))*dfdy
  //      + g11*d2fdx2 + (g22-1/g_22)*d2fdy2 + 2*g12*d2fdxdy)
  //   + g11*dAdx*dfdx + (g22-1/g_22)*dAdy*dfdy + g12*(dAdx*dfdy + dAdy*dfdx)
  //   + B*f

  auto coords = localmesh->getCoordinates(location);

  Field2D coef_dfdy = coords->G2 - DDY(coords->J/coords->g_22)/coords->J;

  for(int x = localmesh->xstart; x <= localmesh->xend; x++) {
    for(int y = localmesh->ystart; y <= localmesh->yend; y++) {
      // stencil entries
      PetscScalar c, xm, xp, ym, yp, xpyp, xpym, xmyp, xmym;

      BoutReal dx = coords->dx(x,y);

      // A*G1*dfdx
      BoutReal val = A(x, y)*coords->G1(x, y)/(2.*dx);
      xp = val;
      xm = -val;

      // A*g11*d2fdx2
      val = A(x, y)*coords->g11(x, y)/SQ(dx);
      xp += val;
      c = -2.*val;
      xm += val;
      // Non-uniform grid correction
      val = A(x, y)*coords->g11(x, y)*coords->d1_dx(x, y)/(2.*dx);
      xp += val;
      xm -= val;

      // g11*dAdx*dfdx
      val = coords->g11(x, y)*(A(x+1, y) - A(x-1, y))/(4.*SQ(dx));
      xp += val;
      xm -= val;

      // B*f
      c += B(x,y);

      // Put values into the preconditioner, X derivatives only
      acoef(y - localmesh->ystart, x - xstart) = xm;
      bcoef(y - localmesh->ystart, x - xstart) = c;
      ccoef(y - localmesh->ystart, x - xstart) = xp;

      if(include_y_derivs) {
        BoutReal dy = coords->dy(x,y);
        BoutReal dAdx = (A(x+1, y) - A(x-1, y))/(2.*dx);
        BoutReal dAdy = (A(x, y+1) - A(x, y-1))/(2.*dy);

        // A*(G2-1/J*d/dy(J/g_22))*dfdy
        val = A(x, y)*coef_dfdy(x, y)/(2.*dy);
        yp = val;
        ym = -val;

        // A*(g22-1/g_22)*d2fdy2
        val = A(x, y)*(coords->g22(x, y) - 1./coords->g_22(x,y))/SQ(dy);
        yp += val;
        c -= 2.*val;
        ym += val;
        // Non-uniform mesh correction
        val = A(x, y)*(coords->g22(x, y) - 1./coords->g_22(x,y))
          *coords->d1_dy(x, y)/(2.*dy);
        yp += val;
        ym -= val;

        // 2*A*g12*d2dfdxdy
        val = A(x, y)*coords->g12(x, y)/(2.*dx*dy);
        xpyp = val;
        xpym = -val;
        xmyp = -val;
        xmym = val;

        // g22*dAdy*dfdy
        val = (coords->g22(x, y) - 1./coords->g_22(x,y))*dAdy/(2.*dy);
        yp += val;
        ym -= val;

        // g12*(dAdx*dfdy + dAdy*dfdx)
        val = coords->g12(x, y)*dAdx/(2.*dy);
        yp += val;
        ym -= val;
        val = coords->g12(x, y)*dAdy/(2.*dx);
        xp += val;
        xm -= val;
      }

      /////////////////////////////////////////////////
      // Now have a 9-point stencil for the Laplacian

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

        // X + 1, Y + 1
        col = globalIndex(x+1, y+1);
        MatSetValues(MatA,1,&row,1,&col,&xpyp,INSERT_VALUES);

        // X + 1, Y - 1
        col = globalIndex(x+1, y-1);
        MatSetValues(MatA,1,&row,1,&col,&xpym,INSERT_VALUES);

        // X - 1, Y + 1
        col = globalIndex(x-1, y+1);
        MatSetValues(MatA,1,&row,1,&col,&xmyp,INSERT_VALUES);

        // X - 1, Y - 1
        col = globalIndex(x-1, y-1);
        MatSetValues(MatA,1,&row,1,&col,&xmym,INSERT_VALUES);
      }
    }
  }
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

const Field2D LaplaceXY::solve(const Field2D& rhs, const Field2D& x0) {
  Timer timer("invert");

  ASSERT1(rhs.getMesh() == localmesh);
  ASSERT1(x0.getMesh() == localmesh);
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  // Load initial guess x0 into xs and rhs into bs

  for (int x = localmesh->xstart; x <= localmesh->xend; x++) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      int ind = globalIndex(x, y);

      PetscScalar val = x0(x, y);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      val = rhs(x, y);
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }
  }

  if (finite_volume) {
    solveFiniteVolume(x0);
  } else {
    solveFiniteDifference(x0);
  }

  // Assemble RHS Vector
  VecAssemblyBegin(bs);
  VecAssemblyEnd(bs);

  // Assemble Trial Solution Vector
  VecAssemblyBegin(xs);
  VecAssemblyEnd(xs);

  // Solve the system
  KSPSolve(ksp, bs, xs);

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  if (reason <= 0) {
    throw BoutException("LaplaceXY failed to converge. Reason {:d}", reason);
  }

  if (save_performance) {
    // Update performance monitoring information
    n_calls++;

    int iterations = 0;
    KSPGetIterationNumber(ksp, &iterations);

    average_iterations = BoutReal(n_calls - 1) / BoutReal(n_calls) * average_iterations
                         + BoutReal(iterations) / BoutReal(n_calls);
  }

  //////////////////////////
  // Copy data into result

  Field2D result;
  result.allocate();
  result.setLocation(location);

  for (int x = localmesh->xstart; x <= localmesh->xend; x++) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      int ind = globalIndex(x, y);

      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val);
      result(x, y) = val;
    }
  }

  // Inner X boundary
  if (localmesh->firstX()) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      int ind = globalIndex(localmesh->xstart - 1, y);
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val);
      for (int x = localmesh->xstart - 1; x >= 0; x--)
        result(x, y) = val;
    }
  }

  // Outer X boundary
  if (localmesh->lastX()) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      int ind = globalIndex(localmesh->xend + 1, y);
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val);
      for (int x = localmesh->xend + 1; x < localmesh->LocalNx; x++)
        result(x, y) = val;
    }
  }

  // Lower Y boundary
  for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    if (
         // Should not go into corner cells, finite-volume LaplaceXY stencil does not include
         // them
         (finite_volume and (it.ind < localmesh->xstart or it.ind > localmesh->xend))

         // Only go into first corner cell for finite-difference
         or (it.ind < localmesh->xstart - 1 or it.ind > localmesh->xend + 1)
       ) {
      continue;
    }
    int ind = globalIndex(it.ind, localmesh->ystart - 1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val);
    for (int y = localmesh->ystart - 1; y >= 0; y--)
      result(it.ind, y) = val;
  }

  // Upper Y boundary
  for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    if (
         // Should not go into corner cells, finite-volume LaplaceXY stencil does not include
         // them
         (finite_volume and (it.ind < localmesh->xstart or it.ind > localmesh->xend))

         // Only go into first corner cell for finite-difference
         or (it.ind < localmesh->xstart - 1 or it.ind > localmesh->xend + 1)
       ) {
      continue;
    }
    int ind = globalIndex(it.ind, localmesh->yend + 1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val);
    for (int y = localmesh->yend + 1; y < localmesh->LocalNy; y++)
      result(it.ind, y) = val;
  }

  return result;
}

void LaplaceXY::solveFiniteVolume(const Field2D& x0) {
  // Use original LaplaceXY implementation of passing boundary values for backward
  // compatibility
  if (localmesh->firstX()) {
    if (x_inner_dirichlet) {
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        int ind = globalIndex(localmesh->xstart - 1, y);

        PetscScalar val = x0(localmesh->xstart - 1, y);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        val = 0.5 * (x0(localmesh->xstart - 1, y) + x0(localmesh->xstart, y));
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    } else {
      // Inner X boundary (Neumann)
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        int ind = globalIndex(localmesh->xstart - 1, y);

        PetscScalar val = x0(localmesh->xstart - 1, y);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        val = 0.0; // x0(localmesh->xstart-1,y) - x0(localmesh->xstart,y);
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    }
  }

  // Outer X boundary (Dirichlet)
  if (localmesh->lastX()) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      int ind = globalIndex(localmesh->xend + 1, y);

      PetscScalar val = x0(localmesh->xend + 1, y);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      val = 0.5 * (x0(localmesh->xend, y) + x0(localmesh->xend + 1, y));
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }
  }

  if (y_bndry == "dirichlet") {
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->ystart - 1);

      PetscScalar val = x0(it.ind, localmesh->ystart - 1);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      val = 0.5 * (x0(it.ind, localmesh->ystart - 1) + x0(it.ind, localmesh->ystart));
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->yend + 1);

      PetscScalar val = x0(it.ind, localmesh->yend + 1);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      val = 0.5 * (x0(it.ind, localmesh->yend + 1) + x0(it.ind, localmesh->yend));
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }
  } else if (y_bndry == "neumann" or y_bndry == "free_o3") {
    // Y boundaries Neumann
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->ystart - 1);

      PetscScalar val = x0(it.ind, localmesh->ystart - 1);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      val = 0.0;
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      // Should not go into corner cells, LaplaceXY stencil does not include them
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->yend + 1);

      PetscScalar val = x0(it.ind, localmesh->yend + 1);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      val = 0.0;
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }
  } else {
    throw BoutException("Unsupported option for y_bndry");
  }
}

void LaplaceXY::solveFiniteDifference(const Field2D& x0) {
  // For finite-difference implementation pass boundary values in the same way as for
  // the 'Laplacian' solvers - the value to use (for Dirichlet boundary conditions) on
  // the boundary (which is half way between grid cell and boundary cell) is passed as
  // the value in the first boundary cell.
  if (localmesh->firstX()) {
    if (x_inner_dirichlet) {
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        int ind = globalIndex(localmesh->xstart - 1, y);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = 2. * x0(localmesh->xstart - 1, y) - x0(localmesh->xstart, y);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Pass the value from boundary cell of x0 as the boundary condition to the rhs
        val = x0(localmesh->xstart - 1, y);
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    } else {
      // Inner X boundary (Neumann)
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        int ind = globalIndex(localmesh->xstart - 1, y);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = x0(localmesh->xstart, y);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    }
  }

  // Outer X boundary (Dirichlet)
  if (localmesh->lastX()) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      int ind = globalIndex(localmesh->xend + 1, y);

      // For the boundary value of the initial guess, use the value that would be set by
      // applying the boundary condition to the initial guess
      PetscScalar val = 2. * x0(localmesh->xend + 1, y) - x0(localmesh->xend, y);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      // Pass the value from boundary cell of x0 as the boundary condition to the rhs
      val = x0(localmesh->xend + 1, y);
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }
  }

  if (y_bndry == "dirichlet") {
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      // Should not go into corner cells, they are treated specially below
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->ystart - 1);

      // For the boundary value of the initial guess, use the value that would be set by
      // applying the boundary condition to the initial guess
      PetscScalar val =
          2. * x0(it.ind, localmesh->ystart - 1) - x0(it.ind, localmesh->ystart);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      // Pass the value from boundary cell of x0 as the boundary condition to the rhs
      val = x0(it.ind, localmesh->ystart - 1);
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      // Should not go into corner cells, they are treated specially below
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->yend + 1);

      // For the boundary value of the initial guess, use the value that would be set by
      // applying the boundary condition to the initial guess
      PetscScalar val =
          2. * x0(it.ind, localmesh->yend + 1) - x0(it.ind, localmesh->yend);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      // Pass the value from boundary cell of x0 as the boundary condition to the rhs
      val = x0(it.ind, localmesh->yend + 1);
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }

    // Use free_o3 for the corner boundary cells
    if (localmesh->hasBndryLowerY()) {
      if (localmesh->firstX()) {
        int ind = globalIndex(localmesh->xstart - 1, localmesh->ystart - 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = 3. * x0(localmesh->xstart - 1, localmesh->ystart)
                          - 3. * x0(localmesh->xstart - 1, localmesh->ystart + 1)
                          + x0(localmesh->xstart - 1, localmesh->ystart + 2);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
      if (localmesh->lastX()) {
        int ind = globalIndex(localmesh->xend + 1, localmesh->ystart - 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = 3. * x0(localmesh->xend + 1, localmesh->ystart)
                          - 3. * x0(localmesh->xend + 1, localmesh->ystart + 1)
                          + x0(localmesh->xend + 1, localmesh->ystart + 2);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    }
    if (localmesh->hasBndryUpperY()) {
      if (localmesh->firstX()) {
        int ind = globalIndex(localmesh->xstart - 1, localmesh->yend + 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = 3. * x0(localmesh->xstart - 1, localmesh->yend)
                          - 3. * x0(localmesh->xstart - 1, localmesh->yend - 1)
                          + x0(localmesh->xstart - 1, localmesh->yend - 2);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
      if (localmesh->lastX()) {
        int ind = globalIndex(localmesh->xend + 1, localmesh->yend + 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = 3. * x0(localmesh->xend + 1, localmesh->yend)
                          - 3. * x0(localmesh->xend + 1, localmesh->yend - 1)
                          + x0(localmesh->xend + 1, localmesh->yend - 2);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    }
  } else if (y_bndry == "neumann") {
    // Y boundaries Neumann
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      // Should not go into corner cells, they are treated specially below
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->ystart - 1);

      // For the boundary value of the initial guess, use the value that would be set by
      // applying the boundary condition to the initial guess
      PetscScalar val = x0(it.ind, localmesh->ystart);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      // Set the value for the rhs of the boundary condition to zero
      val = 0.0;
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }
    if (localmesh->hasBndryLowerY()) {
      if (localmesh->firstX()) {
        int ind = globalIndex(localmesh->xstart - 1, localmesh->ystart - 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = x0(localmesh->xstart - 1, localmesh->ystart);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
      if (localmesh->lastX()) {
        int ind = globalIndex(localmesh->xend + 1, localmesh->ystart - 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = x0(localmesh->xend + 1, localmesh->ystart);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      // Should not go into corner cells, they are treated specially below
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->yend + 1);

      // For the boundary value of the initial guess, use the value that would be set by
      // applying the boundary condition to the initial guess
      PetscScalar val = x0(it.ind, localmesh->yend);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      // Set the value for the rhs of the boundary condition to zero
      val = 0.0;
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }
    if (localmesh->hasBndryUpperY()) {
      if (localmesh->firstX()) {
        int ind = globalIndex(localmesh->xstart - 1, localmesh->yend + 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = x0(localmesh->xstart - 1, localmesh->yend);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
      if (localmesh->lastX()) {
        int ind = globalIndex(localmesh->xend + 1, localmesh->yend + 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = x0(localmesh->xend + 1, localmesh->yend);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    }
  } else if (y_bndry == "free_o3") {
    // Y boundaries free_o3
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      // Should not go into corner cells, they are treated specially below
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->ystart - 1);

      // For the boundary value of the initial guess, use the value that would be set by
      // applying the boundary condition to the initial guess
      PetscScalar val = 3. * x0(it.ind, localmesh->ystart)
                        - 3. * x0(it.ind, localmesh->ystart + 1)
                        + x0(it.ind, localmesh->ystart + 2);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      // Set the value for the rhs of the boundary condition to zero
      val = 0.0;
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }
    if (localmesh->hasBndryLowerY()) {
      if (localmesh->firstX()) {
        int ind = globalIndex(localmesh->xstart - 1, localmesh->ystart - 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = 3. * x0(localmesh->xstart - 1, localmesh->ystart)
                          - 3. * x0(localmesh->xstart - 1, localmesh->ystart + 1)
                          + x0(localmesh->xstart - 1, localmesh->ystart + 2);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
      if (localmesh->lastX()) {
        int ind = globalIndex(localmesh->xend + 1, localmesh->ystart - 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = 3. * x0(localmesh->xend + 1, localmesh->ystart)
                          - 3. * x0(localmesh->xend + 1, localmesh->ystart + 1)
                          + x0(localmesh->xend + 1, localmesh->ystart + 2);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      // Should not go into corner cells, they are treated specially below
      if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
        continue;
      }
      int ind = globalIndex(it.ind, localmesh->yend + 1);

      // For the boundary value of the initial guess, use the value that would be set by
      // applying the boundary condition to the initial guess
      PetscScalar val = 3. * x0(it.ind, localmesh->yend)
                        - 3. * x0(it.ind, localmesh->yend - 1)
                        + x0(it.ind, localmesh->yend - 2);
      VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

      // Set the value for the rhs of the boundary condition to zero
      val = 0.0;
      VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
    }
    if (localmesh->hasBndryUpperY()) {
      if (localmesh->firstX()) {
        int ind = globalIndex(localmesh->xstart - 1, localmesh->yend + 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = 3. * x0(localmesh->xstart - 1, localmesh->yend)
                          - 3. * x0(localmesh->xstart - 1, localmesh->yend - 1)
                          + x0(localmesh->xstart - 1, localmesh->yend - 2);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
      if (localmesh->lastX()) {
        int ind = globalIndex(localmesh->xend + 1, localmesh->yend + 1);

        // For the boundary value of the initial guess, use the value that would be set by
        // applying the boundary condition to the initial guess
        PetscScalar val = 3. * x0(localmesh->xend + 1, localmesh->yend)
                          - 3. * x0(localmesh->xend + 1, localmesh->yend - 1)
                          + x0(localmesh->xend + 1, localmesh->yend - 2);
        VecSetValues(xs, 1, &ind, &val, INSERT_VALUES);

        // Set the value for the rhs of the boundary condition to zero
        val = 0.0;
        VecSetValues(bs, 1, &ind, &val, INSERT_VALUES);
      }
    }
  } else {
    throw BoutException("Unsupported option for y_bndry");
  }
}

/*! Preconditioner
 * NOTE: For generality, this routine does use globalIndex() in the inner loop, although
 * this may be slightly less efficient than incrementing an integer for the global index,
 * the finite-volume and finite-difference implementations have slightly different
 * indexing patterns, so incrementing an integer would be tricky.
 */
int LaplaceXY::precon(Vec input, Vec result) {
  
  for(auto itdwn=localmesh->iterateBndryLowerY(); !itdwn.isDone(); itdwn++) {
    // Should not go into corner cells, LaplaceXY stencil does not include them
    if (itdwn.ind < localmesh->xstart or itdwn.ind > localmesh->xend) {
      continue;
    }
    const int ind = globalIndex(itdwn.ind, localmesh->ystart-1);
    PetscScalar val;
    VecGetValues(input, 1, &ind, &val );
    VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
  }
  for(auto itup=localmesh->iterateBndryUpperY(); !itup.isDone(); itup++) {
    // Should not go into corner cells, LaplaceXY stencil does not include them
    if (itup.ind < localmesh->xstart or itup.ind > localmesh->xend) {
      continue;
    }
    const int ind = globalIndex(itup.ind, localmesh->yend+1);
    PetscScalar val;
    VecGetValues(input, 1, &ind, &val );
    VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
  }
    
  // Load vector x into bvals array
  for(int x=xstart; x<=xend; x++) {
    for(int y=localmesh->ystart; y<=localmesh->yend; y++) {
      const int ind = globalIndex(x, y);
      PetscScalar val;
      VecGetValues(input, 1, &ind, &val );
      bvals(y - localmesh->ystart, x - xstart) = val;
    }
  }
  
  // Solve tridiagonal systems using CR solver
  cr->solve(bvals, xvals);

  // Save result xvals into y array
  for(int x=xstart; x<=xend; x++) {
    for(int y=localmesh->ystart; y<=localmesh->yend; y++) {
      const int ind = globalIndex(x, y);
      PetscScalar val = xvals(y - localmesh->ystart, x - xstart);
      VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
    }
  }
  VecAssemblyBegin(result);
  VecAssemblyEnd(result);
  return 0;
}

///////////////////////////////////////////////////////////////

int LaplaceXY::localSize() {
  
  // Bulk of points
  const int nx = localmesh->xend - localmesh->xstart + 1;
  const int ny = localmesh->yend - localmesh->ystart + 1;
  
  int n = nx * ny;  
  
  // X boundaries
  if(localmesh->firstX())
    n += ny;
  if(localmesh->lastX())
    n += ny;
  
  // Y boundaries
  for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    // Should not go into corner cells, LaplaceXY stencil does not include them
    if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
      continue;
    }
    n++;
  }
  if ((not finite_volume) and localmesh->hasBndryLowerY()) {
    if (localmesh->firstX()) {
      n++;
    }
    if (localmesh->lastX()) {
      n++;
    }
  }
  for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    // Should not go into corner cells, LaplaceXY stencil does not include them
    if (it.ind < localmesh->xstart or it.ind > localmesh->xend) {
      continue;
    }
    n++;
  }
  if ((not finite_volume) and localmesh->hasBndryUpperY()) {
    if (localmesh->firstX()) {
      n++;
    }
    if (localmesh->lastX()) {
      n++;
    }
  }
  
  return n;
}

int LaplaceXY::globalIndex(int x, int y) {
  if( (x < 0) || (x >= localmesh->LocalNx) ||
      (y < 0) || (y >= localmesh->LocalNy) )
    return -1; // Out of range
 
  // Get the index from a Field2D, round to integer
  return static_cast<int>(std::round(indexXY(x, y)));
}

void LaplaceXY::savePerformance(Datafile& output_file, Solver& solver,
                                std::string name) {
  // set flag so that performance monitoring values are calculated
  save_performance = true;

  // add values to be saved to the output
  if (name == "") {
    name = default_prefix;
  }
  output_file.addRepeat(output_average_iterations, name + "_average_iterations");

  // add monitor to reset counters/averages for new output timestep
  // monitor added to back of queue, so that values are reset after being saved
  solver.addMonitor(&monitor, Solver::BACK);
}
#endif // BOUT_HAS_PETSC
