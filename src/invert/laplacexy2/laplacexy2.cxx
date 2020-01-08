
#ifdef BOUT_HAS_PETSC

#include <petscksp.h>

#include <bout/invert/laplacexy2.hxx>

#include <bout/assert.hxx>

#include <boutcomm.hxx>
#include <globals.hxx>
#include <utils.hxx>
#include <bout/sys/timer.hxx>

#include <output.hxx>

#include <cmath>

Ind2D index2d(Mesh* mesh, int x, int y) {
  int ny = mesh->LocalNy;
  return Ind2D(x * ny + y, ny, 1);
}

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
  : matrix(Field2D(m)), localmesh(m==nullptr ? bout::globals::mesh : m), location(loc) {
  Timer timer("invert");

  if (opt == nullptr) {
    // If no options supplied, use default
    opt = &(Options::root()["laplacexy"]);
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

  // Note: This is a significant performance optimisation
  
  //////////////////////////////////////////////////
  // Set up KSP
  
  // Declare KSP Context 
  KSPCreate( comm, &ksp ); 
  
  // Configure Linear Solver
  
  bool direct = (*opt)["direct"].doc("Use a direct LU solver").withDefault(false);
  
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
  
  KSPSetFromOptions( ksp );

  ///////////////////////////////////////////////////
  // Decide boundary condititions
  if (localmesh->periodicY(localmesh->xstart)) {
    // Periodic in Y, so in the core
    opt->get("core_bndry_dirichlet", x_inner_dirichlet, false);
  } else {
    // Non-periodic, so in the PF region
    opt->get("pf_bndry_dirichlet", x_inner_dirichlet, true);
  }
  opt->get("y_bndry_dirichlet", y_bndry_dirichlet, false);

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

void LaplaceXY::setCoefs(const Field2D &A, const Field2D &B) {
  Timer timer("invert");

  ASSERT1(A.getMesh() == localmesh);
  ASSERT1(B.getMesh() == localmesh);
  ASSERT1(A.getLocation() == location);
  ASSERT1(B.getLocation() == location);

  Coordinates *coords = localmesh->getCoordinates(location);
  
  //////////////////////////////////////////////////
  // Set Matrix elements
  // 
  // (1/J) d/dx ( J * g11 d/dx ) + (1/J) d/dy ( J * g22 d/dy )

  for (auto& index : A.getRegion("RGN_NOBNDRY")) {
    // stencil entries
    PetscScalar c, xm, xp, ym, yp;

    auto ind_xp = index.xp();
    auto ind_xm = index.xm();
    
    // XX component
    
    // Metrics on x+1/2 boundary
    BoutReal J = 0.5*(coords->J[index] + coords->J[ind_xp]);
    BoutReal g11 = 0.5*(coords->g11[index] + coords->g11[ind_xp]);
    BoutReal dx = 0.5*(coords->dx[index] + coords->dx[ind_xp]);
    BoutReal Acoef = 0.5*(A[index] + A[ind_xp]);
    
    BoutReal val = Acoef * J * g11 / (coords->J[index] * dx * coords->dx[index]);
    xp = val;
    c  = -val;
    
    // Metrics on x-1/2 boundary
    J = 0.5*(coords->J[index] + coords->J[ind_xm]);
    g11 = 0.5*(coords->g11[index] + coords->g11[ind_xm]);
    dx = 0.5*(coords->dx[index] + coords->dx[ind_xm]);
    Acoef = 0.5*(A[index] + A[ind_xm]);
      
    val = Acoef * J * g11 / (coords->J[index] * dx * coords->dx[index]);
    xm = val;
    c  -= val;
    
    c += B[index];
      
    // Put values into the preconditioner, X derivatives only
    int x = index.x(), y = index.y();
    acoef(y - localmesh->ystart, x - xstart) = xm;
    bcoef(y - localmesh->ystart, x - xstart) = c;
    ccoef(y - localmesh->ystart, x - xstart) = xp;

    if( include_y_derivs ) {
      auto ind_yp = index.yp();
      auto ind_ym = index.ym();
      
      // YY component
      // Metrics at y+1/2
      J = 0.5*(coords->J[index] + coords->J[ind_yp]);
      BoutReal g_22 = 0.5*(coords->g_22[index] + coords->g_22[ind_yp]);
      BoutReal g23  = 0.5*(coords->g23[index] + coords->g23[ind_yp]);
      BoutReal g_23 = 0.5*(coords->g_23[index] + coords->g_23[ind_yp]);
      BoutReal dy   = 0.5*(coords->dy[index] + coords->dy[ind_yp]);
      Acoef = 0.5*(A[ind_yp] + A[index]);
        
      val = -Acoef * J * g23 * g_23 / (g_22 * coords->J[index] * dy * coords->dy[index]);
      yp = val;
      c -= val;
        
      // Metrics at y-1/2
      J    = 0.5*(coords->J[index]    + coords->J[ind_ym]);
      g_22 = 0.5*(coords->g_22[index] + coords->g_22[ind_ym]);
      g23  = 0.5*(coords->g23[index]  + coords->g23[ind_ym]);
      g_23 = 0.5*(coords->g_23[index] + coords->g_23[ind_ym]);
      dy   = 0.5*(coords->dy[index]   + coords->dy[ind_ym]);
      Acoef = 0.5*(A[ind_ym] + A[index]);
      
      val = -Acoef * J * g23 * g_23 / (g_22 * coords->J[index] * dy * coords->dy[index]);
      ym = val;
      c -= val;
    }
      
    /////////////////////////////////////////////////
    // Now have a 5-point stencil for the Laplacian

    matrix(index, index) = c;
    matrix(index, ind_xp) = xp;
    matrix(index, ind_xm) = xm;
    
    if( include_y_derivs ) {
      // Y + 1
      matrix(index, ind_yp) = yp;
      
      // Y - 1
      matrix(index, ind_ym) = ym;
    }
  }
  
  // X boundaries
  if (localmesh->firstX()) {
    if (x_inner_dirichlet) {

      // Dirichlet on inner X boundary
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart, y);
        auto ind_xm = index.xm();
                            
        matrix(ind_xm, index) = 0.5;
        matrix(ind_xm, ind_xm) = 0.5;
        
        // Preconditioner
        bcoef(y - localmesh->ystart, 0) = 0.5;
        ccoef(y - localmesh->ystart, 0) = 0.5;
      }
      
    } else {
      
      // Neumann on inner X boundary
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart, y);
        auto ind_xm = index.xm();

        matrix(ind_xm, index) =  1.0;
        matrix(ind_xm, ind_xm) = -1.0;
        
        // Preconditioner
        bcoef(y - localmesh->ystart, 0) = 1.0;
        ccoef(y - localmesh->ystart, 0) = -1.0;
      }
    }
  }
  if (localmesh->lastX()) {
    // Dirichlet on outer X boundary

    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      auto index = index2d(localmesh, localmesh->xend, y);
      auto ind_xp = index.xp();

      matrix(ind_xp, ind_xp) = 0.5;
      matrix(ind_xp, index) = 0.5;
      
      // Preconditioner
      acoef(y - localmesh->ystart, localmesh->xend + 1 - xstart) = 0.5;
      bcoef(y - localmesh->ystart, localmesh->xend + 1 - xstart) = 0.5;
    }
  }

  if (y_bndry_dirichlet) {
    // Dirichlet on Y boundaries
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart);
      auto ind_ym = index.ym();

      matrix(ind_ym, ind_ym) = 0.5;
      matrix(ind_ym, index) = 0.5;
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {

      auto index = index2d(localmesh, it.ind, localmesh->yend);
      auto ind_yp = index.yp();

      matrix(ind_yp, ind_yp) = 0.5;
      matrix(ind_yp, index) = 0.5;
    }
  } else {
    // Neumann on Y boundaries
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart);
      auto ind_ym = index.ym();

      matrix(ind_ym, ind_ym) = -1.0;
      matrix(ind_ym, index) = 1.0;
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->yend);
      auto ind_yp = index.yp();

      matrix(ind_yp, ind_yp) = 1.0;
      matrix(ind_yp, index) = -1.0;
    }
  }
  
  // Assemble Matrix
  matrix.assemble();
  
  // Set the operator
#if PETSC_VERSION_GE(3,5,0)
  KSPSetOperators(ksp, *matrix.get(), *matrix.get());
#else
  KSPSetOperators(ksp, *matrix.get(), *matrix.get(), DIFFERENT_NONZERO_PATTERN);
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
}

const Field2D LaplaceXY::solve(const Field2D &rhs, const Field2D &x0) {
  Timer timer("invert");
  
  ASSERT1(rhs.getMesh() == localmesh);
  ASSERT1(x0.getMesh() == localmesh);
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  // Load initial guess x0 into xs and rhs into bs
  
  for(int x=localmesh->xstart;x<= localmesh->xend;x++) {
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      int ind = globalIndex(x,y);
      
      PetscScalar val = x0(x,y);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      
      val = rhs(x,y);
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  }

  if(localmesh->firstX()) {
    if(x_inner_dirichlet) {
      for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
        int ind = globalIndex(localmesh->xstart-1,y);
      
        PetscScalar val = x0(localmesh->xstart-1,y);
        VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
        
        val = 0.5*(x0(localmesh->xstart-1,y) + x0(localmesh->xstart,y));
        VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
      }
    }else {
      // Inner X boundary (Neumann)
      for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
        int ind = globalIndex(localmesh->xstart-1,y);
        
        PetscScalar val = x0(localmesh->xstart-1,y);
        VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
        
        val = 0.0; //x0(localmesh->xstart-1,y) - x0(localmesh->xstart,y);
        VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
      }
    }
  }
  
  // Outer X boundary (Dirichlet)
  if(localmesh->lastX()) {
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      int ind = globalIndex(localmesh->xend+1,y);
      
      PetscScalar val = x0(localmesh->xend+1,y);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      
      val = 0.5*(x0(localmesh->xend,y) + x0(localmesh->xend+1,y));
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  }

  if(y_bndry_dirichlet) {
    for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, localmesh->ystart-1);
    
      PetscScalar val = x0(it.ind,localmesh->ystart-1);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
    
      val = 0.5*(x0(it.ind, localmesh->ystart-1) + x0(it.ind, localmesh->ystart));
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  
    for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, localmesh->yend+1);
    
      PetscScalar val = x0(it.ind,localmesh->yend+1);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
    
      val = 0.5*(x0(it.ind, localmesh->yend+1) + x0(it.ind, localmesh->yend));
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  } else {
    // Y boundaries Neumann
    for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, localmesh->ystart-1);
    
      PetscScalar val = x0(it.ind,localmesh->ystart-1);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      
      val = 0.0;
      VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
    }
  
    for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      int ind = globalIndex(it.ind, localmesh->yend+1);
      
      PetscScalar val = x0(it.ind,localmesh->yend+1);
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
    throw BoutException("LaplaceXY failed to converge. Reason {:d}", reason);
  }
  
  //////////////////////////
  // Copy data into result
  
  Field2D result;
  result.allocate();
  result.setLocation(location);
  
  for(int x=localmesh->xstart;x<= localmesh->xend;x++) {
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      int ind = globalIndex(x,y);
      
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      result(x,y) = val;
    }
  }
  
  // Inner X boundary
  if(localmesh->firstX()) {
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      int ind = globalIndex(localmesh->xstart-1,y);
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      for(int x=localmesh->xstart-1; x >= 0; x--)
        result(x,y) = val;
    }
  }

  // Outer X boundary
  if(localmesh->lastX()) {
    for(int y=localmesh->ystart;y<=localmesh->yend;y++) {
      int ind = globalIndex(localmesh->xend+1,y);
      PetscScalar val;
      VecGetValues(xs, 1, &ind, &val );
      for(int x=localmesh->xend+1;x < localmesh->LocalNx;x++)
        result(x,y) = val;
    }
  }
  
  // Lower Y boundary
  for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, localmesh->ystart-1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val );
    for(int y=localmesh->ystart-1;y>=0;y--)
      result(it.ind, y) = val;
  }
  
  // Upper Y boundary
  for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    int ind = globalIndex(it.ind, localmesh->yend+1);
    PetscScalar val;
    VecGetValues(xs, 1, &ind, &val );
    for(int y=localmesh->yend+1;y<localmesh->LocalNy;y++)
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
  
  RangeIterator itdwn=localmesh->iterateBndryLowerY();
  if(!itdwn.isDone()) {
    ind = globalIndex(itdwn.ind, localmesh->ystart-1);
    
    for(; !itdwn.isDone(); itdwn++) {
      PetscScalar val;
      VecGetValues(input, 1, &ind, &val ); 
      VecSetValues(result, 1, &ind, &val, INSERT_VALUES );
      ind++;
    }
  }
  RangeIterator itup=localmesh->iterateBndryUpperY();
  if(!itup.isDone()) {
    if(ind == -1) {
      // No lower boundary
      ind = globalIndex(itup.ind, localmesh->yend+1);
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
    ind = globalIndex(xstart, localmesh->ystart);
  }
    
  int ind0 = ind;
  // Load vector x into bvals array
  for(int x=xstart;x<=xend;x++) {
    for(int y=localmesh->ystart; y<=localmesh->yend;y++) {
      PetscScalar val;
      VecGetValues(input, 1, &ind, &val );
      bvals(y - localmesh->ystart, x - xstart) = val;
      ind++;
    }
  }
  
  // Solve tridiagonal systems using CR solver
  cr->solve(bvals, xvals);

  // Save result xvals into y array
  ind = ind0;
  for(int x=xstart;x<=xend;x++) {
    for(int y=localmesh->ystart; y<=localmesh->yend;y++) {
      PetscScalar val = xvals(y - localmesh->ystart, x - xstart);
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
  int nx = localmesh->xend - localmesh->xstart + 1;
  int ny = localmesh->yend - localmesh->ystart + 1;
  
  int n = nx * ny;  
  
  // X boundaries
  if(localmesh->firstX())
    n += ny;
  if(localmesh->lastX())
    n += ny;
  
  // Y boundaries
  for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    n++;
  }
  for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    n++;
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
#endif // BOUT_HAS_PETSC
