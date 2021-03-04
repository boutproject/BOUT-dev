#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

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

LaplaceXY2::LaplaceXY2(Mesh *m, Options *opt, const CELL_LOC loc)
  : localmesh(m==nullptr ? bout::globals::mesh : m), f2dinit(localmesh), matrix(f2dinit), location(loc) {
  Timer timer("invert");

  if (opt == nullptr) {
    // If no options supplied, use default
    opt = &(Options::root()["laplacexy"]);
  }
  
  // Get MPI communicator
  MPI_Comm comm = BoutComm::get();
  
  //////////////////////////////////////////////////
  // Pre-allocate PETSc storage

  // Note: This is a significant performance optimisation
  
  //////////////////////////////////////////////////
  // Set up KSP
  
  // Declare KSP Context 
  KSPCreate( comm, &ksp ); 
  
  // Configure Linear Solver
  
  bool direct = (*opt)["direct"].doc("Use a direct LU solver").withDefault<bool>(false);
  
  if (direct) {
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);
#if PETSC_VERSION_GE(3,9,0)
    PCFactorSetMatSolverType(pc,"mumps");
#else
    PCFactorSetMatSolverPackage(pc,"mumps");
#endif
  } else {

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
                         .withDefault<bool>(true);

  ///////////////////////////////////////////////////
  // Set the default coefficients
  Field2D one(1., localmesh);
  Field2D zero(0., localmesh);
  one.setLocation(location);
  zero.setLocation(location);
  setCoefs(one, zero);
}

void LaplaceXY2::setCoefs(const Field2D &A, const Field2D &B) {
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

  for (auto& index_const : A.getRegion("RGN_NOBNDRY")) {
    // Note: This is needed for now because PetscMatrix::operator() takes non-const refs
    PetscMatrix<Field2D>::ind_type index = index_const;

    // Index offsets
    auto ind_xp = index.xp();
    auto ind_xm = index.xm();
    
    // XX component
    
    // Metrics on x+1/2 boundary
    BoutReal J = 0.5*(coords->J[index] + coords->J[ind_xp]);
    BoutReal g11 = 0.5*(coords->g11[index] + coords->g11[ind_xp]);
    BoutReal dx = 0.5*(coords->dx[index] + coords->dx[ind_xp]);
    BoutReal Acoef = 0.5*(A[index] + A[ind_xp]);
    
    BoutReal xp = Acoef * J * g11 / (coords->J[index] * dx * coords->dx[index]);
    
    // Metrics on x-1/2 boundary
    J = 0.5*(coords->J[index] + coords->J[ind_xm]);
    g11 = 0.5*(coords->g11[index] + coords->g11[ind_xm]);
    dx = 0.5*(coords->dx[index] + coords->dx[ind_xm]);
    Acoef = 0.5*(A[index] + A[ind_xm]);
      
    BoutReal xm = Acoef * J * g11 / (coords->J[index] * dx * coords->dx[index]);
    
    BoutReal c = B[index] - xp - xm; // Central coefficient

    matrix(index, ind_xp) = xp;
    matrix(index, ind_xm) = xm;
    
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
        
      BoutReal yp = -Acoef * J * g23 * g_23 / (g_22 * coords->J[index] * dy * coords->dy[index]);
      c -= yp;
      matrix(index, ind_yp) = yp;
      
      // Metrics at y-1/2
      J    = 0.5*(coords->J[index]    + coords->J[ind_ym]);
      g_22 = 0.5*(coords->g_22[index] + coords->g_22[ind_ym]);
      g23  = 0.5*(coords->g23[index]  + coords->g23[ind_ym]);
      g_23 = 0.5*(coords->g_23[index] + coords->g_23[ind_ym]);
      dy   = 0.5*(coords->dy[index]   + coords->dy[ind_ym]);
      Acoef = 0.5*(A[ind_ym] + A[index]);
      
      BoutReal ym = -Acoef * J * g23 * g_23 / (g_22 * coords->J[index] * dy * coords->dy[index]);
      c -= ym;
      matrix(index, ind_ym) = ym;
    }
    // Note: The central coefficient is done last because this may be modified
    // if y derivs are/are not included.
    matrix(index, index) = c;
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
      }
      
    } else {
      
      // Neumann on inner X boundary
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart, y);
        auto ind_xm = index.xm();

        matrix(ind_xm, index) =  1.0;
        matrix(ind_xm, ind_xm) = -1.0;
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
}

LaplaceXY2::~LaplaceXY2() {
  PetscBool is_finalised;
  PetscFinalized(&is_finalised);

  if (!is_finalised) {
    // PetscFinalize may already have destroyed this object
    KSPDestroy(&ksp);
  }
}

const Field2D LaplaceXY2::solve(const Field2D &rhs, const Field2D &x0) {
  Timer timer("invert");
  
  ASSERT1(rhs.getMesh() == localmesh);
  ASSERT1(x0.getMesh() == localmesh);
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  // Load initial guess x0 into xs and rhs into bs

  PetscVector<Field2D> xs(x0), bs(rhs);

  if (localmesh->firstX()) {
    if (x_inner_dirichlet) {
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart - 1, y);

        xs(index) = x0[index];
        bs(index) = 0.5 * (x0[index] + x0[index.xp()]);
      }
    } else {
      // Inner X boundary (Neumann)
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart - 1, y);

        xs(index) = x0[index];
        bs(index) = 0.0; // x0[index] - x0[index.xp()];
      }
    }
  }

  // Outer X boundary (Dirichlet)
  if (localmesh->lastX()) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      auto index = index2d(localmesh, localmesh->xend + 1, y);
      
      xs(index) = x0[index];
      bs(index) = 0.5 * (x0[index.xm()] + x0[index]);
    }
  }

  if (y_bndry_dirichlet) {
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart - 1);
      
      xs(index) = x0[index];
      bs(index) = 0.5 * (x0[index] + x0[index.yp()]);
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->yend + 1);

      xs(index) = x0[index];
      bs(index) = 0.5 * (x0[index] + x0[index.xm()]);
    }
  } else {
    // Y boundaries Neumann
    for(RangeIterator it=localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart - 1);
      
      xs(index) = x0[index];
      bs(index) = 0.0;
    }
  
    for(RangeIterator it=localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->yend + 1);
      
      xs(index) = x0[index];
      bs(index) = 0.0;
    }
  }

  // Assemble RHS Vector
  bs.assemble();

  // Assemble Trial Solution Vector
  xs.assemble();
  
  // Solve the system
  KSPSolve( ksp, *bs.get(), *xs.get() );
  
  KSPConvergedReason reason;
  KSPGetConvergedReason( ksp, &reason );

  if (reason <= 0) {
    throw BoutException("LaplaceXY2 failed to converge. Reason {:d}", reason);
  }

  // Convert result into a Field2D
  auto result = xs.toField();

  // Set boundary cells past the first one
  ////////////////////////////////////////

  // Inner X boundary
  if (localmesh->firstX()) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      for (int x = localmesh->xstart - 2; x >= 0; x--)
        result(x, y) = result(localmesh->xstart-1, y);
    }
  }

  // Outer X boundary
  if (localmesh->lastX()) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      for (int x = localmesh->xend + 2; x < localmesh->LocalNx; x++)
        result(x, y) = result(localmesh->xend + 1, y);
    }
  }

  // Lower Y boundary
  for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    for (int y = localmesh->ystart - 2; y >= 0; y--)
      result(it.ind, y) = result(it.ind, localmesh->ystart - 1);
  }

  // Upper Y boundary
  for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    for (int y = localmesh->yend + 2; y < localmesh->LocalNy; y++)
      result(it.ind, y) = result(it.ind, localmesh->yend + 1);
  }

  return result;
}


#endif // BOUT_HAS_PETSC
