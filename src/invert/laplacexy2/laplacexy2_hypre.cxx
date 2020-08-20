
#ifdef BOUT_HAS_HYPRE

#include <bout/invert/laplacexy2_hypre.hxx>

#include <bout/assert.hxx>

#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <globals.hxx>
#include <utils.hxx>

#include <output.hxx>

#include <cmath>

Ind2D index2d(Mesh* mesh, int x, int y) {
  int ny = mesh->LocalNy;
  return Ind2D(x * ny + y, ny, 1);
}

LaplaceXY2Hypre::LaplaceXY2Hypre(Mesh* m, Options* opt, const CELL_LOC loc)
    : localmesh(m == nullptr ? bout::globals::mesh : m), f2dinit(localmesh),
      matrix(f2dinit), location(loc) {
  Timer timer("invert");

  if (opt == nullptr) {
    // If no options supplied, use default
    opt = &(Options::root()["laplacexy"]);
  }

  //////////////////////////////////////////////////
  // Pre-allocate Hypre storage

  // Note: This is a significant performance optimisation

  //////////////////////////////////////////////////
  // Set up KSP

  // Declare KSP Context
  HYPRE_BoomerAMGCreate(&solver);

  // Configure Linear Solver

  // Defaults taken from one of the HYPRE examples (ex5.c)
  const BoutReal tol = (*opt)["tol"].doc("Tolerance").withDefault(1e-7);
  const int relax_type = (*opt)["relax_type"].doc("What smoother to use").withDefault(3);
  const bool cf_relaxation =
      (*opt)["cf_relaxtion"].doc("Use CF-relaxation").withDefault(true);
  const int num_sweeps = (*opt)["num_sweeps"].doc("Number of sweeps").withDefault(1);
  const int max_levels =
      (*opt)["max_levels"].doc("Maximum number of multigrid levels").withDefault(20);

#if CHECKLEVEL >= 1
  HYPRE_BoomerAMGSetPrintLevel(solver, 3);
#endif

  // Falgout coarsening with modified classical interpolaiton
  HYPRE_BoomerAMGSetOldDefault(solver);
  // G-S/Jacobi hybrid relaxation
  HYPRE_BoomerAMGSetRelaxType(solver, relax_type);
  // uses C/F relaxation
  HYPRE_BoomerAMGSetRelaxOrder(solver, cf_relaxation);
  // Sweeps on each level
  HYPRE_BoomerAMGSetNumSweeps(solver, num_sweeps);
  // maximum number of levels
  HYPRE_BoomerAMGSetMaxLevels(solver, max_levels);
  // convergence tolerance
  HYPRE_BoomerAMGSetTol(solver, tol);

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

void LaplaceXY2Hypre::setCoefs(const Field2D& A, const Field2D& B) {
  Timer timer("invert");

  ASSERT1(A.getMesh() == localmesh);
  ASSERT1(B.getMesh() == localmesh);
  ASSERT1(A.getLocation() == location);
  ASSERT1(B.getLocation() == location);

  Coordinates* coords = localmesh->getCoordinates(location);

  //////////////////////////////////////////////////
  // Set Matrix elements
  //
  // (1/J) d/dx ( J * g11 d/dx ) + (1/J) d/dy ( J * g22 d/dy )

  for (auto& index : A.getRegion("RGN_NOBNDRY")) {
    // Index offsets
    auto ind_xp = index.xp();
    auto ind_xm = index.xm();

    // XX component

    // Metrics on x+1/2 boundary
    BoutReal J = 0.5 * (coords->J[index] + coords->J[ind_xp]);
    BoutReal g11 = 0.5 * (coords->g11[index] + coords->g11[ind_xp]);
    BoutReal dx = 0.5 * (coords->dx[index] + coords->dx[ind_xp]);
    BoutReal Acoef = 0.5 * (A[index] + A[ind_xp]);

    BoutReal xp = Acoef * J * g11 / (coords->J[index] * dx * coords->dx[index]);

    // Metrics on x-1/2 boundary
    J = 0.5 * (coords->J[index] + coords->J[ind_xm]);
    g11 = 0.5 * (coords->g11[index] + coords->g11[ind_xm]);
    dx = 0.5 * (coords->dx[index] + coords->dx[ind_xm]);
    Acoef = 0.5 * (A[index] + A[ind_xm]);

    BoutReal xm = Acoef * J * g11 / (coords->J[index] * dx * coords->dx[index]);

    BoutReal c = B[index] - xp - xm; // Central coefficient

    matrix(index, ind_xp) = xp;
    matrix(index, ind_xm) = xm;

    if (include_y_derivs) {
      auto ind_yp = index.yp();
      auto ind_ym = index.ym();

      // YY component
      // Metrics at y+1/2
      J = 0.5 * (coords->J[index] + coords->J[ind_yp]);
      BoutReal g_22 = 0.5 * (coords->g_22[index] + coords->g_22[ind_yp]);
      BoutReal g23 = 0.5 * (coords->g23[index] + coords->g23[ind_yp]);
      BoutReal g_23 = 0.5 * (coords->g_23[index] + coords->g_23[ind_yp]);
      BoutReal dy = 0.5 * (coords->dy[index] + coords->dy[ind_yp]);
      Acoef = 0.5 * (A[ind_yp] + A[index]);

      BoutReal yp =
          -Acoef * J * g23 * g_23 / (g_22 * coords->J[index] * dy * coords->dy[index]);
      c -= yp;
      matrix(index, ind_yp) = yp;

      // Metrics at y-1/2
      J = 0.5 * (coords->J[index] + coords->J[ind_ym]);
      g_22 = 0.5 * (coords->g_22[index] + coords->g_22[ind_ym]);
      g23 = 0.5 * (coords->g23[index] + coords->g23[ind_ym]);
      g_23 = 0.5 * (coords->g_23[index] + coords->g_23[ind_ym]);
      dy = 0.5 * (coords->dy[index] + coords->dy[ind_ym]);
      Acoef = 0.5 * (A[ind_ym] + A[index]);

      BoutReal ym =
          -Acoef * J * g23 * g_23 / (g_22 * coords->J[index] * dy * coords->dy[index]);
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

        matrix(ind_xm, index) = 1.0;
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
  HYPRE_BoomerAMGSetup(solver, matrix.getParallel(), nullptr, nullptr);
}

LaplaceXY2Hypre::~LaplaceXY2Hypre() {
  if (solver != nullptr) {
    HYPRE_BoomerAMGDestroy(solver);
  }
}

const Field2D LaplaceXY2Hypre::solve(const Field2D& rhs, const Field2D& x0) {
  Timer timer("invert");

  ASSERT1(rhs.getMesh() == localmesh);
  ASSERT1(x0.getMesh() == localmesh);
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  // Load initial guess x0 into xs and rhs into bs

  bout::HypreVector<Field2D> xs(x0), bs(rhs);

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
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart - 1);

      xs(index) = x0[index];
      bs(index) = 0.0;
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
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
  HYPRE_BoomerAMGSolve(solver, matrix.getParallel(), bs.getParallel(), xs.getParallel());

  // Convert result into a Field2D
  return xs.toField();
}

#endif // BOUT_HAS_HYPRE
