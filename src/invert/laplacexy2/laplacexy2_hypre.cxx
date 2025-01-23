#include <bout/invert/laplacexy2_hypre.hxx>

#if BOUT_HAS_HYPRE

#include "bout/assert.hxx"
#include "bout/boutcomm.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"
#include "bout/sys/timer.hxx"

#include <cmath>

#if BOUT_HAS_CUDA && defined(__CUDACC__)
#define gpuErrchk(ans) \
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) {
      exit(code);
    }
  }
}
#endif

LaplaceXY2Hypre::LaplaceXY2Hypre(Mesh* m, Options* opt, const CELL_LOC loc)
    : localmesh(m == nullptr ? bout::globals::mesh : m),
      indexConverter(std::make_shared<GlobalIndexer<Field2D>>(
          localmesh, squareStencil<Field2D::ind_type>(localmesh))),
      M(indexConverter), x(indexConverter), b(indexConverter),
      linearSystem(*localmesh, (opt == nullptr) ? Options::root()["laplacexy"] : *opt),
      location(loc) {
  Timer timer("invert");

  if (opt == nullptr) {
    // If no options supplied, use default
    opt = &(Options::root()["laplacexy"]);
  }

  linearSystem.setMatrix(&M);
  linearSystem.setSolutionVector(&x);
  linearSystem.setRHSVector(&b);

  ///////////////////////////////////////////////////
  // Decide boundary condititions
  if (localmesh->periodicY(localmesh->xstart)) {
    // Periodic in Y, so in the core
    x_inner_dirichlet = (*opt)["core_bndry_dirichlet"].withDefault(false);
  } else {
    // Non-periodic, so in the PF region
    x_inner_dirichlet = (*opt)["pf_bndry_dirichlet"].withDefault(true);
  }
  y_bndry_dirichlet = (*opt)["y_bndry_dirichlet"].withDefault(false);

  ///////////////////////////////////////////////////
  // Including Y derivatives?

  include_y_derivs = (*opt)["include_y_derivs"]
                         .doc("Include Y derivatives in operator to invert?")
                         .withDefault<bool>(true);

  print_timing = (*opt)["print_timing"]
                     .doc("Print extra timing information for LaplaceXY2Hypre")
                     .withDefault(false);

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

  if (print_timing) {
    output_info.write("setting up matrix...\n");
  }
  auto start = std::chrono::system_clock::now();
  for (auto& index : indexConverter->getRegionNobndry()) {
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

    M(index, ind_xp) = xp;
    M(index, ind_xm) = xm;

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
      M(index, ind_yp) = yp;
      M(index, ind_ym) = ym;
    }
    // Note: The central coefficient is done last because this may be modified
    // if y derivs are/are not included.
    M(index, index) = c;
  }

  // X boundaries
  if (x_inner_dirichlet) {
    // Dirichlet on inner X boundary
    BOUT_FOR_SERIAL(i, indexConverter->getRegionInnerX()) {
      M(i, i.xp()) = 0.5;
      M(i, i) = 0.5;
    }

  } else {
    // Neumann on inner X boundary
    BOUT_FOR_SERIAL(i, indexConverter->getRegionInnerX()) {
      M(i, i.xp()) = 1.0;
      M(i, i) = -1.0;
    }
  }
  // Dirichlet on outer X boundary
  BOUT_FOR_SERIAL(i, indexConverter->getRegionOuterX()) {
    M(i, i) = 0.5;
    M(i, i.xm()) = 0.5;
  }

  if (y_bndry_dirichlet) {
    // Dirichlet on Y boundaries
    BOUT_FOR_SERIAL(i, indexConverter->getRegionLowerY()) {
      M(i, i) = 0.5;
      M(i, i.yp()) = 0.5;
    }

    BOUT_FOR_SERIAL(i, indexConverter->getRegionUpperY()) {
      M(i, i) = 0.5;
      M(i, i.ym()) = 0.5;
    }
  } else {
    // Neumann on Y boundaries
    BOUT_FOR_SERIAL(i, indexConverter->getRegionLowerY()) {
      M(i, i) = -1.0;
      M(i, i.yp()) = 1.0;
    }

    BOUT_FOR_SERIAL(i, indexConverter->getRegionUpperY()) {
      M(i, i) = 1.0;
      M(i, i.ym()) = -1.0;
    }
  }

  auto end = std::chrono::system_clock::now();
  if (print_timing) {
    const auto dur = end - start;
    output_info.write("*****Matrix set time:  {}\n", dur.count());
  }

  start = std::chrono::system_clock::now();
  M.assemble();

  if (print_timing) {
    end = std::chrono::system_clock::now();
    const auto dur = end - start;
    output_info.write("*****Matrix asm time:  {}\n", dur.count());
  }

  start = std::chrono::system_clock::now();
  linearSystem.setupAMG(&M);

  if (print_timing) {
    end = std::chrono::system_clock::now();
    const auto dur = end - start;
    output_info.write("*****Matrix prec time:  {}\n", dur.count());
  }
}

Field2D LaplaceXY2Hypre::solve(Field2D& rhs, Field2D& x0) {
  Timer timer("invert");

  ASSERT1(rhs.getMesh() == localmesh);
  ASSERT1(x0.getMesh() == localmesh);
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  auto start = std::chrono::system_clock::now();

  x.importValuesFromField(x0);
  b.importValuesFromField(rhs);

  // Set boundary values
  //////////////////////

  if (x_inner_dirichlet) {
    BOUT_FOR_SERIAL(i, indexConverter->getRegionInnerX()) {
      b(i) = 0.5 * (x0[i] + x0[i.xp()]);
    }
  } else {
    // Inner X boundary (Neumann)
    BOUT_FOR_SERIAL(i, indexConverter->getRegionInnerX()) { b(i) = 0.0; }
  }

  // Outer X boundary (Dirichlet)
  BOUT_FOR_SERIAL(i, indexConverter->getRegionOuterX()) {
    b(i) = 0.5 * (x0[i.xm()] + x0[i]);
  }

  if (y_bndry_dirichlet) {
    BOUT_FOR_SERIAL(i, indexConverter->getRegionLowerY()) {
      b(i) = 0.5 * (x0[i] + x0[i.yp()]);
    }

    BOUT_FOR_SERIAL(i, indexConverter->getRegionUpperY()) {
      b(i) = 0.5 * (x0[i] + x0[i.ym()]);
    }
  } else {
    // Y boundaries Neumann
    BOUT_FOR_SERIAL(i, indexConverter->getRegionLowerY()) { b(i) = 0.0; }

    BOUT_FOR_SERIAL(i, indexConverter->getRegionUpperY()) { b(i) = 0.0; }
  }

  x.assemble();
  b.assemble();

  auto form_vec = std::chrono::system_clock::now();

  if (print_timing) {
    const auto formvec_dur = form_vec - start;
    output_info.write("*****Form Vectors time:  {}\n", formvec_dur.count());
  }

  start = std::chrono::system_clock::now();
  linearSystem.solve();

  auto slv = std::chrono::system_clock::now();

  if (print_timing) {
    const auto slv_dur = slv - start;
    output_info.write("*****BoomerAMG solve time:  {}\n", slv_dur.count());
  }

  // Convert result into a Field2D
  start = std::chrono::system_clock::now();
  Field2D sol = x.toField();

  auto formfield = std::chrono::system_clock::now();

  if (print_timing) {
    const auto formfield_dur = formfield - start;
    output_info.write("*****Form field time:  {}\n", formfield_dur.count());
  }

  // Set boundary cells past the first one
  ////////////////////////////////////////
  // Inner X boundary
  if (localmesh->firstX()) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      for (int x = localmesh->xstart - 2; x >= 0; x--) {
        sol(x, y) = sol(localmesh->xstart - 1, y);
      }
    }
  }

  // Outer X boundary
  if (localmesh->lastX()) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      for (int x = localmesh->xend + 2; x < localmesh->LocalNx; x++) {
        sol(x, y) = sol(localmesh->xend + 1, y);
      }
    }
  }

  // Lower Y boundary
  for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    for (int y = localmesh->ystart - 2; y >= 0; y--) {
      sol(it.ind, y) = sol(it.ind, localmesh->ystart - 1);
    }
  }

  // Upper Y boundary
  for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    for (int y = localmesh->yend + 2; y < localmesh->LocalNy; y++) {
      sol(it.ind, y) = sol(it.ind, localmesh->yend + 1);
    }
  }
  return sol;
}

#endif // BOUT_HAS_HYPRE
