
#ifdef BOUT_HAS_HYPRE

#include <bout/invert/laplacexy2_hypre.hxx>

#include <bout/assert.hxx>

#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <globals.hxx>
#include <utils.hxx>

#include <output.hxx>

#include <cmath>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

Ind2D index2d(Mesh* mesh, int x, int y) {
  int ny = mesh->LocalNy;
  return Ind2D(x * ny + y, ny, 1);
}

LaplaceXY2Hypre::LaplaceXY2Hypre(Mesh* m, Options* opt, const CELL_LOC loc)
    : localmesh(m == nullptr ? bout::globals::mesh : m), f2dinit(localmesh),
      location(loc) {
  Timer timer("invert");

  if (opt == nullptr) {
    // If no options supplied, use default
    opt = &(Options::root()["laplacexy"]);
  }

  indexConverter = std::make_shared<GlobalIndexer<Field2D>>(localmesh, squareStencil<Field2D::ind_type>(localmesh));

  linearSystem = new bout::HypreSystem<Field2D>(*localmesh);
  M = new bout::HypreMatrix<Field2D>(indexConverter);
  x = new bout::HypreVector<Field2D>(indexConverter);
  b = new bout::HypreVector<Field2D>(indexConverter);
  linearSystem->setMatrix(M);
  linearSystem->setSolutionVector(x);
  linearSystem->setRHSVector(b);

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

  //const auto& region = f2dinit.getRegion("RGN_NOBNDRY");
  const auto &region = indexConverter->getRegionAll();

  Coordinates* coords = localmesh->getCoordinates(location);

  //////////////////////////////////////////////////
  // Set Matrix elements
  //



  // (1/J) d/dx ( J * g11 d/dx ) + (1/J) d/dy ( J * g22 d/dy )
  std::cout << "setting up matrix..." << std::endl;
  auto start = std::chrono::system_clock::now();  //AARON
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

    (*M)(index, ind_xp) = xp;
    (*M)(index, ind_xm) = xm;

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
      (*M)(index, ind_yp) = yp;
      (*M)(index, ind_ym) = ym;  
    }
    // Note: The central coefficient is done last because this may be modified
    // if y derivs are/are not included.
    (*M)(index, index) = c;
  }

  // X boundaries
  if (localmesh->firstX()) {
    if (x_inner_dirichlet) {

      // Dirichlet on inner X boundary
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart, y);
        auto ind_xm = index.xm();
        (*M)(ind_xm, index) = 0.5;
        (*M)(ind_xm, ind_xm) = 0.5;
      }

    } else {

      // Neumann on inner X boundary
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart, y);
        auto ind_xm = index.xm();
        HYPRE_BigInt I = indexConverter->getGlobal(index);
        HYPRE_BigInt XM = indexConverter->getGlobal(ind_xm);       
        (*M)(ind_xm, index) = 1.0;
        (*M)(ind_xm, ind_xm) = -1.0;
      }
    }
  }
  if (localmesh->lastX()) {
    // Dirichlet on outer X boundary

    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      auto index = index2d(localmesh, localmesh->xend, y);      
      auto ind_xp = index.xp();
      (*M)(ind_xp, ind_xp) = 0.5;
      (*M)(ind_xp, index) = 0.5;   
    }
  }

  if (y_bndry_dirichlet) {
    // Dirichlet on Y boundaries 
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart);
      auto ind_ym = index.ym();
      (*M)(ind_ym, ind_ym) = 0.5;
      (*M)(ind_ym, index) = 0.5;       
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->yend);
      auto ind_yp = index.yp();
      (*M)(ind_yp, ind_yp) = 0.5;
      (*M)(ind_yp, index) = 0.5;
    }
  } else {
    // Neumann on Y boundaries   
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart);
      auto ind_ym = index.ym();
      (*M)(ind_ym, ind_ym) = -1.0;
      (*M)(ind_ym, index) = 1.0;    
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->yend);
      auto ind_yp = index.yp();
      (*M)(ind_yp, ind_yp) = 1.0;
      (*M)(ind_yp, index) = -1.0;   
    }
  }
  auto end = std::chrono::system_clock::now();  //AARON
  std::chrono::duration<double> dur = end-start;  //AARON
  std::cout << "*****Matrix set time:  " << dur.count() << std::endl;    

  start = std::chrono::system_clock::now();
  // Assemble Matrix
  M->assemble();

  end = std::chrono::system_clock::now();  //AARON
  dur = end-start;  //AARON
  std::cout << "*****Matrix asm time:  " << dur.count() << std::endl;    

  start = std::chrono::system_clock::now();
  linearSystem->setupAMG(M);

  end = std::chrono::system_clock::now();  //AARON
  dur = end-start;  //AARON
  std::cout << "*****Matrix prec time:  " << dur.count() << std::endl;
}

LaplaceXY2Hypre::~LaplaceXY2Hypre() {
  delete M;
  delete x;
  delete b;
  delete linearSystem;
}

Field2D LaplaceXY2Hypre::solve(Field2D& rhs, Field2D& x0) {
  Timer timer("invert");

  ASSERT1(rhs.getMesh() == localmesh);
  ASSERT1(x0.getMesh() == localmesh);
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  auto start = std::chrono::system_clock::now();  //AARON

  x->importValuesFromField(x0);
  b->importValuesFromField(rhs);
  x->assemble();
  b->assemble();

  auto form_vec = std::chrono::system_clock::now();  //AARON  
  std::chrono::duration<double> formvec_dur = form_vec-start;  //AARON
  std::cout << "*****Form Vectors time:  " << formvec_dur.count() << std::endl;

  // Solve the system
  start = std::chrono::system_clock::now();  //AARON
  linearSystem->solve();

  auto slv = std::chrono::system_clock::now();  //AARON
  std::chrono::duration<double> slv_dur = slv-start;  //AARON
  std::cout << "*****BoomerAMG solve time:  " << slv_dur.count() << std::endl;

  // Convert result into a Field2D
  start = std::chrono::system_clock::now();  //AARON  
  Field2D sol = x->toField();
  //sol.allocate().setLocation(CELL_LOC::centre);
  //x->exportValuesToField(sol);

  auto formfield = std::chrono::system_clock::now();  //AARON
  std::chrono::duration<double> formfield_dur = formfield-start;  //AARON
  std::cout << "*****Form field time:  " << formfield_dur.count() << std::endl;  

  return sol;
}

#endif // BOUT_HAS_HYPRE
