
#ifdef BOUT_HAS_HYPRE

#include <bout/invert/laplacexy2_hypre.hxx>

#include <bout/assert.hxx>

#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <globals.hxx>
#include <utils.hxx>

#include <output.hxx>

#include <cmath>

#include "_hypre_utilities.h"

Ind2D index2d(Mesh* mesh, int x, int y) {
  int ny = mesh->LocalNy;
  return Ind2D(x * ny + y, ny, 1);
}

LaplaceXY2Hypre::LaplaceXY2Hypre(Mesh* m, Options* opt, const CELL_LOC loc)
    : localmesh(m == nullptr ? bout::globals::mesh : m), f2dinit(localmesh),
      location(loc), matrix(nullptr) {
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
  HYPRE_Init();
  //hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_DEVICE;

  matrix = new bout::HypreMatrix<Field2D>(f2dinit);

  indexConverter = GlobalIndexer::getInstance(f2dinit.getMesh());
  const auto& region = f2dinit.getRegion("RGN_ALL_THIN");

  const auto ilower =
      static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*std::begin(region)));
  const auto iupper =
      static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*--std::end(region)));

  //const MPI_Comm comm = f2dinit.getMesh()->getXcomm();
  const MPI_Comm comm = BoutComm::get();


  std::cout << "Matrix bounds:  " << ilower << ", " << iupper << std::endl;
  HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &ij_matrix);
  HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(ij_matrix);  

  HYPRE_BoomerAMGCreate(&solver);

  // Configure Linear Solver

  // Defaults taken from one of the HYPRE examples (ex5.c)
  const BoutReal tol = (*opt)["tol"].doc("Tolerance").withDefault(1e-7);
  const int relax_type = (*opt)["relax_type"].doc("What smoother to use").withDefault(18);
  const bool cf_relaxation =
      (*opt)["cf_relaxtion"].doc("Use CF-relaxation").withDefault(false);
  const int num_sweeps = (*opt)["num_sweeps"].doc("Number of sweeps").withDefault(1);
  const int max_levels =
      (*opt)["max_levels"].doc("Maximum number of multigrid levels").withDefault(20);

#if CHECKLEVEL >= 1
  HYPRE_BoomerAMGSetPrintLevel(solver, 0);
#endif

  // Falgout coarsening with modified classical interpolaiton
  HYPRE_BoomerAMGSetOldDefault(solver);
  // G-S/Jacobi hybrid relaxation must be 18 or 7 for GPU implementation
  HYPRE_BoomerAMGSetRelaxType(solver, 18);
  // uses C/F relaxation must be false for GPU
  HYPRE_BoomerAMGSetRelaxOrder(solver, false);
  // Coarsening type must be PMIS (8) for GPU implementation
  HYPRE_BoomerAMGSetCoarsenType(solver, 8);
  // Interp type must be 3 or 15 for GPU implementation
  HYPRE_BoomerAMGSetInterpType(solver, 3);  
  // Sweeps on each level
  HYPRE_BoomerAMGSetNumSweeps(solver, num_sweeps);
  // maximum number of levels
  HYPRE_BoomerAMGSetMaxLevels(solver, max_levels);
  // convergence tolerance for 20 iterations for testing
  HYPRE_BoomerAMGSetTol(solver, 0.0);
  HYPRE_BoomerAMGSetMaxIter(solver, 20);
  HYPRE_BoomerAMGSetKeepTranspose(solver, 1);

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

  const auto& region = f2dinit.getRegion("RGN_NOBNDRY");
  HYPRE_BigInt nrows = static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*--std::end(region))) -
                        static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*std::begin(region)));
  HYPRE_BigInt inner_size = 5*nrows;
  HYPRE_BigInt *ncols, *rows, *cols;
  HYPRE_Complex *vals;
  cudaMallocManaged(&ncols, nrows*sizeof(HYPRE_BigInt));
  cudaMallocManaged(&rows, nrows*sizeof(HYPRE_BigInt));
  cudaMallocManaged(&cols, inner_size*sizeof(HYPRE_BigInt));
  cudaMallocManaged(&vals, inner_size*sizeof(HYPRE_Complex));
  for (int i = 0; i < nrows; ++i) {
    ncols[i] = 5;
  }

  Coordinates* coords = localmesh->getCoordinates(location);

  //////////////////////////////////////////////////
  // Set Matrix elements
  //
  // (1/J) d/dx ( J * g11 d/dx ) + (1/J) d/dy ( J * g22 d/dy )
  auto start = std::chrono::system_clock::now();  //AARON
  int rowi = 0;
  for (auto& index : A.getRegion("RGN_NOBNDRY")) {
    // Index offsets
    auto ind_xp = index.xp();
    auto ind_xm = index.xm();
    HYPRE_BigInt I = indexConverter->getGlobal(index);
    HYPRE_BigInt XP = indexConverter->getGlobal(ind_xp);
    HYPRE_BigInt XM = indexConverter->getGlobal(ind_xm);
    rows[rowi] = I;
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

    //matrix(index, ind_xp) = xp;
    //matrix(index, ind_xm) = xm;
    cols[5*rowi] = XP;
    vals[5*rowi] = xp;
    cols[5*rowi+1] = XP;
    vals[5*rowi+1] = xm;

    if (include_y_derivs) {
      auto ind_yp = index.yp();
      auto ind_ym = index.ym();
      HYPRE_BigInt YP = indexConverter->getGlobal(ind_yp);
      HYPRE_BigInt YM = indexConverter->getGlobal(ind_ym);

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
      //matrix(index, ind_yp) = yp;

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
      //matrix(index, ind_yp) = yp;
      //matrix(index, ind_ym) = ym;
      cols[5*rowi+2] = YP;
      vals[5*rowi+2] = yp;
      cols[5*rowi+3] = YM;
      vals[5*rowi+3] = ym;      
    }
    // Note: The central coefficient is done last because this may be modified
    // if y derivs are/are not included.
    //matrix(index, index) = c;
    // manV[0] = c;
    // HYPRE_IJMatrixSetValues(ij_matrix, 1, one, manI, manI, manV);
    cols[5*rowi+4] = I;
    vals[5*rowi+4] = c;
    rowi++;
  }
  HYPRE_IJMatrixSetValues(ij_matrix, nrows, ncols, rows, cols, vals);

  // X boundaries
  if (localmesh->firstX()) {
    nrows = localmesh->yend - localmesh->ystart + 1;
    for (int i = 0; i < nrows; ++i) {
      ncols[i] = 2;
    }
    rowi = 0;
    if (x_inner_dirichlet) {

      // Dirichlet on inner X boundary
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart, y);
        auto ind_xm = index.xm();
        HYPRE_BigInt I = indexConverter->getGlobal(index);
        HYPRE_BigInt XM = indexConverter->getGlobal(ind_xm);
        // matrix(ind_xm, index) = 0.5;
        // matrix(ind_xm, ind_xm) = 0.5;
        rows[rowi] = XM;
        cols[2*rowi] = I;
        vals[2*rowi] = 0.5;
        cols[2*rowi+1] = XM;
        vals[2*rowi+1] = 0.5;
        rowi++;
      }

    } else {

      // Neumann on inner X boundary
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart, y);
        auto ind_xm = index.xm();
        HYPRE_BigInt I = indexConverter->getGlobal(index);
        HYPRE_BigInt XM = indexConverter->getGlobal(ind_xm);
        // matrix(ind_xm, index) = 1.0;
        // matrix(ind_xm, ind_xm) = -1.0;
        rows[rowi] = XM;
        cols[2*rowi] = I;
        vals[2*rowi] = 1.0;
        cols[2*rowi+1] = XM;
        vals[2*rowi+1] = -1.0;
        rowi++;
      }
    }
    HYPRE_IJMatrixSetValues(ij_matrix, nrows, ncols, rows, cols, vals);  
  }
  if (localmesh->lastX()) {
    // Dirichlet on outer X boundary
    nrows = localmesh->yend - localmesh->ystart + 1;
    for (int i = 0; i < nrows; ++i) {
      ncols[i] = 2;
    }
    rowi = 0;
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      auto index = index2d(localmesh, localmesh->xend, y);
      auto ind_xp = index.xp();
      HYPRE_BigInt I = indexConverter->getGlobal(index);
      HYPRE_BigInt XP = indexConverter->getGlobal(ind_xp);        
      // matrix(ind_xp, ind_xp) = 0.5;
      // matrix(ind_xp, index) = 0.5;
      rows[rowi] = XP;
      cols[2*rowi] = XP;
      vals[2*rowi] = 0.5;
      cols[2*rowi+1] = I;
      vals[2*rowi+1] = 0.5;
      rowi++;      
    }
    HYPRE_IJMatrixSetValues(ij_matrix, nrows, ncols, rows, cols, vals);  
  }

  if (y_bndry_dirichlet) {
    // Dirichlet on Y boundaries
    nrows = 0;
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      nrows ++;
    }
    for (int i = 0; i < nrows; ++i) {
      ncols[i] = 2;
    }
    rowi = 0;    
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart);
      auto ind_ym = index.ym();
      HYPRE_BigInt I = indexConverter->getGlobal(index);
      HYPRE_BigInt YM = indexConverter->getGlobal(ind_ym);     
      // matrix(ind_ym, ind_ym) = 0.5;
      // matrix(ind_ym, index) = 0.5;
      rows[rowi] = YM;
      cols[2*rowi] = YM;
      vals[2*rowi] = 0.5;
      cols[2*rowi+1] = I;
      vals[2*rowi+1] = 0.5;
      rowi++;            
    }
    HYPRE_IJMatrixSetValues(ij_matrix, nrows, ncols, rows, cols, vals);  

    nrows = 0;
    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      nrows ++;
    }
    for (int i = 0; i < nrows; ++i) {
      ncols[i] = 2;
    }
    rowi = 0;    
    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->yend);
      auto ind_yp = index.yp();
      HYPRE_BigInt I = indexConverter->getGlobal(index);
      HYPRE_BigInt YP = indexConverter->getGlobal(ind_yp);      
      // matrix(ind_yp, ind_yp) = 0.5;
      // matrix(ind_yp, index) = 0.5;
      rows[rowi] = YP;
      cols[2*rowi] = YP;
      vals[2*rowi] = 0.5;
      cols[2*rowi+1] = I;
      vals[2*rowi+1] = 0.5;
      rowi++;
    }
    HYPRE_IJMatrixSetValues(ij_matrix, nrows, ncols, rows, cols, vals);  
  } else {
    // Neumann on Y boundaries
    nrows = 0;
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      nrows ++;
    }
    for (int i = 0; i < nrows; ++i) {
      ncols[i] = 2;
    }
    rowi = 0;       
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart);
      auto ind_ym = index.ym();
      HYPRE_BigInt I = indexConverter->getGlobal(index);
      HYPRE_BigInt YM = indexConverter->getGlobal(ind_ym);      
      // matrix(ind_ym, ind_ym) = -1.0;
      // matrix(ind_ym, index) = 1.0;
      rows[rowi] = YM;
      cols[2*rowi] = YM;
      vals[2*rowi] = -1.0;
      cols[2*rowi+1] = I;
      vals[2*rowi+1] = 1.0;
      rowi++;      
    }
    HYPRE_IJMatrixSetValues(ij_matrix, nrows, ncols, rows, cols, vals);


    nrows = 0;
    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      nrows ++;
    }
    for (int i = 0; i < nrows; ++i) {
      ncols[i] = 2;
    }
    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->yend);
      auto ind_yp = index.yp();
      HYPRE_BigInt I = indexConverter->getGlobal(index);
      HYPRE_BigInt YP = indexConverter->getGlobal(ind_yp);       
      // matrix(ind_yp, ind_yp) = 1.0;
      // matrix(ind_yp, index) = -1.0;
      rows[rowi] = YP;
      cols[2*rowi] = YP;
      vals[2*rowi] = 1.0;
      cols[2*rowi+1] = I;
      vals[2*rowi+1] = -1.0;
      rowi++;         
    }
    HYPRE_IJMatrixSetValues(ij_matrix, nrows, ncols, rows, cols, vals);
  }
  auto end = std::chrono::system_clock::now();  //AARON
  std::chrono::duration<double> dur = end-start;  //AARON
  std::cout << "*****Matrix set time:  " << dur.count() << std::endl;    


  start = std::chrono::system_clock::now();
  // Assemble Matrix
  HYPRE_IJMatrixAssemble(ij_matrix);
  HYPRE_IJMatrixGetObject(ij_matrix, (void**) &parcsr_matrix);

  end = std::chrono::system_clock::now();  //AARON
  dur = end-start;  //AARON
  std::cout << "*****Matrix asm time:  " << dur.count() << std::endl;    

  start = std::chrono::system_clock::now();
  // Set the operator
  HYPRE_BoomerAMGSetup(solver, parcsr_matrix, nullptr, nullptr);

  end = std::chrono::system_clock::now();  //AARON
  dur = end-start;  //AARON
  std::cout << "*****Matrix prec time:  " << dur.count() << std::endl;

  cudaFree(rows);
  cudaFree(ncols);
  cudaFree(cols);
  cudaFree(vals);
}

LaplaceXY2Hypre::~LaplaceXY2Hypre() {
  if (solver != nullptr) {
    HYPRE_BoomerAMGDestroy(solver);
  }

  if (matrix != nullptr) {
    delete matrix;
  }
}

const Field2D LaplaceXY2Hypre::solve(const Field2D& rhs, const Field2D& x0) {
  Timer timer("invert");

  ASSERT1(rhs.getMesh() == localmesh);
  ASSERT1(x0.getMesh() == localmesh);
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  // Load initial guess x0 into xs and rhs into bs
  HYPRE_IJVector ij_x, ij_b;
  HYPRE_ParVector par_x, par_b;
  const MPI_Comm comm = BoutComm::get();
  //const MPI_Comm comm = f2dinit.getMesh()->getXcomm();
  const auto& region = rhs.getRegion("RGN_ALL_THIN");
  const auto jlower =
        static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*std::begin(region)));
  const auto jupper =
        static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*--std::end(region)));  
  HYPRE_IJVectorCreate(comm, jlower, jupper, &ij_x);
  HYPRE_IJVectorCreate(comm, jlower, jupper, &ij_b);
  HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
  HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(ij_x);
  HYPRE_IJVectorInitialize(ij_b);

  HYPRE_BigInt nrows = jupper - jlower + 1;
  HYPRE_BigInt count;
  HYPRE_BigInt *manI, *rows;
  HYPRE_Complex *manV, *vals_x, *vals_b;
  cudaMallocManaged(&manI, sizeof(HYPRE_BigInt));
  cudaMallocManaged(&manV, sizeof(HYPRE_Complex));
  cudaMallocManaged(&rows, nrows*sizeof(HYPRE_BigInt));
  cudaMallocManaged(&vals_x, nrows*sizeof(HYPRE_Complex));
  cudaMallocManaged(&vals_b, nrows*sizeof(HYPRE_Complex));
  IndexerPtr indexConverter = GlobalIndexer::getInstance(rhs.getMesh());

  auto start = std::chrono::system_clock::now();  //AARON
  count = 0;
  BOUT_FOR_SERIAL(i, region) {
    const auto index = static_cast<HYPRE_BigInt>(indexConverter->getGlobal(i));
    if (index != -1) {
      rows[count] = index;
      vals_x[count] = x0[i];
      vals_b[count] = rhs[i];
      count ++;   
    }
  }
  HYPRE_IJVectorSetValues(ij_x, count, rows, vals_x);
  HYPRE_IJVectorSetValues(ij_b, count, rows, vals_b);

  count = 0;
  if (localmesh->firstX()) {
    if (x_inner_dirichlet) {
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart - 1, y);
        rows[count] = indexConverter->getGlobal(index);
        vals_x[count] = x0[index];
        vals_b[count] = 0.5 * (x0[index] + x0[index.xp()]);
        count ++;       
      }
    } else {
      // Inner X boundary (Neumann)
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        auto index = index2d(localmesh, localmesh->xstart - 1, y);
        rows[count] = indexConverter->getGlobal(index);
        vals_x[count] = x0[index];
        vals_b[count] = 0.0;
        count ++;  
      }
    }
  }
  HYPRE_IJVectorSetValues(ij_x, count, rows, vals_x);
  HYPRE_IJVectorSetValues(ij_b, count, rows, vals_b);  

  // Outer X boundary (Dirichlet)
  if (localmesh->lastX()) {
    count = 0;
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      auto index = index2d(localmesh, localmesh->xend + 1, y);
      rows[count] = indexConverter->getGlobal(index);
      vals_x[count] = x0[index];
      vals_b[count] = 0.5 * (x0[index.xm()] + x0[index]);
      count ++;        
    }
    HYPRE_IJVectorSetValues(ij_x, count, rows, vals_x);
    HYPRE_IJVectorSetValues(ij_b, count, rows, vals_b);  
  }

  count = 0;
  if (y_bndry_dirichlet) {

    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart - 1);
      rows[count] = indexConverter->getGlobal(index);
      vals_x[count] = x0[index];
      vals_b[count] = 0.5 * (x0[index] + x0[index.yp()]);
      count ++;        
    }
  
    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->yend + 1);
      rows[count] = indexConverter->getGlobal(index);
      vals_x[count] = x0[index];
      vals_b[count] = 0.5 * (x0[index] + x0[index.xm()]);
      count ++;        
    }
  } else {
    // Y boundaries Neumann
    for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->ystart - 1);
      rows[count] = indexConverter->getGlobal(index);
      vals_x[count] = x0[index];
      vals_b[count] = 0.0;
      count ++;   
    }

    for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
      auto index = index2d(localmesh, it.ind, localmesh->yend + 1);
      rows[count] = indexConverter->getGlobal(index);
      vals_x[count] = x0[index];
      vals_b[count] = 0.0;
      count ++;         
    }
  }
  HYPRE_IJVectorSetValues(ij_x, count, rows, vals_x);
  HYPRE_IJVectorSetValues(ij_b, count, rows, vals_b);      

  HYPRE_IJVectorAssemble(ij_x);
  HYPRE_IJVectorAssemble(ij_b);
  HYPRE_IJVectorGetObject(ij_x, reinterpret_cast<void**>(&par_x));
  HYPRE_IJVectorGetObject(ij_b, reinterpret_cast<void**>(&par_b));

  auto form_vec = std::chrono::system_clock::now();  //AARON  
  std::chrono::duration<double> formvec_dur = form_vec-start;  //AARON
  std::cout << "*****Form Vectors time:  " << formvec_dur.count() << std::endl;

  // Solve the system
  start = std::chrono::system_clock::now();  //AARON
  HYPRE_BoomerAMGSolve(solver, parcsr_matrix, par_b, par_x);

  auto slv = std::chrono::system_clock::now();  //AARON
  std::chrono::duration<double> slv_dur = slv-start;  //AARON
  std::cout << "*****BoomerAMG solve time:  " << slv_dur.count() << std::endl;

  HYPRE_Int numiter;
  HYPRE_BoomerAMGGetNumIterations(solver, &numiter);
  std::cout << "*****BoomerAMG num iterations:  " << numiter << std::endl;


  // Convert result into a Field2D
  start = std::chrono::system_clock::now();  //AARON  
  Field2D x;
  x.allocate().setLocation(CELL_LOC::centre);
  count = 0;
  BOUT_FOR_SERIAL(i, x0.getRegion("RGN_ALL_THIN")) {
    const auto index = static_cast<HYPRE_BigInt>(indexConverter->getGlobal(i));
    if (index != -1) {
      rows[count] = index;
      count ++;
    }
  }
  HYPRE_IJVectorGetValues(ij_x, count, rows, vals_x);
  count = 0;
  BOUT_FOR_SERIAL(i, x0.getRegion("RGN_ALL_THIN")) {
    const auto index = static_cast<HYPRE_BigInt>(indexConverter->getGlobal(i));
    if (index != -1) {
      x[i] = static_cast<BoutReal>(vals_x[count]);
      count ++;
    }
  }  

  auto formfield = std::chrono::system_clock::now();  //AARON
  std::chrono::duration<double> formfield_dur = formfield-start;  //AARON
  std::cout << "*****Form field time:  " << formfield_dur.count() << std::endl;  

  cudaFree(vals_x);
  cudaFree(vals_b);  
  HYPRE_Finalize();
  return x;
}

#endif // BOUT_HAS_HYPRE
