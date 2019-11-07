/**************************************************************************
 * 3D Laplacian Solver
 *                           Using PETSc Solvers
 *
 **************************************************************************
 * Copyright 2013 J. Buchanan, J.Omotani
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/
#ifdef BOUT_HAS_PETSC

#include "petsc3damg.hxx"

#include <bout/mesh.hxx>
#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <bout/assert.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <bout/petsc_interface.hxx>

constexpr auto KSP_RICHARDSON = "richardson";
constexpr auto KSP_CHEBYSHEV  = "chebyshev";
constexpr auto KSP_CG         = "cg";
constexpr auto KSP_GMRES      = "gmres";
constexpr auto KSP_TCQMR      = "tcqmr";
constexpr auto KSP_BCGS       = "bcgs";
constexpr auto KSP_CGS        = "cgs";
constexpr auto KSP_TFQMR      = "tfqmr";
constexpr auto KSP_CR         = "cr";
constexpr auto KSP_LSQR       = "lsqr";
constexpr auto KSP_BICG       = "bicg";
constexpr auto KSP_PREONLY    = "preonly";

LaplacePetsc3dAmg::LaplacePetsc3dAmg(Options *opt, const CELL_LOC loc, Mesh *mesh_in) :
  Laplacian(opt, loc, mesh_in),
  A(0.0), C1(1.0), C2(1.0), D(1.0), Ex(0.0), Ez(0.0),
  issetD(false), issetC(false), issetE(false), updateRequired(true), 
  operator3D(A), kspInitialised(false)
{

  // Provide basic initialisation of field coefficients, etc.
  // Get relevent options from user input
  // Initialise PETSc objects
  A.setLocation(location);
  C1.setLocation(location);
  C2.setLocation(location);
  D.setLocation(location);
  Ex.setLocation(location);
  Ez.setLocation(location);

  // Get Options in Laplace Section
  if (!opt) {
    opts = Options::getRoot()->getSection("laplace");
  } else {
    opts=opt;
  }

  // Get y boundary flags
  lower_boundary_flags = (*opts)["lower_boundary_flags"].withDefault(0);
  upper_boundary_flags = (*opts)["upper_boundary_flags"].withDefault(0);

  #if CHECK > 0
    // These are the implemented flags
    implemented_flags = INVERT_START_NEW;
    implemented_boundary_flags = INVERT_AC_GRAD
                                 + INVERT_SET
                                 + INVERT_RHS
                                 ;
    // Checking flags are set to something which is not implemented
    // This is done binary (which is possible as each flag is a power of 2)
    if ( global_flags & ~implemented_flags ) {
      if (global_flags&INVERT_4TH_ORDER) output<<"For PETSc based Laplacian inverter, use 'fourth_order=true' instead of setting INVERT_4TH_ORDER flag"<<endl;
      throw BoutException("Attempted to set Laplacian inversion flag that is not implemented in petsc_laplace.cxx");
    }
    if ( inner_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }
    if ( outer_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }
    if ( lower_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }
    if ( upper_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }    
    if(localmesh->periodicX) {
      throw BoutException("LaplacePetsc3dAmg does not work with periodicity in the x direction (localmesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
    }
  #endif

  // Get 4th order solver switch
  fourth_order = (*opts)["fourth_order"].withDefault(false);

  // Get Tolerances for KSP solver
  rtol = (*opts)["rtol"].doc("Relative tolerance for KSP solver").withDefault(1e-5);
  atol = (*opts)["atol"].doc("Absolute tolerance for KSP solver").withDefault(1e-5);
  dtol = (*opts)["dtol"].doc("Divergence tolerance for KSP solver").withDefault(1e5);
  maxits = (*opts)["maxits"].doc("Maximum number of KSP iterations").withDefault(100000);

  richardson_damping_factor = (*opts)["richardson_damping_factor"].withDefault(1.0);
  chebyshev_max = (*opts)["chebyshev_max"].withDefault(100.0);
  chebyshev_min = (*opts)["chebyshev_min"].withDefault(0.01);
  gmres_max_steps = (*opts)["gmres_max_steps"].withDefault(30);

  // Get KSP Solver Type (Generalizes Minimal RESidual is the default)
  ksptype = (*opts)["ksptype"].doc("KSP solver type").withDefault(KSP_GMRES);

  // Get direct solver switch
  direct = (*opts)["direct"].doc("Use direct (LU) solver?").withDefault(false);
  if (direct) {
    output << endl << "Using LU decompostion for direct solution of system" << endl << endl;
  }

  // Set up boundary conditions in operator
  BOUT_FOR(i, localmesh->getRegion3D("RGN_INNER_X_THIN")) {
    if(inner_boundary_flags & INVERT_AC_GRAD) {
      // Neumann on inner X boundary
      operator3D(i, i) = -1./coords->dx[i]/sqrt(coords->g_11[i]);
      operator3D(i, i.xp()) = 1./coords->dx[i]/sqrt(coords->g_11[i]);
    } else {
      // Dirichlet on inner X boundary
      operator3D(i, i) = 0.5;
      operator3D(i, i.xp()) = 0.5;
    }
  }

  BOUT_FOR(i, localmesh->getRegion3D("RGN_OUTER_X_THIN")) {
    if(outer_boundary_flags & INVERT_AC_GRAD) {
      // Neumann on outer X boundary
      operator3D(i, i) = 1./coords->dx[i]/sqrt(coords->g_11[i]);
      operator3D(i, i.xm()) = -1./coords->dx[i]/sqrt(coords->g_11[i]);
    } else {
      // Dirichlet on outer X boundary
      operator3D(i, i) = 0.5;
      operator3D(i, i.xm()) = 0.5;
    }
  }

  BOUT_FOR(i, localmesh->getRegion3D("RGN_LOWER_Y_THIN")) {
    if(lower_boundary_flags & INVERT_AC_GRAD) {
      // Neumann on lower Y boundary
      operator3D(i, i) = -1./coords->dy[i]/sqrt(coords->g_22[i]);
      operator3D(i, i.yp()) = 1./coords->dy[i]/sqrt(coords->g_22[i]);
    } else {
      // Dirichlet on lower Y boundary
      operator3D(i, i) = 0.5;
      operator3D(i, i.yp()) = 0.5;
    }
  }

  BOUT_FOR(i, localmesh->getRegion3D("RGN_UPPER_Y_THIN")) {
    if(upper_boundary_flags & INVERT_AC_GRAD) {
      // Neumann on upper Y boundary
      operator3D(i, i) = 1./coords->dy[i]/sqrt(coords->g_22[i]);
      operator3D(i, i.ym()) = -1./coords->dy[i]/sqrt(coords->g_22[i]);
    } else {
      // Dirichlet on upper Y boundary
      operator3D(i, i) = 0.5;
      operator3D(i, i.ym()) = 0.5;
    }
  }
}


LaplacePetsc3dAmg::~LaplacePetsc3dAmg() {
  if (kspInitialised) KSPDestroy(&ksp);
}


Field3D LaplacePetsc3dAmg::solve(const Field3D &b_in, const Field3D &x0) {
  // If necessary, update the values in the matrix operator and initialise
  // the Krylov solver
  if (updateRequired) updateMatrix3D();
  PetscVector<Field3D> rhs(b_in), guess(x0);
  
  // Adjust vectors to represent boundary conditions
  BOUT_FOR(i, localmesh->getRegion3D("RGN_INNER_X_THIN")) {
    const BoutReal val = (inner_boundary_flags & INVERT_SET) ? x0[i] : 0.;
    if (!(inner_boundary_flags & INVERT_RHS)) {
      rhs(i) = val;
    }
  }

  BOUT_FOR(i, localmesh->getRegion3D("RGN_OUTER_X_THIN")) {
    const BoutReal val = (outer_boundary_flags & INVERT_SET) ? x0[i] : 0.;
    if (!(outer_boundary_flags & INVERT_RHS)) {
      rhs(i) = val;
    }
  }

  BOUT_FOR(i, localmesh->getRegion3D("RGN_LOWER_Y_THIN")) {
    const BoutReal val = (lower_boundary_flags & INVERT_SET) ? x0[i] : 0.;
    if (!(lower_boundary_flags & INVERT_RHS)) {
      rhs(i) = val;
    }
  }

  BOUT_FOR(i, localmesh->getRegion3D("RGN_UPPER_Y_THIN")) {
    const BoutReal val = (upper_boundary_flags & INVERT_SET) ? x0[i] : 0.;
    if (!(upper_boundary_flags & INVERT_RHS)) {
      rhs(i) = val;
    }
  }

  rhs.assemble();
  guess.assemble();

  PetscViewer view, view2;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "rhs.out", &view);
  VecView(*rhs.getVectorPointer(), view);
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "guess.out", &view2);
  VecView(*guess.getVectorPointer(), view2);
  
  // Invoke solver
  { Timer timer("petscsolve");
    KSPSolve(ksp, *rhs.getVectorPointer(), *guess.getVectorPointer());
  }

  // Check for convergence
  KSPConvergedReason reason;
  KSPGetConvergedReason( ksp, &reason );
  if (reason==-3) { // Too many iterations, might be fixed by taking smaller timestep
    throw BoutIterationFail("Petsc3dAmg: too many iterations");
  }
  else if (reason<=0) {
    output << "KSPConvergedReason is " << reason << endl;
    throw BoutException("Petsc3dAmg: inversion failed to converge.");
  }

  // Create field from result
  return guess.toField();
}

Field2D LaplacePetsc3dAmg::solve(const Field2D &b) {
  return Laplacian::solve(b);
}

PetscMatrix<Field3D>& LaplacePetsc3dAmg::getMatrix3D() {
  if (updateRequired) updateMatrix3D();
  return operator3D;
}

void LaplacePetsc3dAmg::updateMatrix3D() {
  localmesh->communicate(C2);
  const Field3D dc_dx = issetC ? DDX(C2) : Field3D(),
    dc_dy = issetC ? DDY(C2) : Field3D(),
    dc_dz = issetC ? DDZ(C2) : Field3D(),
    dJ_dy = DDY(coords->J/coords->g_22);

  // Set up the matrix for the internal points on the grid.
  // Boundary conditions were set in the constructor.
  BOUT_FOR(l, localmesh->getRegion3D("RGN_NOBNDRY")) {
    // Index is called l for "location". It is not called i so as to
    // avoid confusing it with the x-index.

    // Calculate coefficients for the terms in the differential operator
    BoutReal C_df_dx = coords->G1[l], C_df_dz = coords->G3[l];
    if (issetD) {
      C_df_dx *= D[l];
      C_df_dz *= D[l];
    }
    if (issetC) {
      C_df_dx += (coords->g11[l]*dc_dx[l] + coords->g12[l]*dc_dy[l] +
		  coords->g13[l]*dc_dz[l])/C1[l];
      C_df_dz += (coords->g13[l]*dc_dx[l] + coords->g23[l]*dc_dy[l] +
		  coords->g33[l]*dc_dz[l])/C1[l];
    }
    if (issetE) {
      C_df_dx += Ex[l];
      C_df_dz += Ez[l];
    }

    BoutReal C_d2f_dx2 = coords->g11[l],
      C_d2f_dy2 = (coords->g22[l] - 1.0/coords->g_22[l]),
      C_d2f_dz2 = coords->g33[l];
    if (issetD) {
      C_d2f_dx2 *= D[l];
      C_d2f_dy2 *= D[l];
      C_d2f_dz2 *= D[l];
    }
  
    BoutReal C_d2f_dxdz = 2*coords->g13[l];
    if (issetD) {
      C_d2f_dxdz *= D[l];
    }
  
    // Adjust the coefficients to include finite-difference factors
    if (nonuniform) {
      C_df_dx -= C_d2f_dx2*coords->d1_dx[l];
    }
    C_df_dx /= 2*coords->dx[l];
    C_df_dz /= 2*coords->dz;
    
    C_d2f_dx2 /= SQ(coords->dx[l]);
    C_d2f_dy2 /= SQ(coords->dy[l]);
    C_d2f_dz2 /= SQ(coords->dz);
    
    C_d2f_dxdz /= 4*coords->dx[l]*coords->dz;
    
    operator3D(l, l) = -2*(C_d2f_dx2 + C_d2f_dy2 + C_d2f_dz2) + A[l];
    operator3D(l, l.xp()) = C_df_dx + C_d2f_dx2;
    operator3D(l, l.xm()) = -C_df_dx + C_d2f_dx2;
    operator3D(l, l.zp()) = C_df_dz + C_d2f_dz2;
    operator3D(l, l.zm()) = -C_df_dz + C_d2f_dz2;
    operator3D(l, l.xp().zp()) = C_d2f_dxdz;
    operator3D(l, l.xp().zm()) = -C_d2f_dxdz;
    operator3D(l, l.xm().zp()) = -C_d2f_dxdz;
    operator3D(l, l.xm().zm()) = C_d2f_dxdz;
  }
  operator3D.partialAssemble();

  // Must add these (rather than assign) so that elements used in
  // interpolation don't overwrite each other.
  BOUT_FOR(l, localmesh->getRegion3D("RGN_NOBNDRY")) {
    BoutReal C_df_dy = (coords->G2[l] - dJ_dy[l]/coords->J[l]);
    if (issetD) {
      C_df_dy *= D[l];
    }
    if (issetC) {
      C_df_dy += (coords->g12[l]*dc_dx[l] + (coords->g22[l] - 1./coords->g_22[l])*dc_dy[l] +
		  coords->g23[l]*dc_dz[l])/C1[l];
    }

    BoutReal C_d2f_dy2 = (coords->g22[l] - 1.0/coords->g_22[l]);
    if (issetD) {
      C_d2f_dy2 *= D[l];
    }
  
    BoutReal C_d2f_dxdy = 2*coords->g12[l],
      C_d2f_dydz = 2*coords->g23[l];
    if (issetD) {
      C_d2f_dxdy *= D[l];
      C_d2f_dydz *= D[l];
    }
  
    // Adjust the coefficients to include finite-difference factors
    if (nonuniform) {
      C_df_dy -= C_d2f_dy2*coords->d1_dy[l];
    }
    C_df_dy /= 2*coords->dy[l];
    C_d2f_dy2 /= SQ(coords->dy[l]);
    C_d2f_dxdy /= 4*coords->dy[l]; // NOTE: This is unfinished; will also need to
                                   // divide by dx(i +/- 1, j, k) when using
    C_d2f_dydz /= 4*coords->dy[l]*coords->dz;

    operator3D.yup()(l, l.yp()) += C_df_dy + C_d2f_dy2;
    operator3D.ydown()(l, l.ym()) += -C_df_dy + C_d2f_dy2;
    operator3D.yup()(l, l.xp().yp()) += C_d2f_dxdy/coords->dy[l.xp()];
    operator3D.ydown()(l, l.xp().ym()) += -C_d2f_dxdy/coords->dy[l.xp()];
    operator3D.yup()(l, l.xm().yp()) += -C_d2f_dxdy/coords->dy[l.xm()];
    operator3D.ydown()(l, l.xm().ym()) += C_d2f_dxdy/coords->dy[l.xm()];
    operator3D.yup()(l, l.yp().zp()) += C_d2f_dydz;
    operator3D.yup()(l, l.yp().zm()) += -C_d2f_dydz;
    operator3D.ydown()(l, l.ym().zp()) += -C_d2f_dydz;
    operator3D.ydown()(l, l.ym().zm()) += C_d2f_dydz;    
  }
  operator3D.assemble();
  PetscViewer view;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix.out", &view);
  MatView(*operator3D.getMatrixPointer(), view);
  
  // Declare KSP Context (abstract PETSc object that manages all Krylov methods)
  if (kspInitialised) KSPDestroy(&ksp);
  KSPCreate(BoutComm::get(), &ksp);
  kspInitialised = true;
#if PETSC_VERSION_GE(3,5,0)
  KSPSetOperators(ksp, *operator3D.getMatrixPointer(), *operator3D.getMatrixPointer());
#else
  KSPSetOperators(ksp, *operator3D.getMatrixPointer(), *operator3D.getMatrixPointer(),
		  DIFFERENT_NONZERO_PATTERN );
#endif

  PC pc;
  KSPGetPC(ksp,&pc);

  if (direct) {
  // Set the type of the preconditioner
    PCSetType(pc, PCLU);
    KSPSetType(ksp, KSP_PREONLY);
#if PETSC_VERSION_GE(3,9,0)
    PCFactorSetMatSolverType(pc,"mumps");
#else
    PCFactorSetMatSolverPackage(pc,"mumps");
#endif
  } else {
    KSPSetType(ksp, ksptype.c_str()); // Set the type of the solver
  
    if(ksptype == KSPRICHARDSON) KSPRichardsonSetScale( ksp, richardson_damping_factor );
#ifdef KSPCHEBYSHEV
    else if (ksptype == KSPCHEBYSHEV) KSPChebyshevSetEigenvalues(ksp, chebyshev_max, chebyshev_min);
#endif
    else if (ksptype == KSPGMRES) KSPGMRESSetRestart(ksp, gmres_max_steps);
  
    // Set the relative and absolute tolerances
    KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
  
    // If the initial guess is not set to zero
    if (!(global_flags & INVERT_START_NEW)) KSPSetInitialGuessNonzero(ksp, (PetscBool) true);
    
    // Set the relative and absolute tolerances
    KSPSetFromOptions(ksp);
    PCSetType(pc, PCGAMG);
    PCGAMGSetType(pc, PCGAMGAGG); //TODO: DETERMINE IF THIS IS MOST APPROPRIATE TYPE OF SOLVER
  }
  KSPSetPCSide(ksp, PC_LEFT);

  updateRequired = false;
}

#endif // BOUT_HAS_PETSC_3_3
