
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
#include "bout/bout_types.hxx"
#include "bout/build_config.hxx"

#if BOUT_HAS_PETSC

#include "petsc3damg.hxx"

#include <bout/assert.hxx>
#include <bout/boutcomm.hxx>
#include <bout/derivs.hxx>
#include <bout/mesh.hxx>
#include <bout/operatorstencil.hxx>
#include <bout/petsc_interface.hxx>
#include <bout/sys/timer.hxx>
#include <bout/utils.hxx>

using bout::utils::flagSet;

#ifdef PETSC_HAVE_HYPRE
static constexpr auto DEFAULT_PC_TYPE = PCHYPRE;
#else
static constexpr auto DEFAULT_PC_TYPE = PCGAMG;
#endif // PETSC_HAVE_HYPRE

LaplacePetsc3dAmg::LaplacePetsc3dAmg(Options* opt, const CELL_LOC loc, Mesh* mesh_in,
                                     Solver* UNUSED(solver))
    : Laplacian(opt, loc, mesh_in), A(0.0), C1(1.0), C2(1.0), D(1.0), Ex(0.0), Ez(0.0),
      opts(opt == nullptr ? &(Options::root()["laplace"]) : opt),
      lower_boundary_flags((*opts)["lower_boundary_flags"].withDefault(0)),
      upper_boundary_flags((*opts)["upper_boundary_flags"].withDefault(0)),
      ksptype((*opts)["ksptype"].doc("KSP solver type").withDefault(KSPGMRES)),
      pctype((*opts)["pctype"].doc("PC type").withDefault(DEFAULT_PC_TYPE)),
      richardson_damping_factor((*opts)["richardson_damping_factor"].withDefault(1.0)),
      chebyshev_max((*opts)["chebyshev_max"].withDefault(100.0)),
      chebyshev_min((*opts)["chebyshev_min"].withDefault(0.01)),
      gmres_max_steps((*opts)["gmres_max_steps"].withDefault(30)),
      rtol((*opts)["rtol"].doc("Relative tolerance for KSP solver").withDefault(1e-5)),
      atol((*opts)["atol"].doc("Absolute tolerance for KSP solver").withDefault(1e-5)),
      dtol((*opts)["dtol"].doc("Divergence tolerance for KSP solver").withDefault(1e6)),
      maxits(
          (*opts)["maxits"].doc("Maximum number of KSP iterations").withDefault(100000)),
      direct((*opts)["direct"].doc("Use direct (LU) solver?").withDefault(false)),
      lowerY(localmesh->iterateBndryLowerY()), upperY(localmesh->iterateBndryUpperY()),
      indexer(std::make_shared<GlobalIndexer<Field3D>>(
          localmesh, getStencil(localmesh, lowerY, upperY))),
      operator3D(indexer), lib(opts) {

  // Provide basic initialisation of field coefficients, etc.
  // Get relevent options from user input
  // Initialise PETSc objects
  A.setLocation(location);
  C1.setLocation(location);
  C2.setLocation(location);
  D.setLocation(location);
  Ex.setLocation(location);
  Ez.setLocation(location);

#if CHECK > 0
  // Checking flags are set to something which is not implemented
  // This is done binary (which is possible as each flag is a power of 2)
  if (isGlobalFlagSet(INVERT_4TH_ORDER)) {
    output.write("For PETSc based Laplacian inverter, use 'fourth_order=true' instead of "
                 "setting INVERT_4TH_ORDER flag\n");
  }

  if (isGlobalFlagSet(~implemented_flags)) {
    throw BoutException("Attempted to set global Laplacian inversion flag that is not "
                        "implemented in petsc_laplace.cxx");
  }

  auto unimplementedBoundaryFlag = [](int boundary_flag,
                                      const std::string& name) -> void {
    if (flagSet(boundary_flag, ~implemented_boundary_flags)) {
      throw BoutException("Attempted to set Laplacian inversion {} boundary flag "
                          "that is not implemented in petsc3damg.cxx",
                          name);
    }
  };
  unimplementedBoundaryFlag(getInnerBoundaryFlags(), "inner");
  unimplementedBoundaryFlag(getOuterBoundaryFlags(), "outer");
  unimplementedBoundaryFlag(lower_boundary_flags, "lower");
  unimplementedBoundaryFlag(upper_boundary_flags, "upper");

  if (localmesh->periodicX) {
    throw BoutException("LaplacePetsc3dAmg does not work with periodicity in the x "
                        "direction (localmesh->PeriodicX == true). Change boundary "
                        "conditions or use serial-tri or cyclic solver instead");
  }
#endif

  if (direct) {
    output.write("\nUsing LU decompostion for direct solution of system\n\n");
  }

  // Set up boundary conditions in operator
  const bool inner_X_neumann = isInnerBoundaryFlagSet(INVERT_AC_GRAD);
  const auto inner_X_BC = inner_X_neumann ? -1. / coords->dx / sqrt(coords->g_11) : 0.5;
  const auto inner_X_BC_plus = inner_X_neumann ? -inner_X_BC : 0.5;

  BOUT_FOR_SERIAL(i, indexer->getRegionInnerX()) {
    operator3D(i, i) = inner_X_BC[i];
    operator3D(i, i.xp()) = inner_X_BC_plus[i];
  }

  const bool outer_X_neumann = isOuterBoundaryFlagSet(INVERT_AC_GRAD);
  const auto outer_X_BC = outer_X_neumann ? 1. / coords->dx / sqrt(coords->g_11) : 0.5;
  const auto outer_X_BC_minus = outer_X_neumann ? -outer_X_BC : 0.5;

  BOUT_FOR_SERIAL(i, indexer->getRegionOuterX()) {
    operator3D(i, i) = outer_X_BC[i];
    operator3D(i, i.xm()) = outer_X_BC_minus[i];
  }

  const bool lower_Y_neumann = flagSet(lower_boundary_flags, INVERT_AC_GRAD);
  const auto lower_Y_BC = lower_Y_neumann ? -1. / coords->dy / sqrt(coords->g_22) : 0.5;
  const auto lower_Y_BC_plus = lower_Y_neumann ? -lower_Y_BC : 0.5;

  BOUT_FOR_SERIAL(i, indexer->getRegionLowerY()) {
    operator3D(i, i) = lower_Y_BC[i];
    operator3D(i, i.yp()) = lower_Y_BC_plus[i];
  }

  const bool upper_Y_neumann = flagSet(upper_boundary_flags, INVERT_AC_GRAD);
  const auto upper_Y_BC = upper_Y_neumann ? 1. / coords->dy / sqrt(coords->g_22) : 0.5;
  const auto upper_Y_BC_minus = upper_Y_neumann ? -upper_Y_BC : 0.5;

  BOUT_FOR_SERIAL(i, indexer->getRegionUpperY()) {
    operator3D(i, i) = upper_Y_BC[i];
    operator3D(i, i.ym()) = upper_Y_BC_minus[i];
  }
}

LaplacePetsc3dAmg::~LaplacePetsc3dAmg() {
  if (kspInitialised) {
    KSPDestroy(&ksp);
  }
}

void setBC(PetscVector<Field3D>& rhs, const Field3D& b_in,
           const Region<Field3D::ind_type>& region, int boundary_flags,
           const Field3D& x0) {
  if (flagSet(boundary_flags, INVERT_RHS)) {
    BOUT_FOR(index, region) { ASSERT1(std::isfinite(b_in[index])); }
  } else {
    const auto& outer_X_BC = (flagSet(boundary_flags, INVERT_SET)) ? x0 : 0.0;
    BOUT_FOR_SERIAL(index, region) { rhs(index) = outer_X_BC[index]; }
  }
}

Field3D LaplacePetsc3dAmg::solve(const Field3D& b_in, const Field3D& x0) {
  AUTO_TRACE();

  // Timing reported in the log files. Includes any matrix construction.
  // The timing for just the solve phase can be retrieved from the "petscsolve"
  // timer if desired.
  const Timer timer("invert");

  // If necessary, update the values in the matrix operator and initialise
  // the Krylov solver
  if (updateRequired) {
    updateMatrix3D();
  }
  PetscVector<Field3D> rhs(b_in, indexer);
  PetscVector<Field3D> guess(x0, indexer);

  // Adjust vectors to represent boundary conditions and check that
  // boundary cells are finite
  setBC(rhs, b_in, indexer->getRegionInnerX(), getInnerBoundaryFlags(), x0);
  setBC(rhs, b_in, indexer->getRegionOuterX(), getOuterBoundaryFlags(), x0);
  setBC(rhs, b_in, indexer->getRegionLowerY(), lower_boundary_flags, x0);
  setBC(rhs, b_in, indexer->getRegionUpperY(), upper_boundary_flags, x0);

  rhs.assemble();
  guess.assemble();

  // Invoke solver
  {
    const Timer timer("petscsolve");
    KSPSolve(ksp, *rhs.get(), *guess.get());
  }

  // Check for convergence
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  if (reason == KSP_DIVERGED_ITS) {
    // Too many iterations, might be fixed by taking smaller timestep
    throw BoutIterationFail("Petsc3dAmg: too many iterations");
  }
  if (reason <= 0) {
    throw BoutException(
        "Petsc3dAmg: inversion failed to converge. KSPConvergedReason: {} ({})",
        KSPConvergedReasons[reason], static_cast<int>(reason));
  }

  // Create field from result
  Field3D solution = guess.toField();
  localmesh->communicate(solution);
  if (solution.hasParallelSlices()) {
    BOUT_FOR(i, indexer->getRegionLowerY()) { solution.ydown()[i] = solution[i]; }
    BOUT_FOR(i, indexer->getRegionUpperY()) { solution.yup()[i] = solution[i]; }
    for (int boundary = 1; boundary < localmesh->ystart; boundary++) {
      BOUT_FOR(i, indexer->getRegionLowerY()) {
        solution.ydown(boundary)[i.ym(boundary)] = solution[i];
      }
      BOUT_FOR(i, indexer->getRegionUpperY()) {
        solution.yup(boundary)[i.yp(boundary)] = solution[i];
      }
    }
  }

  // Set inner boundary cells by extrapolating
  // from single boundary cell which is set by the solver
  // Note: RegionInnerX is the set of points just outside the domain
  //       (in the first boundary cell) so one boundary cell is already set
  BOUT_FOR(i, indexer->getRegionInnerX()) {
    for (int boundary = 1; boundary < localmesh->xstart; boundary++) {
      solution[i.xm(boundary)] = 3. * solution[i.xm(boundary - 1)]
                                 - 3. * solution[i.xm(boundary - 2)]
                                 + solution[i.xm(boundary - 3)];
    }
  }

  // Set outer boundary cells by extrapolating
  // Note: RegionOuterX is the set of points just outside the domain
  //       (in the first boundary cell) so one boundary cell is already set
  BOUT_FOR(i, indexer->getRegionOuterX()) {
    for (int boundary = 1; boundary < localmesh->xstart; boundary++) {
      solution[i.xp(boundary)] = 3. * solution[i.xp(boundary - 1)]
                                 - 3. * solution[i.xp(boundary - 2)]
                                 + solution[i.xp(boundary - 3)];
    }
  }

  checkData(solution);

  return solution;
}

Field2D LaplacePetsc3dAmg::solve(const Field2D& b) { return Laplacian::solve(b); }

PetscMatrix<Field3D>& LaplacePetsc3dAmg::getMatrix3D() {
  if (updateRequired) {
    updateMatrix3D();
  }
  return operator3D;
}

void LaplacePetsc3dAmg::updateMatrix3D() {
  const Field3D dc_dx = issetC ? DDX(C2) : Field3D();
  const Field3D dc_dy = issetC ? DDY(C2) : Field3D();
  const Field3D dc_dz = issetC ? DDZ(C2) : Field3D();
  const auto dJ_dy = DDY(coords->J / coords->g_22);

  // Set up the matrix for the internal points on the grid.
  // Boundary conditions were set in the constructor.
  BOUT_FOR_SERIAL(l, indexer->getRegionNobndry()) {
    // Index is called l for "location". It is not called i so as to
    // avoid confusing it with the x-index.

    // Calculate coefficients for the terms in the differential operator
    BoutReal C_df_dx = coords->G1[l];
    BoutReal C_df_dz = coords->G3[l];
    if (issetD) {
      C_df_dx *= D[l];
      C_df_dz *= D[l];
    }
    if (issetC) {
      C_df_dx += (coords->g11[l] * dc_dx[l] + coords->g12[l] * dc_dy[l]
                  + coords->g13[l] * dc_dz[l])
                 / C1[l];
      C_df_dz += (coords->g13[l] * dc_dx[l] + coords->g23[l] * dc_dy[l]
                  + coords->g33[l] * dc_dz[l])
                 / C1[l];
    }
    if (issetE) {
      C_df_dx += Ex[l];
      C_df_dz += Ez[l];
    }

    BoutReal C_d2f_dx2 = coords->g11[l];
    BoutReal C_d2f_dy2 = (coords->g22[l] - 1.0 / coords->g_22[l]);
    BoutReal C_d2f_dz2 = coords->g33[l];
    if (issetD) {
      C_d2f_dx2 *= D[l];
      C_d2f_dy2 *= D[l];
      C_d2f_dz2 *= D[l];
    }

    BoutReal C_d2f_dxdz = 2 * coords->g13[l];
    if (issetD) {
      C_d2f_dxdz *= D[l];
    }

    // Adjust the coefficients to include finite-difference factors
    if (nonuniform) {
      C_df_dx += C_d2f_dx2 * coords->d1_dx[l];
    }
    C_df_dx /= 2 * coords->dx[l];
    C_df_dz /= 2 * coords->dz[l];

    C_d2f_dx2 /= SQ(coords->dx[l]);
    C_d2f_dy2 /= SQ(coords->dy[l]);
    C_d2f_dz2 /= SQ(coords->dz[l]);

    C_d2f_dxdz /= 4 * coords->dx[l] * coords->dz[l];

    operator3D(l, l) = -2 * (C_d2f_dx2 + C_d2f_dy2 + C_d2f_dz2) + A[l];
    operator3D(l, l.xp()) = C_df_dx + C_d2f_dx2;
    operator3D(l, l.xm()) = -C_df_dx + C_d2f_dx2;
    operator3D(l, l.zp()) = C_df_dz + C_d2f_dz2;
    operator3D(l, l.zm()) = -C_df_dz + C_d2f_dz2;
    operator3D(l, l.xp().zp()) = C_d2f_dxdz;
    operator3D(l, l.xp().zm()) = -C_d2f_dxdz;
    operator3D(l, l.xm().zp()) = -C_d2f_dxdz;
    operator3D(l, l.xm().zm()) = C_d2f_dxdz;
    // The values stored in the y-boundary are already interpolated
    // up/down, so we don't want the matrix to do any such
    // interpolation there.
    const int yup = (l.y() == localmesh->yend && upperY.intersects(l.x())) ? -1 : 0;
    const int ydown = (l.y() == localmesh->ystart && lowerY.intersects(l.x())) ? -1 : 0;
    operator3D.yup(yup)(l, l.yp()) = 0.0;
    operator3D.ydown(ydown)(l, l.ym()) = 0.0;
    operator3D.yup(yup)(l, l.xp().yp()) = 0.0;
    operator3D.ydown(ydown)(l, l.xp().ym()) = 0.0;
    operator3D.yup(yup)(l, l.xm().yp()) = 0.0;
    operator3D.ydown(ydown)(l, l.xm().ym()) = 0.0;
    operator3D.yup(yup)(l, l.yp().zp()) = 0.0;
    operator3D.yup(yup)(l, l.yp().zm()) = 0.0;
    operator3D.ydown(ydown)(l, l.ym().zp()) = 0.0;
    operator3D.ydown(ydown)(l, l.ym().zm()) = 0.0;
  }
  operator3D.partialAssemble();

  // Must add these (rather than assign) so that elements used in
  // interpolation don't overwrite each other.
  BOUT_FOR_SERIAL(l, indexer->getRegionNobndry()) {
    BoutReal C_df_dy = (coords->G2[l] - dJ_dy[l] / coords->J[l]);
    if (issetD) {
      C_df_dy *= D[l];
    }
    if (issetC) {
      C_df_dy +=
          (coords->g12[l] * dc_dx[l] + (coords->g22[l] - 1. / coords->g_22[l]) * dc_dy[l]
           + coords->g23[l] * dc_dz[l])
          / C1[l];
    }

    BoutReal C_d2f_dy2 = (coords->g22[l] - 1.0 / coords->g_22[l]);
    if (issetD) {
      C_d2f_dy2 *= D[l];
    }

    BoutReal C_d2f_dxdy = 2 * coords->g12[l];
    BoutReal C_d2f_dydz = 2 * coords->g23[l];
    if (issetD) {
      C_d2f_dxdy *= D[l];
      C_d2f_dydz *= D[l];
    }

    // Adjust the coefficients to include finite-difference factors
    if (nonuniform) {
      C_df_dy += C_d2f_dy2 * coords->d1_dy[l];
    }
    C_df_dy /= 2 * coords->dy[l];
    C_d2f_dy2 /= SQ(coords->dy[l]);
    C_d2f_dxdy /=
        4 * coords->dx[l]; // NOTE: This value is not completed here. It needs to
                           // be divide by dx(i +/- 1, j, k) when using to set a
                           // matrix element
    C_d2f_dydz /= 4 * coords->dy[l] * coords->dz[l];

    // The values stored in the y-boundary are already interpolated
    // up/down, so we don't want the matrix to do any such
    // interpolation there.
    const int yup = (l.y() == localmesh->yend && upperY.intersects(l.x())) ? -1 : 0;
    const int ydown = (l.y() == localmesh->ystart && lowerY.intersects(l.x())) ? -1 : 0;

    operator3D.yup(yup)(l, l.yp()) += C_df_dy + C_d2f_dy2;
    operator3D.ydown(ydown)(l, l.ym()) += -C_df_dy + C_d2f_dy2;
    operator3D.yup(yup)(l, l.xp().yp()) += C_d2f_dxdy / coords->dy[l.xp()];
    operator3D.ydown(ydown)(l, l.xp().ym()) += -C_d2f_dxdy / coords->dy[l.xp()];
    operator3D.yup(yup)(l, l.xm().yp()) += -C_d2f_dxdy / coords->dy[l.xm()];
    operator3D.ydown(ydown)(l, l.xm().ym()) += C_d2f_dxdy / coords->dy[l.xm()];
    operator3D.yup(yup)(l, l.yp().zp()) += C_d2f_dydz;
    operator3D.yup(yup)(l, l.yp().zm()) += -C_d2f_dydz;
    operator3D.ydown(ydown)(l, l.ym().zp()) += -C_d2f_dydz;
    operator3D.ydown(ydown)(l, l.ym().zm()) += C_d2f_dydz;
  }
  operator3D.assemble();
  MatSetBlockSize(*operator3D.get(), 1);

  // Declare KSP Context (abstract PETSc object that manages all Krylov methods)
  if (kspInitialised) {
    KSPDestroy(&ksp);
  }
  KSPCreate(BoutComm::get(), &ksp);
  kspInitialised = true;
#if PETSC_VERSION_GE(3, 5, 0)
  KSPSetOperators(ksp, *operator3D.get(), *operator3D.get());
#else
  KSPSetOperators(ksp, *operator3D.get(), *operator3D.get(), DIFFERENT_NONZERO_PATTERN);
#endif

  PC pc = nullptr;
  KSPGetPC(ksp, &pc);

  if (direct) {
    // Set the type of the preconditioner
    PCSetType(pc, PCLU);
    KSPSetType(ksp, KSPPREONLY);
#ifdef PETSC_HAVE_MUMPS
#if PETSC_VERSION_GE(3, 9, 0)
    PCFactorSetMatSolverType(pc, "mumps");
#else
    PCFactorSetMatSolverPackage(pc, "mumps");
#endif
#else
    // MUMPS not available, hope that PETSc has a working default option
#endif // PETSC_HAVE_MUMPS
  } else {
    KSPSetType(ksp, ksptype.c_str()); // Set the type of the solver

    if (ksptype == KSPRICHARDSON) {
      KSPRichardsonSetScale(ksp, richardson_damping_factor);
#ifdef KSPCHEBYSHEV
    } else if (ksptype == KSPCHEBYSHEV) {
      KSPChebyshevSetEigenvalues(ksp, chebyshev_max, chebyshev_min);
#endif
    } else if (ksptype == KSPGMRES) {
      KSPGMRESSetRestart(ksp, gmres_max_steps);
    }

    // Set the relative and absolute tolerances
    KSPSetTolerances(ksp, rtol, atol, dtol, maxits);

    // If the initial guess is not set to zero
    if (!isGlobalFlagSet(INVERT_START_NEW)) {
      KSPSetInitialGuessNonzero(ksp, (PetscBool) true);
    }

    // Set the relative and absolute tolerances
    PCSetType(pc, pctype.c_str());
#if PETSC_VERSION_LT(3, 18, 0)
    PCGAMGSetSymGraph(pc, PETSC_TRUE);
#endif
  }
  lib.setOptionsFromInputFile(ksp);

  updateRequired = false;
}

OperatorStencil<Ind3D> LaplacePetsc3dAmg::getStencil(Mesh* localmesh,
                                                     const RangeIterator& lowerYBound,
                                                     const RangeIterator& upperYBound) {
  OperatorStencil<Ind3D> stencil;

  // Get the pattern used for interpolation. This is assumed to be the
  // same across the whole grid.
  const auto positions_weights =
      localmesh->getCoordinates()->getParallelTransform().getWeightsForYDownApproximation(
          localmesh->xstart, localmesh->ystart + 1, localmesh->zstart);
  std::vector<OffsetInd3D> interpPattern;
  std::transform(
      positions_weights.begin(), positions_weights.end(),
      std::back_inserter(interpPattern),
      [localmesh](ParallelTransform::PositionsAndWeights position) -> OffsetInd3D {
        return {localmesh->xstart - position.i, localmesh->ystart - position.j,
                ((localmesh->LocalNz - position.k) < position.k)
                    ? position.k - localmesh->LocalNz
                    : position.k};
      });

  const OffsetInd3D zero;

  // Add interior cells
  const std::vector<OffsetInd3D> interpolatedUpElements = {
      zero.yp(), zero.xp().yp(), zero.xm().yp(), zero.yp().zp(), zero.yp().zm()};
  const std::vector<OffsetInd3D> interpolatedDownElements = {
      zero.ym(), zero.xp().ym(), zero.xm().ym(), zero.ym().zp(), zero.ym().zm()};
  std::set<OffsetInd3D> interiorStencil = {
      zero,           zero.xp(),      zero.xm(),      zero.zp(),     zero.zm(),
      zero.xp().zp(), zero.xp().zm(), zero.xm().zp(), zero.xm().zm()};
  std::set<OffsetInd3D> lowerEdgeStencil = interiorStencil;
  std::set<OffsetInd3D> upperEdgeStencil = interiorStencil;

  for (const auto& i : interpolatedDownElements) {
    for (auto& j : interpPattern) {
      interiorStencil.insert(i + j);
      upperEdgeStencil.insert(i + j);
    }
    lowerEdgeStencil.insert(i);
  }
  for (const auto& i : interpolatedUpElements) {
    for (auto& j : interpPattern) {
      interiorStencil.insert(i + j);
      lowerEdgeStencil.insert(i + j);
    }
    upperEdgeStencil.insert(i);
  }
  const std::vector<OffsetInd3D> interiorStencilVector(interiorStencil.begin(),
                                                       interiorStencil.end());
  const std::vector<OffsetInd3D> lowerEdgeStencilVector(lowerEdgeStencil.begin(),
                                                        lowerEdgeStencil.end());
  const std::vector<OffsetInd3D> upperEdgeStencilVector(upperEdgeStencil.begin(),
                                                        upperEdgeStencil.end());

  // If there is a lower y-boundary then create a part of the stencil
  // for cells immediately adjacent to it.
  if (lowerYBound.max() - lowerYBound.min() > 0) {
    stencil.add(
        [index = localmesh->ystart, lowerYBound](Ind3D ind) -> bool {
          return index == ind.y() && lowerYBound.intersects(ind.x());
        },
        lowerEdgeStencilVector);
  }

  // If there is an upper y-boundary then create a part of the stencil
  // for cells immediately adjacent to it.
  if (upperYBound.max() - upperYBound.min() > 0) {
    stencil.add(
        [index = localmesh->yend, upperYBound](Ind3D ind) -> bool {
          return index == ind.y() && upperYBound.intersects(ind.x());
        },
        upperEdgeStencilVector);
  }

  // Create a part of the stencil for the interior cells. Although the
  // test here would also pass for the edge-cells immediately adjacent
  // to upper/lower y-boundaries, because those tests are run first
  // the cells will be assigned to those regions.
  stencil.add(
      [localmesh](Ind3D ind) -> bool {
        return (localmesh->xstart <= ind.x() && ind.x() <= localmesh->xend
                && localmesh->ystart <= ind.y() && ind.y() <= localmesh->yend
                && localmesh->zstart <= ind.z() && ind.z() <= localmesh->zend);
      },
      interiorStencilVector);

  // Add Y boundaries before X boundaries so corners are assigned to
  // the former.

  // Note that the tests for whether a cell is in the boundary
  // actually are not exclusive enough. They really only correspond to
  // whether the cell is in the guard regions. However, the tests are
  // much simpler to implement this way and other checks will prevent
  // guard-cells which are not boundaries from having memory
  // pre-allocated.

  // Add lower Y boundary.
  stencil.add(
      [index = localmesh->ystart - 1](Ind3D ind) -> bool { return ind.y() == index; },
      {zero, zero.yp()});
  // Add upper Y boundary
  stencil.add(
      [index = localmesh->yend + 1](Ind3D ind) -> bool { return ind.y() == index; },
      {zero, zero.ym()});
  // Add inner X boundary
  if (localmesh->firstX()) {
    stencil.add(
        [index = localmesh->xstart - 1](Ind3D ind) -> bool { return ind.x() == index; },
        {zero, zero.xp()});
  }
  // Add outer X boundary
  if (localmesh->lastX()) {
    stencil.add(
        [index = localmesh->xend + 1](Ind3D ind) -> bool { return ind.x() == index; },
        {zero, zero.xm()});
  }

  return stencil;
}

#endif // BOUT_HAS_PETSC
