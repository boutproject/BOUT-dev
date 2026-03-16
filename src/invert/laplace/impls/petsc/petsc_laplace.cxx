/**************************************************************************
 * Perpendicular Laplacian inversion.
 *                           Using PETSc Solvers
 *
 **************************************************************************
 * Copyright 2013 - 2025 BOUT++ developers
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "petsc_laplace.hxx"

#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/field2d.hxx>
#include <bout/globalindexer.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/operatorstencil.hxx>
#include <bout/output.hxx>
#include <bout/petsc_interface.hxx>
#include <bout/petsclib.hxx>
#include <bout/region.hxx>
#include <bout/sys/timer.hxx>
#include <bout/utils.hxx>

#include <cmath>
#include <memory>
#include <set>
#include <vector>

namespace {
PetscErrorCode laplacePCapply(PC pc, Vec x, Vec y) {
  PetscFunctionBegin; // NOLINT

  LaplacePetsc* laplace = nullptr;
  CHKERRQ(PCShellGetContext(pc, reinterpret_cast<void**>(&laplace))); // NOLINT

  PetscFunctionReturn(laplace->precon(x, y)); // NOLINT
}

auto set_stencil(const Mesh& localmesh, bool fourth_order) {
  OperatorStencil<IndPerp> stencil;
  IndexOffset<IndPerp> zero;
  // Start with a square stencil 1-point wide
  std::set offsets = {
      // clang-format off
      zero.xm().zp(), zero.zp(), zero.xp().zp(),
      zero.xm(),      zero,      zero.xp(),
      zero.xm().zm(), zero.zm(), zero.xp().zm(),
      // clang-format on
  };

  if (fourth_order) {
    // Add a square stencil 2-points wide
    offsets.insert({
        // clang-format off
        zero.xm(2).zp(2), zero.xm().zp(2), zero.zp(2), zero.xp().zp(2), zero.xp(2).zp(2),
        zero.xm(2).zp(),                                                zero.xp(2).zp(),
        zero.xm(2),                                                     zero.xp(2),
        zero.xm(2).zm(),                                                zero.xp(2).zm(),
        zero.xm(2).zm(2), zero.xm().zm(2), zero.zm(2), zero.xp().zm(2), zero.xp(2).zm(2),
        // clang-format on
    });
  }

  const std::vector offsetsVec(offsets.begin(), offsets.end());
  stencil.add(
      [&localmesh](IndPerp ind) -> bool {
        return (localmesh.xstart <= ind.x() && ind.x() <= localmesh.xend
                and (localmesh.zstart <= ind.z() && ind.z() <= localmesh.zend));
      },
      offsetsVec);

  // Add inner X boundary
  if (localmesh.firstX()) {
    const auto first_boundary = localmesh.xstart - 1;
    const auto second_boundary = localmesh.xstart - 2;

    if (fourth_order) {
      stencil.add(
          [first_boundary, second_boundary](IndPerp ind) -> bool {
            const auto x = ind.x();
            return x == first_boundary or x == second_boundary;
          },
          {zero, zero.xp(1), zero.xp(2), zero.xp(3), zero.xp(4)});
    } else {
      stencil.add(
          [first_boundary](IndPerp ind) -> bool { return ind.x() == first_boundary; },
          {zero, zero.xp()});
    }
  }
  // Add outer X boundary
  if (localmesh.lastX()) {
    const auto first_boundary = localmesh.xend + 1;
    const auto second_boundary = localmesh.xend + 2;

    if (fourth_order) {
      stencil.add(
          [first_boundary, second_boundary](IndPerp ind) -> bool {
            const auto x = ind.x();
            return x == first_boundary or x == second_boundary;
          },
          {zero, zero.xm(1), zero.xm(2), zero.xm(3), zero.xm(4)});
    } else {
      stencil.add(
          [first_boundary](IndPerp ind) -> bool { return ind.x() == first_boundary; },
          {zero, zero.xm()});
    }
  }

  stencil.add([]([[maybe_unused]] IndPerp ind) -> bool { return true; }, {zero});
  return stencil;
}
} // namespace

LaplacePetsc::LaplacePetsc(Options* opt, const CELL_LOC loc, Mesh* mesh_in,
                           [[maybe_unused]] Solver* solver)
    : Laplacian(opt, loc, mesh_in), A(0.0, mesh_in), C1(1.0, mesh_in), C2(1.0, mesh_in),
      D(1.0, mesh_in), Ex(0.0, mesh_in), Ez(0.0, mesh_in), issetD(false), issetC(false),
      issetE(false), comm(localmesh->getXZcomm()),
      opts(opt == nullptr ? &(Options::root()["laplace"]) : opt),
      // WARNING: only a few of these options actually make sense: see the
      // PETSc documentation to work out which they are (possibly
      // pbjacobi, sor might be useful choices?)
      ksptype((*opts)["ksptype"].doc("KSP solver type").withDefault(KSPGMRES)),
      pctype((*opts)["pctype"]
                 .doc("Preconditioner type. See the PETSc documentation for options")
                 .withDefault("none")),
      richardson_damping_factor((*opts)["richardson_damping_factor"].withDefault(1.0)),
      chebyshev_max((*opts)["chebyshev_max"].withDefault(100)),
      chebyshev_min((*opts)["chebyshev_min"].withDefault(0.01)),
      gmres_max_steps((*opts)["gmres_max_steps"].withDefault(30)),
      rtol((*opts)["rtol"].doc("Relative tolerance for KSP solver").withDefault(1e-5)),
      atol((*opts)["atol"].doc("Absolute tolerance for KSP solver").withDefault(1e-50)),
      dtol((*opts)["dtol"].doc("Divergence tolerance for KSP solver").withDefault(1e5)),
      maxits(
          (*opts)["maxits"].doc("Maximum number of KSP iterations").withDefault(100000)),
      direct((*opts)["direct"].doc("Use direct (LU) solver?").withDefault(false)),
      fourth_order(
          (*opts)["fourth_order"].doc("Use fourth order stencil").withDefault(false)),
      indexer(std::make_shared<GlobalIndexer<FieldPerp>>(
          localmesh, set_stencil(*localmesh, fourth_order))),
      operator2D(indexer), lib(opts) {

  A.setLocation(location);
  C1.setLocation(location);
  C2.setLocation(location);
  D.setLocation(location);
  Ex.setLocation(location);
  Ez.setLocation(location);

#if CHECK > 0
  // Checking flags are set to something which is not implemented
  checkFlags();

  if (localmesh->periodicX) {
    throw BoutException("LaplacePetsc does not work with periodicity in the x direction "
                        "(localmesh->PeriodicX == true). Change boundary conditions or "
                        "use serial-tri or cyclic solver instead");
  }
#endif

  // Let "user" be a synonym for "shell"
  if (pctype == "user") {
    pctype = PCSHELL;
  }

  if (direct) {
    output.write("\nUsing LU decompostion for direct solution of system\n");
  }

  if (pctype == PCSHELL) {
    rightprec = (*opts)["rightprec"].doc("Right preconditioning?").withDefault(true);

    // Options for preconditioner are in a subsection
    pcsolve = Laplacian::create(opts->getSection("precon"));
  }
}

LaplacePetsc::~LaplacePetsc() {
  if (ksp_initialised) {
    KSPDestroy(&ksp);
  }
}

FieldPerp LaplacePetsc::solve(const FieldPerp& b) { return solve(b, b); }

/*!
 * Solves Ax=b for x given a b and an initial guess for x (x0)
 *
 * This function will:
 *      1. Set the matrix element of the matrix A, used to solve Ax=b
 *         (this includes setting the values for the bounary condition)
 *      2. Solve the matrix Ax = b
 *
 * \param[in] b     The RHS of the equation Ax=b.
 *                  This is an y-slice of the original field. The field wil be
 *                  flattened to an 1D array in order to write the equation on
 *                  the form Ax=b
 * \param[in] x0    The initial guess for the solver.
 *                  May also contain the  boundary condition if flag 32 - INVERT_SET is set
 *
 * \returns sol     The solution x of the problem Ax=b.
 */
FieldPerp LaplacePetsc::solve(const FieldPerp& b, const FieldPerp& x0) {

  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

#if CHECK > 0
  checkFlags();
#endif

  // Set member variable so that we can pass through to shell preconditioner if
  // required
  yindex = b.getIndex();
  {
    const Timer timer("petscsetup");

    const bool inner_X_neumann = isInnerBoundaryFlagSet(INVERT_AC_GRAD);
    const bool outer_X_neumann = isOuterBoundaryFlagSet(INVERT_AC_GRAD);

    // Set the operator matrix
    if (fourth_order) {
      setFourthOrderMatrix(yindex, inner_X_neumann, outer_X_neumann);
    } else {
      setSecondOrderMatrix(yindex, inner_X_neumann, outer_X_neumann);
    }

    operator2D.assemble();
    MatSetBlockSize(*operator2D.get(), 1);

    // Declare KSP Context (abstract PETSc object that manages all Krylov methods)
    if (ksp_initialised) {
      KSPDestroy(&ksp);
    }
    KSPCreate(comm, &ksp);

    // Configure Linear Solver
#if PETSC_VERSION_GE(3, 5, 0)
    KSPSetOperators(ksp, *operator2D.get(), *operator2D.get());
#else
    KSPSetOperators(ksp, *operator2D.get(), *operator2D.get(), DIFFERENT_NONZERO_PATTERN);
#endif
    PC pc = nullptr; // The preconditioner option

    if (direct) { // If a direct solver has been chosen
      // Get the preconditioner
      KSPGetPC(ksp, &pc);
      // Set the preconditioner
      PCSetType(pc, PCLU);
      // Set the solver type
#if PETSC_VERSION_GE(3, 9, 0)
      PCFactorSetMatSolverType(pc, "mumps");
#else
      PCFactorSetMatSolverPackage(pc, "mumps");
#endif
    } else {                            // If a iterative solver has been chosen
      KSPSetType(ksp, ksptype.c_str()); // Set the type of the solver

      if (ksptype == KSPRICHARDSON) {
        KSPRichardsonSetScale(ksp, richardson_damping_factor);
      }
#ifdef KSPCHEBYSHEV
      else if (ksptype == KSPCHEBYSHEV) {
        KSPChebyshevSetEigenvalues(ksp, chebyshev_max, chebyshev_min);
      }
#endif
      else if (ksptype == KSPGMRES) {
        KSPGMRESSetRestart(ksp, gmres_max_steps);
      }

      // Set the relative and absolute tolerances
      KSPSetTolerances(ksp, rtol, atol, dtol, maxits);

      // If the initial guess is not set to zero
      if (!isGlobalFlagSet(INVERT_START_NEW)) {
        KSPSetInitialGuessNonzero(ksp, static_cast<PetscBool>(true));
      }

      // Get the preconditioner
      KSPGetPC(ksp, &pc);

      // Set the type of the preconditioner
      PCSetType(pc, pctype.c_str());

      // If pctype = user in BOUT.inp, it will be translated to PCSHELL upon
      // construction of the object
      if (pctype == PCSHELL) {
        // User-supplied preconditioner function
        PCShellSetApply(pc, laplacePCapply);
        PCShellSetContext(pc, this);
        if (rightprec) {
          KSPSetPCSide(ksp, PC_RIGHT); // Right preconditioning
        } else {
          KSPSetPCSide(ksp, PC_LEFT); // Left preconditioning
        }
      }

      lib.setOptionsFromInputFile(ksp);
    }
  }

  PetscVector<FieldPerp> rhs(b, indexer);
  PetscVector<FieldPerp> guess(x0, indexer);

  // Set boundary conditions
  if (!isInnerBoundaryFlagSet(INVERT_RHS)) {
    BOUT_FOR_SERIAL(index, indexer->getRegionInnerX()) {
      rhs(index) = isInnerBoundaryFlagSet(INVERT_SET) ? x0[index] : 0.0;
    }
  }
  if (!isOuterBoundaryFlagSet(INVERT_RHS)) {
    BOUT_FOR_SERIAL(index, indexer->getRegionOuterX()) {
      rhs(index) = isInnerBoundaryFlagSet(INVERT_SET) ? x0[index] : 0.0;
    }
  }

  rhs.assemble();
  guess.assemble();

  // Call the actual solver
  {
    const Timer timer("petscsolve");
    KSPSolve(ksp, *rhs.get(), *guess.get());
  }

  KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
  KSPGetConvergedReason(ksp, &reason);
  if (reason == -3) { // Too many iterations, might be fixed by taking smaller timestep
    throw BoutIterationFail("petsc_laplace: too many iterations");
  }
  if (reason <= 0) {
    throw BoutException(
        "petsc_laplace: inversion failed to converge. KSPConvergedReason: {} ({})",
        KSPConvergedReasons[reason], static_cast<int>(reason));
  }

  auto sol = guess.toField();
  sol.setIndex(yindex);
  checkData(sol);

  // Return the solution
  return sol;
}

/*!
 * Set the matrix components of A in Ax=b, solving
 * D*Laplace_perp(x) + (1/C1)Grad_perp(C2)*Grad_perp(x) + Ax = B
 *
 * \note "A" in the equation above is not added here.
 * For calculations of the coefficients, please refer to the user manual.
 *
 * \param[in] x The current x index
 * \param[in] y The current y index
 * \param[in] z The current y index
 * \param[in] coef1  Placeholder for convenient variable used to set matrix
 *                   (see manual for details)
 * \param[in] coef2  Convenient variable used to set matrix
 *                   (see manual for details)
 * \param[in] coef3  Placeholder for convenient variable used to set matrix
 *                   (see manual for details)
 * \param[in] coef4  Placeholder for convenient variable used to set matrix
 *                   (see manual for details)
 * \param[in] coef5  Placeholder for convenient variable used to set matrix
 *                   (see manual for details)
 *
 * \param[out] coef1    Convenient variable used to set matrix
 *                      (see manual for details)
 * \param[out] coef2    Convenient variable used to set matrix
 *                      (see manual for details)
 * \param[out] coef3    Convenient variable used to set matrix
 *                      (see manual for details)
 * \param[out] coef4    Convenient variable used to set matrix
 *                      (see manual for details)
 * \param[out] coef5    Convenient variable used to set matrix
 *                      (see manual for details)
 */
LaplacePetsc::CoeffsA LaplacePetsc::Coeffs(Ind3D i) {
  const auto x = i.x();

  BoutReal coef1 = coords->g11[i];      // X 2nd derivative coefficient
  BoutReal coef2 = coords->g33[i];      // Z 2nd derivative coefficient
  BoutReal coef3 = 2. * coords->g13[i]; // X-Z mixed derivative coefficient

  BoutReal coef4 = 0.0;
  BoutReal coef5 = 0.0;
  // If global flag all_terms are set (true by default)
  if (all_terms) {
    coef4 = coords->G1[i]; // X 1st derivative
    coef5 = coords->G3[i]; // Z 1st derivative

    ASSERT3(std::isfinite(coef4));
    ASSERT3(std::isfinite(coef5));
  }

  if (nonuniform) {
    // non-uniform mesh correction
    if ((x != 0) && (x != (localmesh->LocalNx - 1))) {
      coef4 -= 0.5 * ((coords->dx[i.xp()] - coords->dx[i.xm()]) / SQ(coords->dx[i]))
               * coef1; // BOUT-06 term
    }
  }

  if (localmesh->IncIntShear) {
    // d2dz2 term
    coef2 += coords->g11[i] * coords->IntShiftTorsion[i] * coords->IntShiftTorsion[i];
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }

  if (issetD) {
    coef1 *= D[i];
    coef2 *= D[i];
    coef3 *= D[i];
    coef4 *= D[i];
    coef5 *= D[i];
  }

  // A second/fourth order derivative term
  if (issetC) {
    if ((x > 1) && (x < (localmesh->LocalNx - 2))) {
      BoutReal ddx_C = BoutNaN;
      BoutReal ddz_C = BoutNaN;

      if (fourth_order) {
        // Fourth order discretization of C in x
        ddx_C = (-C2[i.xpp()] + (8. * C2[i.xp()]) - (8. * C2[i.xm()]) + C2[i.xmm()])
                / (12. * coords->dx[i] * (C1[i]));
        // Fourth order discretization of C in z
        ddz_C = (-C2[i.zpp()] + (8. * C2[i.zp()]) - (8. * C2[i.zm()]) + C2[i.zmm()])
                / (12. * coords->dz[i] * (C1[i]));
      } else {
        // Second order discretization of C in x
        ddx_C = (C2[i.xp()] - C2[i.xm()]) / (2. * coords->dx[i] * (C1[i]));
        // Second order discretization of C in z
        ddz_C = (C2[i.zp()] - C2[i.zm()]) / (2. * coords->dz[i] * (C1[i]));
      }

      coef4 += (coords->g11[i] * ddx_C) + (coords->g13[i] * ddz_C);
      coef5 += (coords->g13[i] * ddx_C) + (coords->g33[i] * ddz_C);
    }
  }

  /* Ex and Ez
   * Additional 1st derivative terms to allow for solution field to be
   * components of a vector
   *
   * NB multiply by D or Grad_perp(C)/C as appropriate before passing to
   * setCoefEx()/setCoefEz() because (in principle) both are needed and we
   * don't know how to split them up here
   */
  if (issetE) {
    // These coefficients are 0 by default
    coef4 += Ex[i];
    coef5 += Ez[i];
  }

  return {coef1, coef2, coef3, coef4, coef5};
}

void LaplacePetsc::setSecondOrderMatrix(int y, bool inner_X_neumann,
                                        bool outer_X_neumann) {
  // Set the boundaries
  if (inner_X_neumann) {
    const auto dx = sliceXZ(coords->dx, y);
    const auto g11 = sliceXZ(coords->g11, y);

    BOUT_FOR_SERIAL(i, indexer->getRegionInnerX()) {
      const auto factor = 1. / dx[i] / std::sqrt(g11[i]);
      operator2D(i, i) = -factor;
      operator2D(i, i.xp()) = factor;
    }
  } else {
    BOUT_FOR_SERIAL(i, indexer->getRegionInnerX()) {
      operator2D(i, i) = 0.5;
      operator2D(i, i.xp()) = 0.5;
    }
  }
  if (outer_X_neumann) {
    const auto dx = sliceXZ(coords->dx, y);
    const auto g11 = sliceXZ(coords->g11, y);

    BOUT_FOR_SERIAL(i, indexer->getRegionOuterX()) {
      const auto factor = 1. / dx[i] / std::sqrt(g11[i]);
      operator2D(i, i) = factor;
      operator2D(i, i.xm()) = -factor;
    }
  } else {
    BOUT_FOR_SERIAL(i, indexer->getRegionOuterX()) {
      operator2D(i, i) = 0.5;
      operator2D(i, i.xm()) = 0.5;
    }
  }

  // Set the interior region
  BOUT_FOR_SERIAL(l, indexer->getRegionNobndry()) {
    const auto i = localmesh->indPerpto3D(l, y);

    // NOTE: Only A0 is the A from setCoefA ()
    const BoutReal A0 = A[i];

    ASSERT3(std::isfinite(A0));

    // Set the matrix coefficients
    const auto [A1, A2, A3, A4, A5] = Coeffs(i);

    ASSERT3(std::isfinite(A1));
    ASSERT3(std::isfinite(A2));
    ASSERT3(std::isfinite(A3));
    ASSERT3(std::isfinite(A4));
    ASSERT3(std::isfinite(A5));

    const BoutReal dx = coords->dx[i];
    const BoutReal dx2 = SQ(dx);
    const BoutReal dz = coords->dz[i];
    const BoutReal dz2 = SQ(dz);
    const BoutReal dxdz = dx * dz;
    operator2D(l, l) = A0 - (2.0 * ((A1 / dx2) + (A2 / dz2)));
    operator2D(l, l.xm().zm()) = A3 / (4.0 * dxdz);
    operator2D(l, l.xm()) = (A1 / dx2) - (A4 / (2.0 * dx));
    operator2D(l, l.xm().zp()) = -1.0 * A3 / (4.0 * dxdz);
    operator2D(l, l.zm()) = (A2 / dz2) - (A5 / (2.0 * dz));
    operator2D(l, l.zp()) = (A2 / dz2) + (A5 / (2.0 * dz));
    operator2D(l, l.xp().zm()) = -1.0 * A3 / (4.0 * dxdz);
    operator2D(l, l.xp()) = (A1 / dx2) + (A4 / (2.0 * dx));
    operator2D(l, l.xp().zp()) = A3 / (4.0 * dxdz);
  }
}

void LaplacePetsc::setFourthOrderMatrix(int y, bool inner_X_neumann,
                                        bool outer_X_neumann) {

  // Set boundaries
  if (inner_X_neumann) {
    const auto dx = sliceXZ(coords->dx, y);
    const auto g11 = sliceXZ(coords->g11, y);

    BOUT_FOR_SERIAL(i, indexer->getRegionInnerX()) {
      const auto factor = 1. / dx[i] / std::sqrt(g11[i]);
      operator2D(i, i) = (-25.0 / 12.0) * factor;
      operator2D(i, i.xp(1)) = 4.0 * factor;
      operator2D(i, i.xp(2)) = -3.0 * factor;
      operator2D(i, i.xp(3)) = (4.0 / 3.0) * factor;
      operator2D(i, i.xp(4)) = (-1.0 / 4.0) * factor;
    }
  } else {
    BOUT_FOR_SERIAL(i, indexer->getRegionInnerX()) {
      operator2D(i, i) = 1.0;
      operator2D(i, i.xp(1)) = 0.0;
      operator2D(i, i.xp(2)) = 0.0;
      operator2D(i, i.xp(3)) = 0.0;
      operator2D(i, i.xp(4)) = 0.0;
    }
  }

  if (outer_X_neumann) {
    const auto dx = sliceXZ(coords->dx, y);
    const auto g11 = sliceXZ(coords->g11, y);

    BOUT_FOR_SERIAL(i, indexer->getRegionOuterX()) {
      const auto factor = 1. / dx[i] / std::sqrt(g11[i]);
      operator2D(i, i) = (25.0 / 12.0) * factor;
      operator2D(i, i.xm(1)) = -4.0 * factor;
      operator2D(i, i.xm(2)) = 3.0 * factor;
      operator2D(i, i.xm(3)) = (-4.0 / 3.0) * factor;
      operator2D(i, i.xm(4)) = (1.0 / 4.0) * factor;
    }
  } else {
    BOUT_FOR_SERIAL(i, indexer->getRegionOuterX()) {
      operator2D(i, i) = 1.0;
      operator2D(i, i.xm(1)) = 0.0;
      operator2D(i, i.xm(2)) = 0.0;
      operator2D(i, i.xm(3)) = 0.0;
      operator2D(i, i.xm(4)) = 0.0;
    }
  }

  // Set interior region
  BOUT_FOR_SERIAL(l, indexer->getRegionNobndry()) {
    const auto i = localmesh->indPerpto3D(l, y);

    // NOTE: Only A0 is the A from setCoefA ()
    const BoutReal A0 = A[i];

    ASSERT3(std::isfinite(A0));

    // Set the matrix coefficients
    const auto [A1, A2, A3, A4, A5] = Coeffs(i);

    ASSERT3(std::isfinite(A1));
    ASSERT3(std::isfinite(A2));
    ASSERT3(std::isfinite(A3));
    ASSERT3(std::isfinite(A4));
    ASSERT3(std::isfinite(A5));

    const BoutReal dx = coords->dx[i];
    const BoutReal dx2 = SQ(dx);
    const BoutReal dz = coords->dz[i];
    const BoutReal dz2 = SQ(dz);
    const BoutReal dxdz = dx * dz;

    operator2D(l, l) = A0 - ((5.0 / 2.0) * ((A1 / dx2) + (A2 / dz2)));
    operator2D(l, l.xmm().zmm()) = A3 / (144.0 * dxdz);
    operator2D(l, l.xmm().zm()) = -1.0 * A3 / (18.0 * dxdz);
    operator2D(l, l.xmm()) = (1.0 / 12.0) * ((-1.0 * A1 / dx2) + (A4 / dx));
    operator2D(l, l.xmm().zp()) = A3 / (18.0 * dxdz);
    operator2D(l, l.xmm().zpp()) = -1.0 * A3 / (144.0 * dxdz);
    operator2D(l, l.xm().zmm()) = -1.0 * A3 / (18.0 * dxdz);
    operator2D(l, l.xm().zm()) = 4.0 * A3 / (9.0 * dxdz);
    operator2D(l, l.xm()) = (4.0 * A1 / (3.0 * dx2)) - (2.0 * A4 / (3.0 * dx));
    operator2D(l, l.xm().zp()) = -4.0 * A3 / (9.0 * dxdz);
    operator2D(l, l.xm().zpp()) = A3 / (18.0 * dxdz);
    operator2D(l, l.zmm()) = (1.0 / 12.0) * ((-1.0 * A2 / dz2) + (A5 / dz));
    operator2D(l, l.zm()) = (4.0 * A2 / (3.0 * dz2)) - (2.0 * A5 / (3.0 * dz));
    operator2D(l, l.zp()) = (4.0 * A2 / (3.0 * dz2)) + (2.0 * A5 / (3.0 * dz));
    operator2D(l, l.zpp()) = (-1.0 / 12.0) * ((A2 / dz2) + (A5 / dz));
    operator2D(l, l.xp().zmm()) = A3 / (18.0 * dxdz);
    operator2D(l, l.xp().zm()) = -4.0 * A3 / (9.0 * dxdz);
    operator2D(l, l.xp()) = (4.0 * A1 / (3.0 * dx2)) + (2.0 * A4 / (3.0 * dx));
    operator2D(l, l.xp().zp()) = 4.0 * A3 / (9.0 * dxdz);
    operator2D(l, l.xp().zpp()) = -1.0 * A3 / (18.0 * dxdz);
    operator2D(l, l.xpp().zmm()) = -1.0 * A3 / (144.0 * dxdz);
    operator2D(l, l.xpp().zm()) = A3 / (18.0 * dxdz);
    operator2D(l, l.xpp()) = (-1.0 / 12.0) * ((A1 / dx2) + (A4 / dx));
    operator2D(l, l.xpp().zp()) = -1.0 * A3 / (18.0 * dxdz);
    operator2D(l, l.xpp().zpp()) = A3 / (144.0 * dxdz);
  }
}

/// Preconditioner function
int LaplacePetsc::precon(Vec x, Vec y) {
  FieldPerp xfield(indexer->getMesh(), location, yindex);
  xfield = 0.0;

  BOUT_FOR_SERIAL(i, indexer->getRegionAll()) {
    const auto ind = indexer->getGlobal(i);
    PetscScalar val = BoutNaN;
    VecGetValues(x, 1, &ind, &val);
    xfield[i] = val;
  }

  // Call the preconditioner solver
  const FieldPerp yfield = pcsolve->solve(xfield);

  VecCopy(*PetscVector{yfield, indexer}.get(), y);

  return 0;
}

void LaplacePetsc::checkFlags() {
  if (isGlobalFlagSet(~implemented_flags)) {
    if (isGlobalFlagSet(INVERT_4TH_ORDER)) {
      output_error.write(
          "For PETSc based Laplacian inverter, use 'fourth_order=true' instead of "
          "setting INVERT_4TH_ORDER flag\n");
    }
    throw BoutException("Attempted to set Laplacian inversion flag that is not "
                        "implemented in petsc_laplace.cxx");
  }
  if (isInnerBoundaryFlagSet(~implemented_boundary_flags)) {
    throw BoutException("Attempted to set Laplacian inversion boundary flag that is not "
                        "implemented in petsc_laplace.cxx");
  }
  if (isOuterBoundaryFlagSet(~implemented_boundary_flags)) {
    throw BoutException("Attempted to set Laplacian inversion boundary flag that is not "
                        "implemented in petsc_laplace.cxx");
  }
}

#endif // BOUT_HAS_PETSC_3_3
