/**************************************************************************
 * Interface to PETSc solver
 *
 **************************************************************************
 * Copyright 2010 - 2025 BOUT++ contributors
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

#include "petsc.hxx"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/globals.hxx>
#include <bout/msg_stack.hxx>
#include <bout/output.hxx>
#include <bout/petsc_interface.hxx>
#include <bout/solver.hxx>
#include <bout/unused.hxx>
#include <bout/utils.hxx>

#include <petsc.h>
#include <petscerror.h>
#include <petscistypes.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscsnes.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petscts.h>
#include <petscvec.h>

#ifndef PETSC_UNLIMITED
// Introduced in PETSc 3.22
#define PETSC_UNLIMITED (-3)
#endif

class ColoringStencil {
private:
  bool static isInSquare(int const i, int const j, int const n_square) {
    return std::abs(i) <= n_square && std::abs(j) <= n_square;
  }
  bool static isInCross(int const i, int const j, int const n_cross) {
    if (i == 0) {
      return std::abs(j) <= n_cross;
    }
    if (j == 0) {
      return std::abs(i) <= n_cross;
    }
    return false;
  }
  bool static isInTaxi(int const i, int const j, int const n_taxi) {
    return std::abs(i) + std::abs(j) <= n_taxi;
  }

public:
  auto static getOffsets(int n_square, int n_taxi, int n_cross) {
    ASSERT2(n_square >= 0 && n_cross >= 0 && n_taxi >= 0
            && n_square + n_cross + n_taxi > 0);
    auto inside = [&](int i, int j) {
      return isInSquare(i, j, n_square) || isInTaxi(i, j, n_taxi)
             || isInCross(i, j, n_cross);
    };
    std::vector<std::pair<int, int>> xy_offsets;
    auto loop_bound = std::max({n_square, n_taxi, n_cross});
    for (int i = -loop_bound; i <= loop_bound; ++i) {
      for (int j = -loop_bound; j <= loop_bound; ++j) {
        if (inside(i, j)) {
          xy_offsets.emplace_back(i, j);
        }
      }
    }
    return xy_offsets;
  }
};

namespace {
// PETSc callback function for matrix-free preconditioner
PetscErrorCode snesPCapply(PC pc, Vec x, Vec y) {
  // Get the context
  void* ctx = nullptr;
  PetscCall(PCShellGetContext(pc, &ctx));
  // Run the preconditioner
  PetscFunctionReturn(static_cast<PetscSolver*>(ctx)->pre(x, y));
}

// PETSc callback function, that evaluates the nonlinear
// function being integrated by TS.
PetscErrorCode solver_rhs(TS UNUSED(ts), BoutReal t, Vec globalin, Vec globalout,
                          void* f_data) {
  PetscFunctionBegin;
  auto* s = static_cast<PetscSolver*>(f_data);
  s->rhs(t, globalin, globalout, false);
  PetscFunctionReturn(0);
}

// Form function for use with SUNDIALS.
// This is needed because SNESTSFormFunction is not available
// PETSc error: No method snesfunction for TS of type sundials
PetscErrorCode solver_form_function(void* UNUSED(dummy), Vec U, Vec F, void* f_data) {
  PetscFunctionBegin;
  auto* s = static_cast<PetscSolver*>(f_data);
  s->formFunction(U, F);
  PetscFunctionReturn(0);
}
} // namespace

// Compute IJacobian = dF/dU + a dF/dUdot
// This is a dummy matrix that saves the shift.
// The shift is later used in the matrix-free preconditioner
PetscErrorCode solver_ijacobian(TS, BoutReal, Vec, Vec, PetscReal shift, Mat J, Mat Jpre,
                                void* ctx) {
  auto* solver = static_cast<PetscSolver*>(ctx);
  solver->shift = shift;

  PetscCall(MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY));
  if (J != Jpre) {
    PetscCall(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));
  }

  PetscFunctionReturn(0);
}

PetscErrorCode solver_ijacobian_color(TS ts, PetscReal t, Vec U, Vec Udot,
                                      PetscReal shift, Mat J, Mat B, void* ctx) {
  auto* solver = static_cast<PetscSolver*>(ctx);
  solver->shift = shift;

  return TSComputeIJacobianDefaultColor(ts, t, U, Udot, shift, J, B, nullptr);
}

// This function is called by the TS object every internal timestep
// It is responsible for triggering interpolation and output at
// the output time interval.
PetscErrorCode PetscMonitor(TS ts, PetscInt UNUSED(step), PetscReal t, Vec X, void* ctx) {
  PetscFunctionBegin;

  auto* s = static_cast<PetscSolver*>(ctx);
  if (t < s->next_output) {
    // Not reached output time yet => return
    PetscFunctionReturn(0);
  }

  PetscReal tfinal;
  static int i = 0;

#if PETSC_VERSION_GE(3, 8, 0)
  PetscCall(TSGetMaxTime(ts, &tfinal));
#else
  PetscCall(TSGetDuration(ts, nullptr, &tfinal));
#endif

  // Duplicate the solution vector X into a work vector
  Vec interpolatedX;
  PetscCall(VecDuplicate(X, &interpolatedX));

  // The internal timestepper may have stepped over multiple output times
  while (s->next_output <= t && s->next_output <= tfinal) {
    BoutReal output_time = t;
    if (s->interpolate) {
      int ierr = TSInterpolate(ts, s->next_output, interpolatedX);
      if (ierr != PETSC_SUCCESS) {
        throw BoutException("This PETSc TS does not support interpolation. Use a "
                            "different method or set solver:interpolate=false");
      }
      output_time = s->next_output;
    }

    // Place the interpolated values into the global variables
    const PetscScalar* x;
    PetscCall(VecGetArrayRead(interpolatedX, &x));
    s->load_vars(const_cast<BoutReal*>(x));
    PetscCall(VecRestoreArrayRead(interpolatedX, &x));

    if (s->call_monitors(output_time, i++, s->getNumberOutputSteps()) != 0) {
      PetscFunctionReturn(1);
    }

    s->next_output = output_time + s->getOutputTimestep();
  }

  // Done with vector, so destroy it
  PetscCall(VecDestroy(&interpolatedX));

  PetscFunctionReturn(0);
}

PetscSolver::PetscSolver(Options* opts)
    : Solver(opts), interpolate((*options)["interpolate"]
                                    .doc("Interpolate to regular output times?")
                                    .withDefault(true)),
      diagnose(
          (*options)["diagnose"].doc("Enable some diagnostic output").withDefault(false)),
      user_precon((*options)["user_precon"]
                      .doc("Use user-supplied preconditioning function?")
                      .withDefault(false)),
      atol((*options)["atol"].doc("Absolute tolerance").withDefault(1.0e-12)),
      rtol((*options)["rtol"].doc("Relative tolerance").withDefault(1.0e-5)),
      stol((*options)["stol"]
               .doc("Convergence tolerance in terms of the norm of the change in "
                    "the solution between steps")
               .withDefault(1e-8)),
      maxnl((*options)["max_nonlinear_iterations"]
                .doc("Maximum number of nonlinear iterations per SNES solve")
                .withDefault(50)),
      maxf((*options)["maxf"]
               .doc("Maximum number of function evaluations per SNES solve")
               .withDefault(10000)),
      maxl((*options)["maxl"].doc("Maximum number of linear iterations").withDefault(20)),
      ts_type((*options)["ts_type"].doc("PETSc time integrator type").withDefault("bdf")),
      adapt_type((*options)["adapt_type"]
                     .doc("PETSc TSAdaptType timestep adaptation method")
                     .withDefault("basic")),
      snes_type((*options)["snes_type"]
                    .doc("PETSc nonlinear solver method to use")
                    .withDefault("newtonls")),
      ksp_type((*options)["ksp_type"]
                   .doc("Linear solver type. By default let PETSc decide (gmres)")
                   .withDefault("default")),
      pc_type(
          (*options)["pc_type"]
              .doc("Preconditioner type. By default lets PETSc decide (ilu or bjacobi)")
              .withDefault("default")),
      pc_hypre_type((*options)["pc_hypre_type"]
                        .doc("hypre preconditioner type: euclid, pilut, parasails, "
                             "boomeramg, ams, ads")
                        .withDefault("pilut")),
      line_search_type((*options)["line_search_type"]
                           .doc("Line search type: basic, bt, l2, cp, nleqerr")
                           .withDefault("default")),
      matrix_free((*options)["matrix_free"]
                      .doc("Use matrix free Jacobian?")
                      .withDefault<bool>(false)),
      matrix_free_operator((*options)["matrix_free_operator"]
                               .doc("Use matrix free Jacobian-vector operator?")
                               .withDefault<bool>(false)),
      lag_jacobian((*options)["lag_jacobian"]
                       .doc("Re-use the Jacobian this number of SNES iterations")
                       .withDefault(50)),
      use_coloring((*options)["use_coloring"]
                       .doc("Use matrix coloring to calculate Jacobian?")
                       .withDefault<bool>(true)),
      kspsetinitialguessnonzero((*options)["kspsetinitialguessnonzero"]
                                    .doc("Set the initial guess to be non-zero")
                                    .withDefault<bool>(false)),
      start_timestep((*options)["start_timestep"]
                         .doc("Initial internal timestep (defaults to output timestep)")
                         .withDefault(getOutputTimestep())),
      mxstep(
          (*options)["mxstep"].doc("Number of steps between outputs").withDefault(500)) {}

PetscSolver::~PetscSolver() {
  VecDestroy(&u);
  MatDestroy(&Jfd);
  MatDestroy(&Jmf);
  MatFDColoringDestroy(&fdcoloring);
  TSDestroy(&ts);
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int PetscSolver::init() {

  TRACE("Initialising PETSc-dev solver");

  Solver::init();

  int nlocal = getLocalN(); // Number of evolving variables on this processor

  // Get total problem size
  int neq;
  if (bout::globals::mpi->MPI_Allreduce(&nlocal, &neq, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())
      != 0) {
    throw BoutException("MPI_Allreduce failed!");
  }

  output_info.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n",
                    n3Dvars(), n2Dvars(), neq, nlocal);

  // Create a vector to contain the state
  PetscCall(VecCreate(BoutComm::get(), &u));
  PetscCall(VecSetSizes(u, nlocal, PETSC_DECIDE));
  PetscCall(VecSetFromOptions(u));

  // Save initial state to PETSc Vec
  BoutReal* udata; // Pointer to data array in vector u.
  PetscCall(VecGetArray(u, &udata));
  save_vars(udata);
  PetscCall(VecRestoreArray(u, &udata));

  // Create timestepper
  PetscCall(TSCreate(BoutComm::get(), &ts));
  PetscCall(TSSetProblemType(ts, TS_NONLINEAR));
  PetscCall(TSSetType(ts, ts_type.c_str()));
  PetscCall(TSSetApplicationContext(ts, this));

  // Set user provided RHSFunction
  // Need to duplicate the solution vector for the residual
  Vec rhs_vec;
  PetscCall(VecDuplicate(u, &rhs_vec));
  PetscCall(TSSetRHSFunction(ts, rhs_vec, solver_rhs, this));

  // Set up adaptive time-stepping
  TSAdapt adapt;
  PetscCall(TSGetAdapt(ts, &adapt));
  PetscCall(TSAdaptSetType(adapt, adapt_type.c_str()));

  // Set default absolute/relative tolerances
  // Note: Vector atol and rtol not given
  PetscCall(TSSetTolerances(ts, atol, nullptr, rtol, nullptr));
  if (ts_type == TSSUNDIALS) {
#if PETSC_HAVE_SUNDIALS2
    // The PETSc interface to SUNDIALS' CVODE
    TSSundialsSetType(ts, SUNDIALS_BDF);
    TSSundialsSetTolerance(ts, atol, rtol);
#else
    throw BoutException("PETSc was not built with SUNDIALS. Reconfigure and build PETSc "
                        "with --download-sundials");
#endif
  }

  // Initial time of the simulation state
  PetscCall(TSSetTime(ts, simtime));
  next_output = simtime;

  PetscCall(TSSetTimeStep(ts, start_timestep));

  // Total number of steps
  PetscInt total_steps = mxstep * getNumberOutputSteps();
  // Final output time
  PetscReal tfinal = simtime + (getNumberOutputSteps() * getOutputTimestep());
  output.write("\tSet total_steps {:d}, tfinal {:g}, simtime {:g}\n", total_steps, tfinal,
               simtime);

#if PETSC_VERSION_GE(3, 8, 0)
  PetscCall(TSSetMaxSteps(ts, total_steps));
  PetscCall(TSSetMaxTime(ts, tfinal));
#else
  PetscCall(TSSetDuration(ts, total_steps, tfinal));
#endif

  // Allow TS to recover from SNES failures
  PetscCall(TSSetMaxSNESFailures(ts, PETSC_UNLIMITED));

  // Recover from step rejections
  PetscCall(TSSetMaxStepRejections(ts, PETSC_UNLIMITED));

  // Set the current solution
  PetscCall(TSSetSolution(ts, u));
  // Allow TS to step over the final time
  // Note: This does not affect intermediate outputs, that
  //       are always interpolated (in PetscMonitor)
  if (interpolate) {
    PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_INTERPOLATE));
  } else {
    PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
  }
  PetscCall(TSMonitorSet(ts, PetscMonitor, this, nullptr));

  if (ts_type != TSSUNDIALS) {
    // Get and configure the SNES nonlinear solver
    // Note: SUNDIALS does not use PETSc SNES or KSP
    // https://petsc.org/release/manual/ts/#using-sundials-from-petsc

    PetscCall(TSGetSNES(ts, &snes));
    PetscCall(SNESSetType(snes, snes_type.c_str()));

    // Line search
    if (line_search_type != "default") {
      SNESLineSearch linesearch;
      SNESGetLineSearch(snes, &linesearch);
      SNESLineSearchSetType(linesearch, line_search_type.c_str());
    }

    // Set tolerances
    // Note: TS should set SNES convergence tolerances
    SNESSetTolerances(snes,
                      PETSC_DETERMINE, // atol
                      PETSC_DETERMINE, // rtol
                      PETSC_DETERMINE, // stol
                      maxnl, maxf);

    // Force SNES to take at least one nonlinear iteration.
    // This may prevent the solver from getting stuck in false steady state conditions
#if PETSC_VERSION_GE(3, 8, 0)
    SNESSetForceIteration(snes, PETSC_TRUE);
#endif

    // Re-use Jacobian
    // Note: If the 'Amat' Jacobian is matrix free, SNESComputeJacobian
    //       always updates its reference 'u' vector every nonlinear iteration
    SNESSetLagJacobian(snes, lag_jacobian);
    SNESSetLagJacobianPersists(snes, PETSC_TRUE);

    SNESSetLagPreconditionerPersists(snes, PETSC_TRUE);
    SNESSetLagPreconditioner(snes, 1); // Rebuild when Jacobian is rebuilt

    // Get and configure the KSP linear solver
    SNESGetKSP(snes, &ksp);

    if (ksp_type != "default") {
      KSPSetType(ksp, ksp_type.c_str());
    }

    if (kspsetinitialguessnonzero) {
      // Set the initial guess to be non-zero
      KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    }

    KSPSetTolerances(ksp,
                     PETSC_DEFAULT, // rtol
                     PETSC_DEFAULT, // abstol
                     PETSC_DEFAULT, // dtol (divergence tolerance)
                     maxl);         // Maximum number of iterations

    // Set up the Jacobian
    // Note: We can have both matrix-free Jacobian-vector operator (Jmf)
    //       and a finite-difference Jacobian matrix (Jfd)
    if (matrix_free or matrix_free_operator) {
      /*
        PETSc SNES matrix free Jacobian, using a different
        operator for differencing.

        See PETSc examples
        http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tests/ex7.c.html
        and this thread:
        http://lists.mcs.anl.gov/pipermail/petsc-users/2014-January/020075.html

      */
      MatCreateSNESMF(snes, &Jmf);

      // Set a function to be called for differencing
      // This can be a linearised form of the SNES function
      // Note: Unsure if setting this will interfere with TS object
      //       because the SNES Jacobian depends on the timestep
      //MatMFFDSetFunction(Jmf, solver_rhs_differencing, this);
    }
  }

  // Configure preconditioner
  PC pc;
  if (ts_type == TSSUNDIALS) {
#if PETSC_HAVE_SUNDIALS2
    TSSundialsGetPC(ts, &pc);
#endif
  } else {
    // Get PC context from KSP
    KSPGetPC(ksp, &pc);
  }
  if (matrix_free) {
    // Matrix-free preconditioner

    if (user_precon) {
      output_info.write("\tUsing user-supplied preconditioner\n");
      if (!hasPreconditioner()) {
        throw BoutException("Model does not define a preconditioner");
      }
      // Set a Shell (matrix-free) preconditioner type
      PCSetType(pc, PCSHELL);

      // Specify the preconditioner function
      PCShellSetApply(pc, snesPCapply);

      // Context used to supply object pointer
      PCShellSetContext(pc, this);

      // Set name of preconditioner
      PetscCall(PCShellSetName(pc, "PhysicsPreconditioner"));
    } else {
      // Can't use preconditioner because no Jacobian matrix available
      PCSetType(pc, PCNONE);
    }

    // Callback to get the shift parameter
    // Use matrix free for both operator and preconditioner
    TSSetIJacobian(ts, Jmf, Jmf, solver_ijacobian, this);

  } else {
    // Set PC type from input
    if (pc_type != "default") {
      PCSetType(pc, pc_type.c_str());

      if (pc_type == "hypre") {
#if PETSC_HAVE_HYPRE
        // Set the type of hypre preconditioner
        PCHYPRESetType(pc, pc_hypre_type.c_str());
#else
        throw BoutException("PETSc was not configured with Hypre.");
#endif
      }
    }

    // Calculate a Jacobian matrix using finite differences. The finite
    // difference Jacobian (Jfd) may be used for both operator and
    // preconditioner or, if matrix_free_operator, in only the
    // preconditioner.

    if (use_coloring) {
      // Use matrix coloring.  This greatly reduces the number of
      // times the rhs() function needs to be evaluated when
      // calculating the Jacobian, by identifying which quantities may
      // be simultaneously perturbed.

      // Use global mesh for now
      Mesh* mesh = bout::globals::mesh;

      //////////////////////////////////////////////////
      // Get the local indices by starting at 0
      Field3D index = globalIndex(0);

      //////////////////////////////////////////////////
      // Pre-allocate PETSc storage

      output_progress.write("Setting Jacobian matrix sizes\n");

      const int n2d = f2d.size();
      const int n3d = f3d.size();

      // Set size of Matrix on each processor to nlocal x nlocal
      MatCreate(BoutComm::get(), &Jfd);
      MatSetOption(Jfd, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetSizes(Jfd, nlocal, nlocal, PETSC_DETERMINE, PETSC_DETERMINE);
      MatSetFromOptions(Jfd);
      // Determine which row/columns of the matrix are locally owned
      int Istart, Iend;
      MatGetOwnershipRange(Jfd, &Istart, &Iend);
      // Convert local into global indices
      // Note: Not in the boundary cells, to keep -1 values
      for (const auto& i : mesh->getRegion3D("RGN_NOBNDRY")) {
        index[i] += Istart;
      }
      // Now communicate to fill guard cells
      mesh->communicate(index);

      // Non-zero elements on this processor
      std::vector<PetscInt> d_nnz;
      std::vector<PetscInt> o_nnz;
      auto n_square = (*options)["stencil:square"]
                          .doc("Extent of stencil (square)")
                          .withDefault<int>(0);
      auto n_taxi = (*options)["stencil:taxi"]
                        .doc("Extent of stencil (taxi-cab norm)")
                        .withDefault<int>(0);
      auto n_cross = (*options)["stencil:cross"]
                         .doc("Extent of stencil (cross)")
                         .withDefault<int>(0);
      // Set n_taxi 2 if nothing else is set
      // Probably a better way to do this
      if (n_square == 0 && n_taxi == 0 && n_cross == 0) {
        output_info.write("Setting solver:stencil:taxi = 2\n");
        n_taxi = 2;
      }

      auto const xy_offsets = ColoringStencil::getOffsets(n_square, n_taxi, n_cross);
      {
        // This is ugly but can't think of a better and robust way to
        // count the non-zeros for some arbitrary stencil
        // effectively the same loop as the one that sets the non-zeros below
        std::vector<std::set<int>> d_nnz_map2d(nlocal);
        std::vector<std::set<int>> o_nnz_map2d(nlocal);
        std::vector<std::set<int>> d_nnz_map3d(nlocal);
        std::vector<std::set<int>> o_nnz_map3d(nlocal);
        // Loop over every element in 2D to count the *unique* non-zeros
        for (int x = mesh->xstart; x <= mesh->xend; x++) {
          for (int y = mesh->ystart; y <= mesh->yend; y++) {

            const int ind0 = ROUND(index(x, y, 0)) - Istart;

            // 2D fields
            for (int i = 0; i < n2d; i++) {
              const PetscInt row = ind0 + i;
              // Loop through each point in the stencil
              for (const auto& [x_off, y_off] : xy_offsets) {
                const int xi = x + x_off;
                const int yi = y + y_off;
                if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx)
                    || (yi >= mesh->LocalNy)) {
                  continue;
                }

                const int ind2 = ROUND(index(xi, yi, 0));
                if (ind2 < 0) {
                  continue; // A boundary point
                }

                // Depends on all variables on this cell
                for (int j = 0; j < n2d; j++) {
                  const PetscInt col = ind2 + j;
                  if (col >= Istart && col < Iend) {
                    d_nnz_map2d[row].insert(col);
                  } else {
                    o_nnz_map2d[row].insert(col);
                  }
                }
              }
            }
            // 3D fields
            for (int z = 0; z < mesh->LocalNz; z++) {
              const int ind = ROUND(index(x, y, z)) - Istart;

              for (int i = 0; i < n3d; i++) {
                PetscInt row = ind + i;
                if (z == 0) {
                  row += n2d;
                }

                // Depends on 2D fields
                for (int j = 0; j < n2d; j++) {
                  const PetscInt col = ind0 + j;
                  if (col >= Istart && col < Iend) {
                    d_nnz_map2d[row].insert(col);
                  } else {
                    o_nnz_map2d[row].insert(col);
                  }
                }

                // Star pattern
                for (const auto& [x_off, y_off] : xy_offsets) {
                  const int xi = x + x_off;
                  const int yi = y + y_off;

                  if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx)
                      || (yi >= mesh->LocalNy)) {
                    continue;
                  }

                  int ind2 = ROUND(index(xi, yi, 0));
                  if (ind2 < 0) {
                    continue; // Boundary point
                  }

                  if (z == 0) {
                    ind2 += n2d;
                  }

                  // 3D fields on this cell
                  for (int j = 0; j < n3d; j++) {
                    const PetscInt col = ind2 + j;
                    if (col >= Istart && col < Iend) {
                      d_nnz_map3d[row].insert(col);
                    } else {
                      o_nnz_map3d[row].insert(col);
                    }
                  }
                }
              }
            }
          }
        }

        d_nnz.reserve(nlocal);
        d_nnz.reserve(nlocal);

        for (int i = 0; i < nlocal; ++i) {
          // Assume all elements in the z direction are potentially coupled
          d_nnz.emplace_back((d_nnz_map3d[i].size() * mesh->LocalNz)
                             + d_nnz_map2d[i].size());
          o_nnz.emplace_back((o_nnz_map3d[i].size() * mesh->LocalNz)
                             + o_nnz_map2d[i].size());
        }
      }

      output_progress.write("Pre-allocating Jacobian\n");
      // Pre-allocate
      MatMPIAIJSetPreallocation(Jfd, 0, d_nnz.data(), 0, o_nnz.data());
      MatSeqAIJSetPreallocation(Jfd, 0, d_nnz.data());
      MatSetUp(Jfd);
      MatSetOption(Jfd, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

      //////////////////////////////////////////////////
      // Mark non-zero entries

      output_progress.write("Marking non-zero Jacobian entries\n");
      PetscScalar const val = 1.0;
      for (int x = mesh->xstart; x <= mesh->xend; x++) {
        for (int y = mesh->ystart; y <= mesh->yend; y++) {

          const int ind0 = ROUND(index(x, y, 0));

          // 2D fields
          for (int i = 0; i < n2d; i++) {
            const PetscInt row = ind0 + i;

            // Loop through each point in the stencil
            for (const auto& [x_off, y_off] : xy_offsets) {
              const int xi = x + x_off;
              const int yi = y + y_off;
              if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx)
                  || (yi >= mesh->LocalNy)) {
                continue;
              }

              int const ind2 = ROUND(index(xi, yi, 0));
              if (ind2 < 0) {
                continue; // A boundary point
              }

              // Depends on all variables on this cell
              for (int j = 0; j < n2d; j++) {
                PetscInt const col = ind2 + j;
                PetscCall(MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES));
              }
            }
          }
          // 3D fields
          for (int z = 0; z < mesh->LocalNz; z++) {
            int const ind = ROUND(index(x, y, z));

            for (int i = 0; i < n3d; i++) {
              PetscInt row = ind + i;
              if (z == 0) {
                row += n2d;
              }

              // Depends on 2D fields
              for (int j = 0; j < n2d; j++) {
                PetscInt const col = ind0 + j;
                PetscCall(MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES));
              }

              // Star pattern
              for (const auto& [x_off, y_off] : xy_offsets) {
                int xi = x + x_off;
                int yi = y + y_off;

                if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx)
                    || (yi >= mesh->LocalNy)) {
                  continue;
                }
                for (int zi = 0; zi < mesh->LocalNz; ++zi) {
                  int ind2 = ROUND(index(xi, yi, zi));
                  if (ind2 < 0) {
                    continue; // Boundary point
                  }

                  if (z == 0) {
                    ind2 += n2d;
                  }

                  // 3D fields on this cell
                  for (int j = 0; j < n3d; j++) {
                    PetscInt const col = ind2 + j;
                    int ierr = MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES);

                    if (ierr != 0) {
                      output.write("ERROR: {} {} : ({}, {}) -> ({}, {}) : {} -> {}\n",
                                   row, x, y, xi, yi, ind2, ind2 + n3d - 1);
                    }
                    CHKERRQ(ierr);
                  }
                }
              }
            }
          }
        }
      }

      // Finished marking non-zero entries

      output_progress.write("Assembling Jacobian matrix\n");

      // Assemble Matrix
      MatAssemblyBegin(Jfd, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(Jfd, MAT_FINAL_ASSEMBLY);

      {
        // Test if the matrix is symmetric
        // Values are 0 or 1 so tolerance (1e-5) shouldn't matter
        PetscBool symmetric;
        PetscCall(MatIsSymmetric(Jfd, 1e-5, &symmetric));
        if (!symmetric) {
          output_warn.write("Jacobian pattern is not symmetric\n");
        }
      }

      // The above can miss entries around the X-point branch cut:
      // The diagonal terms are complicated because moving in X then Y
      // is different from moving in Y then X at the X-point.
      // Making sure the colouring matrix is symmetric does not
      // necessarily give the correct stencil but may help.
      if ((*options)["force_symmetric_coloring"]
              .doc("Modifies coloring matrix to force it to be symmetric")
              .withDefault<bool>(false)) {
        Mat Jfd_T;
        MatCreateTranspose(Jfd, &Jfd_T);
        MatAXPY(Jfd, 1, Jfd_T, DIFFERENT_NONZERO_PATTERN);
      }

      output_progress.write("Creating Jacobian coloring\n");
      updateColoring();

    } else {
      // Brute force calculation
      // There is usually no reason to use this, except as a check of
      // the coloring calculation.

      MatCreateAIJ(
          BoutComm::get(), nlocal, nlocal,  // Local sizes
          PETSC_DETERMINE, PETSC_DETERMINE, // Global sizes
          3, // Number of nonzero entries in diagonal portion of local submatrix
          nullptr,
          0, // Number of nonzeros per row in off-diagonal portion of local submatrix
          nullptr, &Jfd);

      if (matrix_free_operator) {
        SNESSetJacobian(snes, Jmf, Jfd, SNESComputeJacobianDefault, this);
      } else {
        SNESSetJacobian(snes, Jfd, Jfd, SNESComputeJacobianDefault, this);
      }

      MatSetOption(Jfd, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
  }

  // Get runtime options
  lib.setOptionsFromInputFile(ts);

  {
    // Some reporting
    TSType tstype;
    TSGetType(ts, &tstype);
    output_info.write("TS Type : {}\n", tstype);
    TSAdaptType adapttype;
    TSAdaptGetType(adapt, &adapttype);
    output_info.write("TS Adapt Type : {}\n", adapttype);
    if (ts_type != TSSUNDIALS) {
      SNESType snestype;
      SNESGetType(snes, &snestype);
      output_info.write("SNES Type : {}\n", snestype);
      KSPType ksptype;
      KSPGetType(ksp, &ksptype);
      if (ksptype) {
        output_info.write("KSP Type  : {}\n", ksptype);
      }
    }
    PCType pctype;
    PCGetType(pc, &pctype);
    if (pctype) {
      output_info.write("PC Type   : {}\n", pctype);
    }
  }
  return 0;
}

// Starts the time integrator
PetscErrorCode PetscSolver::run() {
  // Set when the next call to monitor is desired
  next_output = simtime + getOutputTimestep();

  PetscFunctionBegin;

  // Run the PETSc time integrator
  // The PetscMonitor function is responsible for regular outputs
  CHKERRQ(TSSolve(ts, u));

  PetscFunctionReturn(0);
}

// Evaluate the user-supplied RHS function
// This method is called from the static PETSc callback functions
PetscErrorCode PetscSolver::rhs(BoutReal t, Vec udata, Vec dudata, bool linear) {
  TRACE("Running RHS: PetscSolver::rhs({:e})", t);

  simtime = t;

  PetscFunctionBegin;

  // Load state from PETSc
  const BoutReal* udata_array;
  VecGetArrayRead(udata, &udata_array);
  load_vars(const_cast<BoutReal*>(udata_array));
  VecRestoreArrayRead(udata, &udata_array);

  // Call RHS function
  run_rhs(t, linear);

  // Save derivatives to PETSc
  BoutReal* dudata_array;
  VecGetArray(dudata, &dudata_array);
  save_derivs(dudata_array);
  VecRestoreArray(dudata, &dudata_array);

  PetscFunctionReturn(0);
}

PetscErrorCode PetscSolver::formFunction(Vec U, Vec F) {
  PetscFunctionBegin;
  PetscCall(rhs(simtime, U, F, true)); // Can be linearised for coloring

  // shift * U - dU/dt
  PetscCall(VecAXPBY(F, shift, -1.0, U));
  PetscFunctionReturn(0);
}

// Matrix-free preconditioner function
PetscErrorCode PetscSolver::pre(Vec x, Vec y) {
  TRACE("PetscSolver::pre()");

  BoutReal* data;

  // Load state
  VecGetArray(state, &data);
  load_vars(data);
  VecRestoreArray(state, &data);

  // Load vector to be inverted into F_vars
  VecGetArray(x, &data);
  load_derivs(data);
  VecRestoreArray(x, &data);

  // Call the preconditioner
  runPreconditioner(ts_time, 1. / shift, 0.0);

  // Save the solution from time derivatives
  VecGetArray(y, &data);
  save_derivs(data);
  VecRestoreArray(y, &data);

  // Petsc's definition of Jacobian differs by a factor from SUNDIALS'
  // PETSc solves (scale + J)^-1
  // SUNDIALS solves (I + gamma J)^-1
  PetscCall(VecScale(y, shift));

  return 0;
}

void PetscSolver::updateColoring() {
  // Re-calculate the coloring
  MatColoring coloring{nullptr};
  MatColoringCreate(Jfd, &coloring);
  MatColoringSetType(coloring, MATCOLORINGGREEDY);
  MatColoringSetFromOptions(coloring);

  // Calculate new index sets
  ISColoring iscoloring{nullptr};
  MatColoringApply(coloring, &iscoloring);
  MatColoringDestroy(&coloring);

  // Replace the old coloring with the new one
  MatFDColoringDestroy(&fdcoloring);
  MatFDColoringCreate(Jfd, iscoloring, &fdcoloring);
  if (ts_type != TSSUNDIALS) {
    // Use the SNES function that is defined by the TS method
    // SNESTSFormFunction is defined in PETSc ts.c
    // The ctx pointer should be the TS object
    //
    // Note: The cast is horrible but the function signature
    //       varies between PETSc versions in ways that break
    //       reinterpret_cast and a C cast to MatFDColoringFn.
    //       This is the cast that PETSc examples use.
    MatFDColoringSetFunction(fdcoloring, (PetscErrorCode(*)(void))SNESTSFormFunction, ts);
  } else {
    // SNESTSFormFunction is not available for SUNDIALS.
    // This solver_form_function needs to know the shift
    // (SUNDIALS' gamma) that we capture in solver_ijacobian_color.
    MatFDColoringSetFunction(fdcoloring, (PetscErrorCode(*)(void))solver_form_function,
                             this);
  }
  MatFDColoringSetFromOptions(fdcoloring);
  MatFDColoringSetUp(Jfd, iscoloring, fdcoloring);
  ISColoringDestroy(&iscoloring);

  // Replace the CTX pointer in SNES Jacobian
  if (ts_type == TSSUNDIALS) {
#if PETSC_HAVE_SUNDIALS2
    // The SUNDIALS interface calls TSGetIJacobian
    // https://www.mcs.anl.gov/petsc/petsc-3.14/src/ts/impls/implicit/sundials/sundials.c
    // This sets the "TSMatFDColoring" property on the Jacobian, that is used in TSComputeIJacobianDefaultColor
    PetscObjectCompose((PetscObject)Jfd, "TSMatFDColoring", (PetscObject)fdcoloring);
    // Call a wrapper function that stores the shift in the PetscSolver
    TSSetIJacobian(ts, Jfd, Jfd, solver_ijacobian_color, this);
#endif
  } else {
    if (matrix_free_operator) {
      // Use matrix-free calculation for operator, finite difference for preconditioner
      SNESSetJacobian(snes, Jmf, Jfd, SNESComputeJacobianDefaultColor, fdcoloring);
    } else {
      SNESSetJacobian(snes, Jfd, Jfd, SNESComputeJacobianDefaultColor, fdcoloring);
    }
  }
}

#endif
