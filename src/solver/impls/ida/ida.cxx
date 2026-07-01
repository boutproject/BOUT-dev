/**************************************************************************
 * Interface to SUNDIALS IDA
 *
 * IdaSolver for DAE systems (so can handle constraints)
 *
 * NOTE: Only one solver can currently be compiled in
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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

#include "bout/build_defines.hxx"

#include "ida.hxx"

#if BOUT_HAS_IDA

#include "../../sundials_nvector_interface.hxx"

#include "bout/bout_types.hxx"
#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/mpi_wrapper.hxx"
#include "bout/msg_stack.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"
#include "bout/solver.hxx"
#include "bout/sundials_backports.hxx"
#include "bout/unused.hxx"

#include <ida/ida.h>
#include <ida/ida_bbdpre.h>
#include <ida/ida_ls.h>

#include <iterator>
#include <numeric>
#include <vector>

// NOLINTBEGIN(readability-identifier-length)
namespace {
int idares(BoutReal t, N_Vector u, N_Vector du, N_Vector rr, void* user_data);
int ida_bbd_res(sunindextype Nlocal, BoutReal t, N_Vector u, N_Vector du, N_Vector rr,
                void* user_data);

int ida_pre(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector rvec,
            N_Vector zvec, BoutReal cj, BoutReal delta, void* user_data);
} // namespace
// NOLINTEND(readability-identifier-length)

IdaSolver::IdaSolver(Options* opts)
    : Solver(opts),
      abstol((*options)["atol"].doc("Absolute tolerance").withDefault(1.0e-12)),
      reltol((*options)["rtol"].doc("Relative tolerance").withDefault(1.0e-5)),
      mxsteps((*options)["mxstep"]
                  .doc("Maximum number of steps to take between outputs")
                  .withDefault(500)),
      use_precon((*options)["use_precon"]
                     .doc("Use user-supplied preconditioner")
                     .withDefault(false)),
      correct_start((*options)["correct_start"]
                        .doc("Correct the initial values")
                        .withDefault(true)),
      nvector_type((*options)["nvector"]
                       .doc("N_Vector backend to use: sundials or manyvector")
                       .withDefault(NVectorType::Sundials)),
      suncontext(createSUNContext(BoutComm::get())) {
  has_constraints = true; // This solver has constraints
}

IdaSolver::~IdaSolver() {
  if (initialised) {
    N_VDestroy(uvec);
    N_VDestroy(duvec);
    N_VDestroy(id);
    IDAFree(&idamem);
    SUNLinSolFree(sun_solver);
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int IdaSolver::init() {
  const auto backend = nvector_backend();

  Solver::init();
  backend.ensure_manyvector_available();

  output.write("Initialising IDA solver\n");

  // Calculate number of variables
  const int n2d = f2d.size();
  const int n3d = f3d.size();
  const int local_N = getLocalN();

  // Get total problem size
  int neq;
  if (bout::globals::mpi->MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    output_error.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }

  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n", n3d, n2d,
               neq, local_N);
  output.write("\tUsing {} N_Vector backend\n", backend.backend_name());

  // Allocate memory
  uvec = backend.create_state_vector(local_N, neq);
  duvec = backend.clone_vector_like(uvec, "SUNDIALS memory allocation failed");
  id = backend.clone_vector_like(uvec, "SUNDIALS memory allocation failed");

  // Get the starting time derivative
  run_rhs(simtime);

  // Put the time-derivatives into duvec
  backend.copy_deriv_to_vector(duvec);

  // Set the equation type in id(Differential or Algebraic. This is optional)
  backend.fill_id_vector(id);

  // Call IDACreate to initialise
  idamem = callWithSUNContext(IDACreate, suncontext);
  if (idamem == nullptr) {
    throw BoutException("IDACreate failed\n");
  }

  // For callbacks, need pointer to solver object
  if (IDASetUserData(idamem, this) != IDA_SUCCESS) {
    throw BoutException("IDASetUserData failed\n");
  }

  if (IDASetId(idamem, id) != IDA_SUCCESS) {
    throw BoutException("IDASetID failed\n");
  }

  if (IDAInit(idamem, idares, simtime, uvec, duvec) != IDA_SUCCESS) {
    throw BoutException("IDAInit failed\n");
  }

  if (IDASStolerances(idamem, reltol, abstol) != IDA_SUCCESS) {
    throw BoutException("IDASStolerances failed\n");
  }

  if (IDASetMaxNumSteps(idamem, mxsteps) != IDA_SUCCESS) {
    throw BoutException("IDASetMaxNumSteps failed\n");
  }

  // Call IDASpgmr to specify the IDA linear solver IDASPGMR
  const auto maxl = (*options)["maxl"].withDefault(6 * n3d);
  sun_solver = callWithSUNContext(SUNLinSol_SPGMR, suncontext, uvec, SUN_PREC_NONE, maxl);
  if (sun_solver == nullptr) {
    throw BoutException("Creating SUNDIALS linear solver failed\n");
  }
  if (IDASetLinearSolver(idamem, sun_solver, nullptr) != IDALS_SUCCESS) {
    throw BoutException("IDASetLinearSolver failed\n");
  }

  if (use_precon) {
    if (hasPreconditioner()) {
      output.write("\tUsing user-supplied preconditioner\n");
      if (IDASetPreconditioner(idamem, nullptr, ida_pre) != IDALS_SUCCESS) {
        throw BoutException("IDASetPreconditioner failed\n");
      }
    } else {
      output.write("\tUsing BBD preconditioner\n");
      /// Get options
      // Compute band_width_default from actually added fields, to allow for multiple Mesh
      // objects
      //
      // Previous implementation was equivalent to:
      //   int MXSUB = mesh->xend - mesh->xstart + 1;
      //   int band_width_default = n3Dvars()*(MXSUB+2);
      const int band_width_default = std::accumulate(
          begin(f3d), end(f3d), 0, [](int a, const VarStr<Field3D>& fvar) {
            Mesh* localmesh = fvar.var->getMesh();
            return a + localmesh->xend - localmesh->xstart + 3;
          });

      const auto mudq = (*options)["mudq"].withDefault(band_width_default);
      const auto mldq = (*options)["mldq"].withDefault(band_width_default);
      const auto mukeep = (*options)["mukeep"].withDefault(n3d);
      const auto mlkeep = (*options)["mlkeep"].withDefault(n3d);
      if (IDABBDPrecInit(idamem, local_N, mudq, mldq, mukeep, mlkeep, 0.0, ida_bbd_res,
                         nullptr)
          != IDALS_SUCCESS) {
        throw BoutException("IDABBDPrecInit failed\n");
      }
    }
  }

  // Call IDACalcIC (with default options) to correct the initial values
  if (correct_start) {
    if (IDACalcIC(idamem, IDA_YA_YDP_INIT, 1e-6) != IDA_SUCCESS) {
      throw BoutException("IDACalcIC failed\n");
    }
  }

  return 0;
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int IdaSolver::run() {

  if (!initialised) {
    throw BoutException("IdaSolver not initialised\n");
  }

  for (int i = 1; i <= getNumberOutputSteps(); i++) {

    /// Run the solver for one output timestep
    simtime = run(simtime + getOutputTimestep());

    /// Check if the run succeeded
    if (simtime < 0.0) {
      // Step failed
      throw BoutException("SUNDIALS IDA timestep failed\n");
    }

    if (call_monitors(simtime, i, getNumberOutputSteps())) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

BoutReal IdaSolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run({:e})", tout);
  const auto backend = nvector_backend();

  if (!initialised) {
    throw BoutException("Running IDA solver without initialisation\n");
  }

  pre_Wtime = 0.0;
  pre_ncalls = 0;

  const int flag = IDASolve(idamem, tout, &simtime, uvec, duvec, IDA_NORMAL);

  // Copy variables
  backend.copy_state_from_vector(uvec);

  // Call rhs function to get extra variables at this time
  run_rhs(simtime);

  if (flag < 0) {
    output_error.write("ERROR IDA solve failed at t = {:e}, flag = {:d}\n", simtime,
                       flag);
    return -1.0;
  }

  return simtime;
}

/**************************************************************************
 * Residual function F(t, u, du)
 **************************************************************************/

void IdaSolver::res(BoutReal t, N_Vector u, N_Vector du, N_Vector rr) {
  TRACE("Running RHS: IdaSolver::res({:e})", t);
  const auto backend = nvector_backend();

  // Load state from udata
  backend.copy_state_from_vector(u);

  // Call RHS function
  run_rhs(t);

  // Save derivatives to rdata (residual)
  backend.copy_state_to_vector(u);
  backend.copy_deriv_to_vector(rr);
  backend.subtract_differential_term(rr, du, id);
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void IdaSolver::pre(BoutReal t, BoutReal cj, BoutReal delta, N_Vector u, N_Vector rvec,
                    N_Vector zvec) {
  TRACE("Running preconditioner: IdaSolver::pre({:e})", t);
  const auto backend = nvector_backend();

  const BoutReal tstart = bout::globals::mpi->MPI_Wtime();

  if (!hasPreconditioner()) {
    // Identity (but should never happen)
    N_VScale(1.0, rvec, zvec);
    return;
  }

  // Load state from udata (as with res function)
  backend.copy_state_from_vector(u);

  // Load vector to be inverted into F_vars
  backend.copy_deriv_from_vector(rvec);

  runPreconditioner(t, cj, delta);

  // Save the solution from F_vars
  backend.copy_state_to_vector(u);
  backend.copy_deriv_to_vector(zvec);

  pre_Wtime += bout::globals::mpi->MPI_Wtime() - tstart;
  pre_ncalls++;
}

/**************************************************************************
 * IDA res function
 **************************************************************************/

// NOLINTBEGIN(readability-identifier-length)
namespace {
int idares(BoutReal t, N_Vector u, N_Vector du, N_Vector rr, void* user_data) {
  auto* s = static_cast<IdaSolver*>(user_data);

  // Calculate residuals
  s->res(t, u, du, rr);

  return 0;
}

/// Residual function for BBD preconditioner
int ida_bbd_res(sunindextype UNUSED(Nlocal), BoutReal t, N_Vector u, N_Vector du,
                N_Vector rr, void* user_data) {
  return idares(t, u, du, rr, user_data);
}

// Preconditioner function
int ida_pre(BoutReal t, N_Vector yy, N_Vector UNUSED(yp), N_Vector UNUSED(rr),
            N_Vector rvec, N_Vector zvec, BoutReal cj, BoutReal delta, void* user_data) {
  auto* s = static_cast<IdaSolver*>(user_data);

  // Calculate residuals
  s->pre(t, cj, delta, yy, rvec, zvec);

  return 0;
}
} // namespace
// NOLINTEND(readability-identifier-length)

#endif
