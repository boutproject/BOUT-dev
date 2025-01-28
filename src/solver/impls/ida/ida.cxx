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

#include <algorithm>
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

  TRACE("Initialising IDA solver");

  Solver::init();

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

  // Allocate memory
  uvec = callWithSUNContext(N_VNew_Parallel, suncontext, BoutComm::get(), local_N, neq);
  if (uvec == nullptr) {
    throw BoutException("SUNDIALS memory allocation failed\n");
  }
  duvec = N_VClone(uvec);
  if (duvec == nullptr) {
    throw BoutException("SUNDIALS memory allocation failed\n");
  }
  id = N_VClone(uvec);
  if (id == nullptr) {
    throw BoutException("SUNDIALS memory allocation failed\n");
  }

  // Put the variables into uvec
  save_vars(N_VGetArrayPointer(uvec));

  // Get the starting time derivative
  run_rhs(simtime);

  // Put the time-derivatives into duvec
  save_derivs(N_VGetArrayPointer(duvec));

  // Set the equation type in id(Differential or Algebraic. This is optional)
  set_id(N_VGetArrayPointer(id));

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
  TRACE("IDA IdaSolver::run()");

  if (!initialised) {
    throw BoutException("IdaSolver not initialised\n");
  }

  for (int i = 0; i < getNumberOutputSteps(); i++) {

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

  if (!initialised) {
    throw BoutException("Running IDA solver without initialisation\n");
  }

  pre_Wtime = 0.0;
  pre_ncalls = 0;

  const int flag = IDASolve(idamem, tout, &simtime, uvec, duvec, IDA_NORMAL);

  // Copy variables
  load_vars(N_VGetArrayPointer(uvec));

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

void IdaSolver::res(BoutReal t, BoutReal* udata, BoutReal* dudata, BoutReal* rdata) {
  TRACE("Running RHS: IdaSolver::res({:e})", t);

  // Load state from udata
  load_vars(udata);

  // Call RHS function
  run_rhs(t);

  // Save derivatives to rdata (residual)
  save_derivs(rdata);

  // If a differential equation, subtract dudata
  const auto length = N_VGetLocalLength_Parallel(id);
  const BoutReal* idd = N_VGetArrayPointer(id);
  for (int i = 0; i < length; i++) {
    if (idd[i] > 0.5) { // 1 -> differential, 0 -> algebraic
      rdata[i] -= dudata[i];
    }
  }
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void IdaSolver::pre(BoutReal t, BoutReal cj, BoutReal delta, BoutReal* udata,
                    BoutReal* rvec, BoutReal* zvec) {
  TRACE("Running preconditioner: IdaSolver::pre({:e})", t);

  const BoutReal tstart = bout::globals::mpi->MPI_Wtime();

  if (!hasPreconditioner()) {
    // Identity (but should never happen)
    const auto length = N_VGetLocalLength_Parallel(id);
    std::copy(rvec, rvec + length, zvec);
    return;
  }

  // Load state from udata (as with res function)
  load_vars(udata);

  // Load vector to be inverted into F_vars
  load_derivs(rvec);

  runPreconditioner(t, cj, delta);

  // Save the solution from F_vars
  save_derivs(zvec);

  pre_Wtime += bout::globals::mpi->MPI_Wtime() - tstart;
  pre_ncalls++;
}

/**************************************************************************
 * IDA res function
 **************************************************************************/

// NOLINTBEGIN(readability-identifier-length)
namespace {
int idares(BoutReal t, N_Vector u, N_Vector du, N_Vector rr, void* user_data) {
  BoutReal* udata = N_VGetArrayPointer(u);
  BoutReal* dudata = N_VGetArrayPointer(du);
  BoutReal* rdata = N_VGetArrayPointer(rr);

  auto* s = static_cast<IdaSolver*>(user_data);

  // Calculate residuals
  s->res(t, udata, dudata, rdata);

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
  BoutReal* udata = N_VGetArrayPointer(yy);
  BoutReal* rdata = N_VGetArrayPointer(rvec);
  BoutReal* zdata = N_VGetArrayPointer(zvec);

  auto* s = static_cast<IdaSolver*>(user_data);

  // Calculate residuals
  s->pre(t, cj, delta, udata, rdata, zdata);

  return 0;
}
} // namespace
// NOLINTEND(readability-identifier-length)

#endif
