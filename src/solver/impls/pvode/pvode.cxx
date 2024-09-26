/**************************************************************************
 * Interface to PVODE solver
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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

#include "bout/build_config.hxx"

#include "pvode.hxx"

#if BOUT_HAS_PVODE

#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/mesh.hxx>
#include <bout/msg_stack.hxx>
#include <bout/output.hxx>
#include <bout/sys/timer.hxx>

#include "bout/unused.hxx"

#include <pvode/cvspgmr.h>  // use CVSPGMR linear solver each internal step
#include <pvode/iterativ.h> // contains the enum for types of preconditioning
#include <pvode/pvbbdpre.h> // band preconditioner function prototypes

#include <string>

using namespace pvode;

void solver_f(integer N, BoutReal t, N_Vector u, N_Vector udot, void* f_data);
void solver_gloc(integer N, BoutReal t, BoutReal* u, BoutReal* udot, void* f_data);
void solver_cfn(integer N, BoutReal t, N_Vector u, void* f_data);

namespace {
// local only
void pvode_load_data_f3d(const std::vector<bool>& evolve_bndrys,
                         std::vector<Field3D>& ffs, const BoutReal* udata) {
  int p = 0;
  Mesh* mesh = ffs[0].getMesh();
  const int nz = mesh->LocalNz;
  for (const auto& bndry : {true, false}) {
    for (const auto& i2d : mesh->getRegion2D(bndry ? "RGN_BNDRY" : "RGN_NOBNDRY")) {
      for (int jz = 0; jz < nz; jz++) {
        // Loop over 3D variables
        std::vector<bool>::const_iterator evolve_bndry = evolve_bndrys.begin();
        for (std::vector<Field3D>::iterator ff = ffs.begin(); ff != ffs.end(); ++ff) {
          if (bndry && !*evolve_bndry) {
            continue;
          }
          (*ff)[mesh->ind2Dto3D(i2d, jz)] = udata[p];
          p++;
        }
        ++evolve_bndry;
      }
    }
  }
}
} // namespace

const BoutReal ZERO = 0.0;

long int iopt[OPT_SIZE];
BoutReal ropt[OPT_SIZE];

PvodeSolver::PvodeSolver(Options* opts)
    : Solver(opts), use_precon((*options)["use_precon"]
                                   .doc("Use user-supplied preconditioner")
                                   .withDefault(false)),
      precon_dimens(
          (*options)["precon_dimens"].doc("Maximum Krylov dimension").withDefault(50)),
      precon_tol((*options)["precon_tol"]
                     .doc("Tolerance for preconditioner")
                     .withDefault(1.0e-4)),
      pvode_mxstep(
          (*options)["mxstep"].doc("Maximum number of internal steps").withDefault(500)),
      abstol((*options)["atol"].doc("Absolute tolerance").withDefault(1.0e-12)),
      reltol((*options)["rtol"].doc("Relative tolerance").withDefault(1.0e-5)) {
  has_constraints = false; ///< This solver doesn't have constraints
}

PvodeSolver::~PvodeSolver() {
  if (pvode_initialised) {
    // Free CVODE memory

    N_VFree(u);
    PVBBDFree(pdata);
    CVodeFree(cvode_mem);
    PVecFreeMPI(machEnv);
  }
}

#if SUNDIALS_VERSION_MAJOR >= 6
#else
#define SUN_MODIFIED_GS MODIFIED_GS
#endif

/**************************************************************************
 * Initialise
 **************************************************************************/

int PvodeSolver::init() {
  TRACE("Initialising PVODE solver");

  int mudq, mldq, mukeep, mlkeep;
  boole optIn;
  int i;

  int n2d = n2Dvars(); // Number of 2D variables
  int n3d = n3Dvars(); // Number of 3D variables

  Solver::init();
  output.write("Initialising PVODE solver\n");

  int local_N = getLocalN();

  if (local_N == 0) {
    throw BoutException("No local evolving variables");
  }

  // Get total problem size
  int neq;
  if (bout::globals::mpi->MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    throw BoutException("\tERROR: MPI_Allreduce failed!\n");
  }

  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n", n3d, n2d,
               neq, local_N);

  // Set machEnv block
  machEnv =
      static_cast<machEnvType>(PVecInitMPI(BoutComm::get(), local_N, neq, pargc, pargv));

  if (machEnv == nullptr) {
    throw BoutException("\tError: PVecInitMPI failed\n");
  }

  // Allocate memory, and set problem data, initial values, tolerances

  u = N_VNew(neq, machEnv);

  ///////////// GET OPTIONS /////////////

  // Compute band_width_default from actually added fields, to allow for multiple Mesh objects
  //
  // Previous implementation was equivalent to:
  //   int MXSUB = mesh->xend - mesh->xstart + 1;
  //   int band_width_default = n3Dvars()*(MXSUB+2);
  int band_width_default = 0;
  for (const auto& fvar : f3d) {
    Mesh* localmesh = fvar.var->getMesh();
    band_width_default += localmesh->xend - localmesh->xstart + 3;
  }

  options->get("mudq", mudq, band_width_default);
  options->get("mldq", mldq, band_width_default);
  options->get("mukeep", mukeep, 0);
  options->get("mlkeep", mlkeep, 0);

  pdata = PVBBDAlloc(local_N, mudq, mldq, mukeep, mlkeep, ZERO, solver_gloc, solver_cfn,
                     static_cast<void*>(this));

  if (pdata == nullptr) {
    throw BoutException("\tError: PVBBDAlloc failed.\n");
  }

  ////////// SAVE DATA TO CVODE ///////////

  // Set pointer to data array in vector u.
  BoutReal* udata = N_VDATA(u);
  save_vars(udata);

  /* Call CVodeMalloc to initialize CVODE: 
     
     neq     is the problem size = number of equations
     f       is the user's right hand side function in y'=f(t,y)
     T0      is the initial time
     u       is the initial dependent variable vector
     BDF     specifies the Backward Differentiation Formula
     NEWTON  specifies a Newton iteration
     SS      specifies scalar relative and absolute tolerances
     &reltol and &abstol are pointers to the scalar tolerances
     data    is the pointer to the user-defined data block
     NULL    is the pointer to the error message file
     FALSE   indicates there are no optional inputs in iopt and ropt
     iopt, ropt  communicate optional integer and BoutReal input/output

     A pointer to CVODE problem memory is returned and stored in cvode_mem.  */

  optIn = TRUE;
  for (i = 0; i < OPT_SIZE; i++) {
    iopt[i] = 0;
  }
  for (i = 0; i < OPT_SIZE; i++) {
    ropt[i] = ZERO;
  }
  /* iopt[MXSTEP] : maximum number of internal steps to be taken by *
   *                the solver in its attempt to reach tout.        *
   *                Optional input. (Default = 500).                */
  iopt[MXSTEP] = pvode_mxstep;

  {
    /* ropt[H0]      : initial step size. Optional input.             */

    /* ropt[HMAX]    : maximum absolute value of step size allowed.   *
     *                 Optional input. (Default is infinity).         */
    const BoutReal hmax(
        (*options)["max_timestep"].doc("Maximum internal timestep").withDefault(-1.));
    if (hmax > 0) {
      ropt[HMAX] = hmax;
    }
    /* ropt[HMIN]    : minimum absolute value of step size allowed.   *
     *                 Optional input. (Default is 0.0).              */
    const BoutReal hmin(
        (*options)["min_timestep"].doc("Minimum internal timestep").withDefault(-1.));
    if (hmin > 0) {
      ropt[HMIN] = hmin;
    }
    /* iopt[MAXORD] : maximum lmm order to be used by the solver.     *
     *                Optional input. (Default = 12 for ADAMS, 5 for  *
     *                BDF).                                           */
    const int maxOrder((*options)["max_order"].doc("Maximum order").withDefault(-1));
    if (maxOrder > 0) {
      iopt[MAXORD] = maxOrder;
    }
  }
  const bool use_adam((*options)["adams_moulton"]
                          .doc("Use Adams Moulton solver instead of BDF")
                          .withDefault(false));

  debug_on_failure =
      (*options)["debug_on_failure"]
          .doc("Run an aditional rhs if the solver fails with extra tracking")
          .withDefault(false);

  cvode_mem = CVodeMalloc(neq, solver_f, simtime, u, use_adam ? ADAMS : BDF, NEWTON, SS,
                          &reltol, &abstol, this, nullptr, optIn, iopt, ropt, machEnv);

  if (cvode_mem == nullptr) {
    throw BoutException("\tError: CVodeMalloc failed.\n");
  }

  /* Call CVSpgmr to specify the CVODE linear solver CVSPGMR with
     left preconditioning, modified Gram-Schmidt orthogonalization,
     default values for the maximum Krylov dimension maxl and the tolerance
     parameter delt, preconditioner setup and solve routines from the
     PVBBDPRE module, and the pointer to the preconditioner data block.    */

  if (use_precon) {
    CVSpgmr(cvode_mem, LEFT, SUN_MODIFIED_GS, precon_dimens, precon_tol, PVBBDPrecon,
            PVBBDPSol, pdata);
  } else {
    CVSpgmr(cvode_mem, NONE, SUN_MODIFIED_GS, 10, ZERO, PVBBDPrecon, PVBBDPSol, pdata);
  }

  // PvodeSolver is now initialised fully
  pvode_initialised = true;

  return (0);
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int PvodeSolver::run() {
  TRACE("PvodeSolver::run()");

  if (!pvode_initialised) {
    throw BoutException("PvodeSolver not initialised\n");
  }

  for (int i = 0; i < getNumberOutputSteps(); i++) {

    /// Run the solver for one output timestep
    simtime = run(simtime + getOutputTimestep());

    /// Check if the run succeeded
    if (simtime < 0.0) {
      // Step failed
      output.write("Timestep failed. Aborting\n");

      throw BoutException("PVODE timestep failed\n");
    }

    /// Call the monitor function

    if (call_monitors(simtime, i, getNumberOutputSteps())) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

BoutReal PvodeSolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run({})", tout);

  BoutReal* udata;

  // Set pointer to data array in vector u.
  udata = N_VDATA(u);

  // Run CVODE
  int flag;
  if (!monitor_timestep) {
    // Run in normal mode
    flag = CVode(cvode_mem, tout, u, &simtime, NORMAL);
  } else {
    // Run in single step mode, to call timestep monitors
    BoutReal internal_time = static_cast<CVodeMem>(cvode_mem)->cv_tn;
    //CvodeGetCurrentTime(cvode_mem, &internal_time);

    while (internal_time < tout) {
      // Run another step
      BoutReal last_time = internal_time;
      flag = CVode(cvode_mem, tout, u, &internal_time, ONE_STEP);
      if (flag < 0) {
        output_error.write("ERROR CVODE solve failed at t = {:e}, flag = {:d}\n",
                           internal_time, flag);
        return -1.0;
      }

      // Call timestep monitor
      call_timestep_monitors(internal_time, internal_time - last_time);
    }
    // Get output at the desired time
    flag = CVodeDky(cvode_mem, tout, 0, u);
    simtime = tout;
  }

  // Copy variables
  load_vars(udata);

  // Call rhs function to get extra variables at this time
  run_rhs(simtime);

  // Check return flag
  if (flag != SUCCESS) {
    output_error.write("ERROR CVODE step failed, flag = {:d}\n", flag);
    if (debug_on_failure) {
      CVodeMemRec* cv_mem = (CVodeMem)cvode_mem;
      if (f2d.empty() and v2d.empty() and v3d.empty()) {
        Options debug{};
        using namespace std::string_literals;
        Mesh* mesh{};
        for (const auto& prefix : {"pre_"s, "residuum_"s}) {
          std::vector<Field3D> list_of_fields{};
          std::vector<bool> evolve_bndrys{};
          for (const auto& f : f3d) {
            mesh = f.var->getMesh();
            Field3D to_load{0., mesh};
            to_load.allocate();
            to_load.setLocation(f.location);
            debug[fmt::format("{:s}{:s}", prefix, f.name)] = to_load;
            list_of_fields.push_back(to_load);
            evolve_bndrys.push_back(f.evolve_bndry);
          }
          pvode_load_data_f3d(evolve_bndrys, list_of_fields,
                              prefix == "pre_"s ? udata : N_VDATA(cv_mem->cv_acor));
        }

        for (auto& f : f3d) {
          f.F_var->enableTracking(fmt::format("ddt_{:s}", f.name), debug);
          setName(*f.var, f.name);
        }
        run_rhs(simtime);
        modelOutputVars(debug);

        for (auto& f : f3d) {
          debug[f.name] = *f.var;
	  if (f.var->hasParallelSlices()) {
	    saveParallel(debug, f.name, *f.var);
	  }
        }

        if (mesh != nullptr) {
          mesh->outputVars(debug);
          debug["BOUT_VERSION"].force(bout::version::as_double);
        }

        const std::string outname =
            fmt::format("{}/BOUT.debug.{}.nc",
                        Options::root()["datadir"].withDefault<std::string>("data"),
                        BoutComm::rank());

        bout::OptionsIO::create(outname)->write(debug);
        MPI_Barrier(BoutComm::get());
      } else {
        output_warn.write("debug_on_failure is currently only supported for Field3Ds");
      }
    }
    return (-1.0);
  }

  return simtime;
}

/**************************************************************************
 * RHS function
 **************************************************************************/

void PvodeSolver::rhs(int UNUSED(N), BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: PvodeSolver::rhs({})", t);

  // Get current timestep
  hcur = 0.0; //((CVodeMemRec*) cvode_mem)->cv_h;

  // Load state from CVODE
  load_vars(udata);

  // Call function
  run_rhs(t);

  // Save derivatives to CVODE
  save_derivs(dudata);
}

void PvodeSolver::gloc(int UNUSED(N), BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: PvodeSolver::gloc({})", t);

  Timer timer("rhs");

  // Load state from CVODE
  load_vars(udata);

  // Call function
  run_rhs(t);

  // Save derivatives to CVODE
  save_derivs(dudata);
}

/**************************************************************************
 * CVODE rhs function
 **************************************************************************/

void solver_f(integer N, BoutReal t, N_Vector u, N_Vector udot, void* f_data) {
  BoutReal *udata, *dudata;
  PvodeSolver* s;

  udata = N_VDATA(u);
  dudata = N_VDATA(udot);

  s = static_cast<PvodeSolver*>(f_data);

  s->rhs(N, t, udata, dudata);
}

// Preconditioner RHS
void solver_gloc(integer N, BoutReal t, BoutReal* u, BoutReal* udot, void* f_data) {
  PvodeSolver* s;

  s = static_cast<PvodeSolver*>(f_data);

  s->gloc(N, t, u, udot);
}

// Preconditioner communication function
void solver_cfn(integer UNUSED(N), BoutReal UNUSED(t), N_Vector UNUSED(u),
                void* UNUSED(f_data)) {
  // doesn't do anything at the moment
}

#endif
