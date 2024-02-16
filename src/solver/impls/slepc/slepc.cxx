#/**************************************************************************
 * Interface to SLEPc solver
 *
 **************************************************************************
 * Copyright 2014 B.D.Dudson, D. Dickinson
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

#include "bout/build_config.hxx"

#if BOUT_HAS_SLEPC

#include "slepc.hxx"
#include <bout/boutcomm.hxx>
#include <bout/globals.hxx>
#include <bout/interpolation.hxx>
#include <bout/msg_stack.hxx>
#include <bout/output.hxx>

#include <cstdlib>

#if BOUT_USE_SIGFPE
#include <fenv.h>
#endif

namespace {
/// Disable floating-point exceptions in a scope, reenable them on exit
struct QuietFPE {
#if BOUT_USE_SIGFPE
  QuietFPE() : flags(fegetexcept()) { fedisableexcept(flags); }
  ~QuietFPE() { feenableexcept(flags); }

private:
  int flags;
#endif
};
} // namespace

std::string formatEig(BoutReal reEig, BoutReal imEig);

// The callback function for the shell matrix-multiply operation
// A simple wrapper around the SlepcSolver advanceStep routine
PetscErrorCode advanceStepWrapper(Mat matOperator, Vec inData, Vec outData) {
  PetscFunctionBegin;
  SlepcSolver* ctx;
  // Here we set the ctx pointer to the solver instance
  MatShellGetContext(matOperator, reinterpret_cast<void**>(&ctx));
  // Actually advance
  PetscFunctionReturn(ctx->advanceStep(matOperator, inData, outData));
}

// The callback function for the eigenvalue comparison
// A simple wrapper around the SlepcSolver compareEigs routine
PetscErrorCode compareEigsWrapper(PetscScalar ar, PetscScalar ai, PetscScalar br,
                                  PetscScalar bi, PetscInt* res, void* ctx) {
  PetscFunctionBegin;

  // Cast context as SlepcSolver and call the actual compare routine
  SlepcSolver* myCtx;
  myCtx = static_cast<SlepcSolver*>(ctx);
  myCtx->compareState = myCtx->compareEigs(ar, ai, br, bi);

  *res = myCtx->compareState;
  PetscFunctionReturn(0);
}

// The callback function for the monitor
// A simple wrapper around the SlepcSolver compareEigs routine
PetscErrorCode monitorWrapper(EPS UNUSED(eps), PetscInt its, PetscInt nconv,
                              PetscScalar* eigr, PetscScalar* eigi, PetscReal* errest,
                              PetscInt nest, void* mctx) {
  PetscFunctionBegin;
  // Cast context as SlepcSolver and call the actual compare routine
  SlepcSolver* myCtx;
  myCtx = static_cast<SlepcSolver*>(mctx);
  myCtx->monitor(its, nconv, eigr, eigi, errest, nest);
  PetscFunctionReturn(0);
}

// The callback function for applying the shell spectral transformation
PetscErrorCode stApplyWrapper(ST st, Vec vecIn, Vec vecOut) {
  PetscFunctionBegin;

  // First get the context of the st object, cast to correct type
  SlepcSolver* myCtx;
  STShellGetContext(st, reinterpret_cast<void**>(&myCtx));

  // Do the matrix vector multiply -- same as STSHIFT with zero shift
  // Use the advanceStepWrapper so any mods made in advanceStep are
  // taken into account
  advanceStepWrapper(myCtx->shellMat, vecIn, vecOut);

  // Note an alternative  approach, which would allow shellMat to remain
  // private, would be to extract the ST object during createEPS and set
  // the operator their. We could then use STGetOperators to obtain a
  // reference to the shell mat.
  PetscFunctionReturn(0);
}

// The callback function for transforming the eigenvalues in the
// custom shell spectral transformation
PetscErrorCode stBackTransformWrapper(ST st, PetscInt nEig, PetscScalar* eigr,
                                      PetscScalar* eigi) {
  PetscFunctionBegin;
  // First get the context of the st object and cast to correct type
  SlepcSolver* myCtx;
  STShellGetContext(st, reinterpret_cast<void**>(&myCtx));

  // Convert to bout eigenvalue
  BoutReal tmpR, tmpI;
  for (PetscInt iEig = 0; iEig < nEig; iEig++) {
    myCtx->slepcToBout(eigr[iEig], eigi[iEig], tmpR, tmpI, true);
    eigr[iEig] = tmpR;
    eigi[iEig] = tmpI;
  };
  PetscFunctionReturn(0);
}

// Helper function
std::string formatEig(BoutReal reEig, BoutReal imEig) {
  const std::string rePad = (reEig < 0) ? "-" : " ";
  const std::string imPad = (imEig < 0) ? "-" : "+";

  std::stringstream tmp;
  tmp.precision(5);
  tmp << std::scientific;
  // Note we use abs here and put the -/+ into pads to cope with
  // the case where the Eig is -0.000/0.000 which require different
  // padding but evaluate as the same number when doing comparisons
  tmp << rePad << std::abs(reEig) << imPad << std::abs(imEig) << "i";
  return tmp.str();
}

SlepcSolver::SlepcSolver(Options* options) {
  has_constraints = false;
  initialised = false;
  stIsShell = PETSC_FALSE;

  // Slepc settings in the Solver section
  auto& options_ref = *options;

  nEig = options_ref["nEig"]
             .doc("Number of eigenvalues to compute. 0 means keep the current value, "
                  "i.e. autoset")
             .withDefault(0);

  tol = options_ref["tol"].doc("SLEPc tolerance").withDefault(1.0e-6);
  maxIt = options_ref["maxIt"]
              .doc("Maximum iterations")
              .withDefault(static_cast<int>(PETSC_DEFAULT));

  mpd = options_ref["mpd"]
            .doc("Maximum dimension allowed for the projected problem")
            .withDefault(static_cast<int>(PETSC_DEFAULT));

  ddtMode = options_ref["ddtMode"].withDefault(true);

  targRe = options_ref["targRe"]
               .doc("Target frequency when using user eig comparison")
               .withDefault(0.0);
  targIm = options_ref["targIm"]
               .doc("Target growth rate when using user eig comparison")
               .withDefault(0.0);

  // Convert bout targs to slepc
  bool userWhichDefault = false;
  if (targRe == 0.0 && targIm == 0.0) {
    target = 999.0;
  } else {
    PetscScalar slepcRe, slepcIm;
    boutToSlepc(targRe, targIm, slepcRe, slepcIm);
    dcomplex tmp(slepcRe, slepcIm);
    target = std::abs(tmp);

    // If we've set a target then we change the default
    // for the userWhich variable as targets work best with this
    // Note this means we only use target in the case where we
    // specify targRe/targIm *and* explicitly set userWhich=false
    userWhichDefault = true;
  }

  target = options_ref["target"]
               .doc("If 999 we don't set the target. This is SLEPc eig target")
               .withDefault(target);

  userWhich = options_ref["userWhich"].withDefault(userWhichDefault);

  // Generic settings
  useInitial = options_ref["useInitial"].withDefault(!ddtMode);
  debugMonitor = options_ref["debugMonitor"].withDefault(false);

  selfSolve = options_ref["selfSolve"]
                  .doc("Solver to advance the state of the system")
                  .withDefault(false);

  if (ddtMode && !selfSolve) {
    // We need to ensure this so that we don't try to use
    // advanceSolver elsewhere. The other option would be to
    // create advanceSolver below in ddtMode but we just don't
    // use it.
    output << "Overridding selfSolve as ddtMode = true\n";
    selfSolve = true;
  }
  eigenValOnly = options_ref["eigenValOnly"].withDefault(false);

  if (!selfSolve && !ddtMode) {
    // Use a sub-section called "advance"
    advanceSolver = SolverFactory::getInstance().create(options->getSection("advance"));
  }
}

SlepcSolver::~SlepcSolver() {
  if (initialised) {
    // Free memory
    if (eps) {
      EPSDestroy(&eps);
    };
    if (shellMat) {
      MatDestroy(&shellMat);
    };
    initialised = false;
  }
}

int SlepcSolver::init() {

  TRACE("Initialising SLEPc solver");

  // Report initialisation
  output.write("Initialising SLEPc solver\n");
  if (selfSolve) {
    Solver::init();

    // If no advanceSolver then can only advance one step at a time
    setNumberOutputSteps(1);
  }

  // Read options
  comm = PETSC_COMM_WORLD;

  // Initialise advanceSolver if not self
  if (!selfSolve && !ddtMode) {
    advanceSolver->init();
  }

  // Calculate grid sizes
  localSize = getLocalN();

  // Also create vector for derivs etc. if SLEPc in charge of solving
  if (selfSolve && !ddtMode) {
    // Allocate memory
    f0.reallocate(localSize);
    f1.reallocate(localSize);
  }

  // Get total problem size
  int neq;
  if (bout::globals::mpi->MPI_Allreduce(&localSize, &neq, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed in SlepcSolver::init");
  }

  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n", n3Dvars(),
               n2Dvars(), neq, localSize);

  // Create EPS solver
  createEPS();

  // Return ok
  return 0;
}

int SlepcSolver::run() {
  // Now the basic idea with slepc is that:
  // Whilst n_eig_converged<n_eig_desired do
  // 1. Let slepc set the initial fields
  // 2. Use the advanceSolver to evolve fields over a certain time
  // 3. Package up the evolved fields for slepc to analyse
  // Once converged data found:
  // 1. Get converged eigenvalues/eigenvectors
  // 2. Write to file
  // 3. Clean up
  //--> The first section is handled by calling EPSSolve(eps) with appropriate shellMat
  //--> The second section has to be handled by this solver

  // Find the eigenvalues
  EPSSolve(eps);

  // Analyse and dump to file
  if (!eigenValOnly) {
    analyseResults();
  }
  return 0;
}

// This routine takes a Vec type object of length localSize and
// unpacks it into the local fields
void SlepcSolver::vecToFields(Vec& inVec) {
  /*
    Whilst this routine does indeed populate the field variables
    most (/all?) solvers overwrite this data on call to run() as
    they actually evolve a single vector which they unpack (not
    unlike this routine). The typical run work flow is something
    like:
      1. Unpack internal vector A into fields.
      2. Use phys_run to calculate ddt of fields.
      3. Pack ddt into an internal vector B.
      4. Use B to evolve A.
    Unfortunately each solver implementation uses a different
    internal vector approach which means it's a implementation
    specific thing which makes it hard/impossible to actually
    set from here. Possible solutions are:
      1. Add a virtual function to solver, which each implementation
         overrides, which just copies field data into appropriate
         internal vector or copies out of internal vector into
         fields (for fieldsToVec).
      2. ?
   */

  // Get pointer to data
  const BoutReal* point;
  VecGetArrayRead(inVec, &point);

  // Copy data from point into fields
  load_vars(const_cast<BoutReal*>(point));
  // Note as the solver instances only have pointers to the
  // fields we can use the SlepcSolver load_vars even if we're
  // not using selfSolve=True

  if (!selfSolve && !ddtMode) {
    // Solver class used must support this procedure which resets any internal state
    // data such that it now holds the same data as the fields
    advanceSolver->resetInternalFields();
  }

  // Restore array
  VecRestoreArrayRead(inVec, &point);
}

// This routine packs the local fields into a vector
void SlepcSolver::fieldsToVec(Vec& outVec) {
  // Get pointer to data
  PetscScalar* point;
  VecGetArray(outVec, &point);

  // Copy fields into point
  if (!ddtMode) {
    save_vars(point);
  } else {
    save_derivs(point);
  };
  // Note as the solver instances only have pointers to the
  // fields we can use the SlepcSolver save_vars even if we're
  // not using selfSolve=True

  // Restore array
  VecRestoreArray(outVec, &point);
}

// Create a shell matrix operator
void SlepcSolver::createShellMat() {
  output << "Creating shellMat with local size : " << localSize << "\n";

  // Create the shell matrix
  // Note we pass the this reference as the matrix context.
  // This allows us to access the SlepcSolver internals from within
  // routines called directly by Petsc/Slepc (i.e. our wrapper functions)
  MatCreateShell(comm, localSize, localSize, PETSC_DETERMINE, PETSC_DETERMINE, this,
                 &shellMat);
  // Define the mat_mult operation --> Define what routine returns M.x, where M
  // is the time advance operator and x are the initial field conditions
  MatShellSetOperation(shellMat, MATOP_MULT,
                       reinterpret_cast<void (*)()>(&advanceStepWrapper));

  // The above function callback can cause issues as member functions have a hidden "this"
  // argument which means if Slepc calls this->advanceStep(Mat,Vec,Vec) this is actually
  // this->advanceStep(this,Mat,Vec,Vec) meaing "this" gets redefined to Mat,
  // and the two Vecs are mangled, this messes up memory and the code crashes.
  // Alternatives include:
  //  1. Making advanceStep a static function
  //  2. Making advanceStep a non-member function
  // These alternatives generally divorce the advanceStep from the SlepcSolver class
  // which might make it difficult to access required data.
  // We've therefore gone with a third option where we define the callback to be a
  // non-member function which then gets the context pointer attached to shellMat
  // which points to the SlepcSolver instance. This can then be used to call the
  // advanceStep member function as if it were called within "this".
}

// Create an EPS Solver
void SlepcSolver::createEPS() {
  // First need to create shell matrix
  createShellMat();

  // Now construct EPS
  EPSCreate(comm, &eps);
  EPSSetOperators(eps, shellMat, nullptr);
  EPSSetProblemType(eps, EPS_NHEP); // Non-hermitian

  // Probably want to read options and set EPS properties
  // at this point.
  EPSSetDimensions(eps, nEig, PETSC_DECIDE, mpd);
  EPSSetTolerances(eps, tol, maxIt);
  if (!(target == 999)) {
    EPSSetTarget(eps, target);
  }

  // Set the user comparison function
  if (userWhich) {
    EPSSetEigenvalueComparison(eps, compareEigsWrapper, this);
    EPSSetWhichEigenpairs(eps, EPS_WHICH_USER);
  }

  // Update options from command line
  EPSSetFromOptions(eps);

  // Register a monitor
  EPSMonitorSet(eps, &monitorWrapper, this, nullptr);

  // Initialize shell spectral transformation if selected by user
  // Note currently the only way to select this is with the
  //"-st_type shell" command line option, should really add a
  // BOUT input flag to force it
  EPSGetST(eps, &st);
  PetscObjectTypeCompare(reinterpret_cast<PetscObject>(st), STSHELL, &stIsShell);
  if (stIsShell) {
    // Set the user-defined routine for applying the operator
    STShellSetApply(st, &stApplyWrapper);

    // Set the STShell context to be the slepcSolver so we can access
    // the solver internals from within the spectral transform routines
    STShellSetContext(st, this);

    // Set the routine to transform the eigenvalues back
    STShellSetBackTransform(st, stBackTransformWrapper);

    // Define the transformations name (optional)
    PetscObjectSetName(reinterpret_cast<PetscObject>(st), "Exponential Linear ST");
  };

  // Should probably call a routine here which interrogates eps
  // to determine the important settings that have been used and dump
  // the settings to screen/file/dmp?
  // I think there may be a Slepc flag which will do this (to screen)
  // but not sure if we can force this is the code (without messing with argv).

  // Set initial space i.e. first guess
  if (useInitial) { // Doesn't seem to help the ddtMode very much so recommend off
    Vec initVec, rightVec;
    bool ddtModeBackup = ddtMode;

#if PETSC_VERSION_LT(3, 6, 0)
    MatGetVecs(shellMat, &rightVec, &initVec);
#else
    MatCreateVecs(shellMat, &rightVec, &initVec);
#endif
    ddtMode = false; // Temporarily disable as initial ddt values not set
    fieldsToVec(initVec);
    ddtMode = ddtModeBackup; // Restore state
    EPSSetInitialSpace(eps, 1, &initVec);
    VecDestroy(&initVec);
    VecDestroy(&rightVec);
  };
}

// This routine takes initial conditions provided by SLEPc, uses this to set the fields,
// advances them with the attached solver and then returns the evolved fields in a slepc
// structure.
// Note: Hidden "this" argument prevents Slepc calling this routine directly
int SlepcSolver::advanceStep(Mat& UNUSED(matOperator), Vec& inData, Vec& outData) {

  // First unpack input into fields
  vecToFields(inData);

  // Now advance
  int retVal;

  if (ddtMode) {
    // In ddtMode we just want the time derivative of the fields
    retVal = run_rhs(0.0);
  } else {
    // Here we actually advance one (big) step
    if (selfSolve) {
      // If we don't have an external solver then we have to advance the solution
      // ourself. This is currently done using Euler and only advances by a small step
      // Not recommended!
      retVal = run_rhs(0.0);
      // Here we add dt*ddt(Fields) to fields to advance solution (cf. Euler)
      save_vars(std::begin(f0));
      save_derivs(std::begin(f1));
      for (int iVec = 0; iVec < localSize; iVec++) {
        f0[iVec] += f1[iVec] * getOutputTimestep();
      }
      load_vars(std::begin(f0));
    } else {
      // Here we exploit one of the built in solver implementations to advance the
      // prescribed fields by a big (tstep*nstep) step.
      retVal = advanceSolver->run();
    }
  }

  // Now pack evolved fields into output
  fieldsToVec(outData);

  // Return
  return retVal;
}

// This routine can be used by Slepc to decide which of two eigenvalues is "preferred"
// Allows us to look at real and imaginary components seperately which is not possible
// with Slepc built in comparisons when Slepc is compiled without native complex support
// (required) Note must be wrapped by non-member function to be called by Slepc
int SlepcSolver::compareEigs(PetscScalar ar, PetscScalar ai, PetscScalar br,
                             PetscScalar bi) {
  BoutReal arBout, aiBout, brBout, biBout;

  // First convert to BOUT values
  slepcToBout(ar, ai, arBout, aiBout);
  slepcToBout(br, bi, brBout, biBout);

  // Now we calculate the distance between eigenvalues and target.
  const auto da = sqrt(pow(arBout - targRe, 2) + pow(aiBout - targIm, 2));
  const auto db = sqrt(pow(brBout - targRe, 2) + pow(biBout - targIm, 2));

  // Now we decide which eigenvalue is preferred.
  int retVal;

  // Smallest distance from complex target
  // If prefer B we return +ve
  if (da > db) {
    retVal = 1;
    // If prefer A we return -ve
  } else if (db > da) {
    retVal = -1;
    // If we don't prefer either we return 0
  } else {
    retVal = 0;
  };

  return retVal;
}

// This is an example of a custom monitor which Slepc can call (not directly) to report
// the current status of the run. Note we could see how many new eigenpairs have been
// found since last called and then write their data to file so that we get progressive
// output rather than waiting until the end to write everything. Unfortunately it seems
// that currently SLEPc does not support using EPSGetEigenvector or EPSGetEigenpair before
// EPSSolve has finished. As such it's not possible to write out the eigenvectors from
// this monitor routine. It should still be possible to write the eigenvalues here, but to
// then get the eigenvectors at the correct time indices later would require resetting the
// time index. I'm not sure if the Datafile object supports this. Note must be wrapped by
// non-member function to be called by Slepc
void SlepcSolver::monitor(PetscInt its, PetscInt nconv, PetscScalar eigr[],
                          PetscScalar eigi[], PetscReal errest[], PetscInt UNUSED(nest)) {
  static int nConvPrev = 0;

  // Disable floating-point exceptions for the duration of this function
  [[maybe_unused]] QuietFPE quiet_fpe{};

  // No output until after first iteration
  if (its < 1) {
    return;
  }

  static bool first = true;
  if (eigenValOnly && first) {
    first = false;
    resetIterationCounter();
  }

  // Temporary eigenvalues, converted from the SLEPc eigenvalues
  BoutReal reEigBout, imEigBout;

  // Only report unconverged eigenvalues if we don't have all the requested ones
  if (nconv < nEig) {
    slepcToBout(eigr[nconv], eigi[nconv], reEigBout, imEigBout);

    // This line more or less replicates the normal slepc output (when using -eps_monitor)
    // but reports Bout eigenvalues rather than the Slepc values. Note we haven't changed
    // error estimate.
    output.write(" {} nconv={}\t first unconverged value (error) {}\t({})\n", its, nconv,
                 formatEig(reEigBout, imEigBout), errest[nconv]);
  }

  // The following can be quite noisy so may want to add a flag to disable/enable.
  const int newConv = nconv - nConvPrev;
  if (newConv > 0) {
    output.write("Found {} new converged eigenvalues:\n", newConv);
    for (PetscInt i = nConvPrev; i < nconv; i++) {
      slepcToBout(eigr[i], eigi[i], reEigBout, imEigBout);
      output.write("\t{}\t: {} --> {}\n", i, formatEig(eigr[i], eigi[i]),
                   formatEig(reEigBout, imEigBout));
      if (eigenValOnly) {
        // Silence the default monitor
        WithQuietOutput progress{output_progress};
        // Call monitors so fields get written
        call_monitors(reEigBout, incrementIterationCounter(), getNumberOutputSteps());
        call_monitors(imEigBout, incrementIterationCounter(), getNumberOutputSteps());
      }
    }
  }

  // Update the number of converged modes already investigated.
  nConvPrev = nconv;
}

// Convert a slepc eigenvalue to a BOUT one
void SlepcSolver::slepcToBout(PetscScalar& reEigIn, PetscScalar& imEigIn,
                              BoutReal& reEigOut, BoutReal& imEigOut, bool force) {
  // If not stIsShell then the slepc eigenvalue is actually
  // Exp(-i*Eig_Bout*tstep) for ddtMode = false
  //-i*Eig_Bout for ddtMode = true
  // where Eig_Bout is the actual eigenvalue and tstep is the time step
  // the solution is evolved over.
  // This routine returns Eig_Bout

  // The optional input force is used by the shell spectral transform
  // in the back transform to force a conversion, which allows us to
  // otherwise skip any transformation when in shell mode (i.e. the back
  // transform routine is the only place we deal with the raw slepc eigenvalue
  // in shell ST mode).

  // If shellST and not forcing we just set the input and output eigenvalues
  // equal and return.
  if (stIsShell && !force) {
    reEigOut = reEigIn;
    imEigOut = imEigIn;
    return;
  }

  const dcomplex slepcEig(reEigIn, imEigIn);
  const dcomplex ci(0.0, 1.0);

  // Protect against the 0,0 trivial eigenvalue
  if (ddtMode and std::abs(slepcEig) < 1.0e-10) {
    reEigOut = 0.0;
    imEigOut = 0.0;
    return;
  }

  const dcomplex boutEig =
      ddtMode ? slepcEig * ci
              : ci * log(slepcEig) / (getOutputTimestep() * getNumberOutputSteps());

  // Set return values
  reEigOut = boutEig.real();
  imEigOut = boutEig.imag();
}

// Convert a BOUT++ eigenvalue to a Slepc one
void SlepcSolver::boutToSlepc(BoutReal& reEigIn, BoutReal& imEigIn, PetscScalar& reEigOut,
                              PetscScalar& imEigOut, bool force) {

  // If shellST and not forcing we just set the input and output eigenvalues
  // equal and return.
  if (stIsShell && !force) {
    reEigOut = reEigIn;
    imEigOut = imEigIn;
    return;
  }

  const dcomplex boutEig(reEigIn, imEigIn);
  const dcomplex ci(0.0, 1.0);
  const dcomplex slepcEig =
      ddtMode ? -ci * boutEig
              : exp(-ci * boutEig * (getOutputTimestep() * getNumberOutputSteps()));

  // Set return values
  reEigOut = slepcEig.real();
  imEigOut = slepcEig.imag();
}

// Interrogate eps to find out how many eigenvalues we've found etc.
void SlepcSolver::analyseResults() {
  PetscInt nEigFound;

  // Find how many eigenvalues have been found
  EPSGetConverged(eps, &nEigFound);

  if (nEigFound < 1) {
    output << "Warning : No converged eigenvalues found!\n";
    return;
  }

  // Now loop over each converged eigenpair and output eigenvalue

  output << "Converged eigenvalues :\n"
            "\tIndex\tSlepc eig (mag.)\t\t\tBOUT eig (mag.)\n";

  resetIterationCounter();

  // Declare and create vectors to store eigenfunctions
  Vec vecReal, vecImag;
#if PETSC_VERSION_LT(3, 6, 0)
  MatGetVecs(shellMat, &vecReal, &vecImag);
#else
  MatCreateVecs(shellMat, &vecReal, &vecImag);
#endif

  for (PetscInt iEig = 0; iEig < nEigFound; iEig++) {
    // Get slepc eigenvalue
    PetscScalar reEig, imEig;
    EPSGetEigenvalue(eps, iEig, &reEig, &imEig);
    dcomplex slepcEig(reEig, imEig);

    // Report
    output << "\t" << iEig << "\t" << formatEig(reEig, imEig) << "\t("
           << std::abs(slepcEig) << ")";

    // Get BOUT eigenvalue
    BoutReal reEigBout, imEigBout;
    slepcToBout(reEig, imEig, reEigBout, imEigBout);
    dcomplex boutEig(reEigBout, imEigBout);

    // Report
    output << "\t" << formatEig(reEigBout, imEigBout) << "\t(" << std::abs(boutEig)
           << ")\n";

    // Get eigenvector
    EPSGetEigenvector(eps, iEig, vecReal, vecImag);

    // Write real part of eigen data
    // First dump real part to fields
    vecToFields(vecReal);

    // Run the rhs in order to calculate aux fields
    run_rhs(0.0);

    // Silence the default monitor
    WithQuietOutput progress{output_progress};
    // Call monitors so fields get written
    call_monitors(reEigBout, incrementIterationCounter(), getNumberOutputSteps());

    // Now write imaginary part of eigen data
    // First dump imag part to fields
    vecToFields(vecImag);
    call_monitors(imEigBout, incrementIterationCounter(), getNumberOutputSteps());
  }

  // Destroy vectors
  VecDestroy(&vecReal);
  VecDestroy(&vecImag);
}

#endif // BOUT_HAS_SLEPC
