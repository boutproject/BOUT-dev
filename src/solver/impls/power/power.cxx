#include "mpi.h"

#include "power.hxx"

#include <boutcomm.hxx>
#include <bout/sys/timer.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include <output.hxx>

int PowerSolver::init(int nout, BoutReal tstep) {
  TRACE("Initialising Power solver");
  
  /// Call the generic initialisation first
  if(Solver::init(nout, tstep))
    return 1;
  
  output << "\n\tPower eigenvalue solver\n";
  
  nsteps = nout; // Save number of output steps
  
  // Get options
  OPTION(options, curtime, 0.0);

  // Calculate number of variables
  nlocal = getLocalN();
  
  // Get total problem size
  if (bout::globals::mpi->MPI_Allreduce(&nlocal, &nglobal, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed in EulerSolver::init");
  }
  
  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n",
	       n3Dvars(), n2Dvars(), nglobal, nlocal);
  
  // Allocate memory
  f0.reallocate(nlocal);

  eigenvalue = 0.0;
  
  // Put starting values into f0
  save_vars(std::begin(f0));

  return 0;
}

int PowerSolver::run() {
  TRACE("PowerSolver::run()");
  
  // Make sure that f0 has a norm of 1
  divide(f0, norm(f0));
  
  for(int s=0;s<nsteps;s++) {

    load_vars(std::begin(f0));
    run_rhs(curtime);
    save_derivs(std::begin(f0));

    // Estimate eigenvalue
    eigenvalue = norm(f0);
    
    // Normalise
    divide(f0, eigenvalue);
    
    /// Call the monitor function. The eigenvalue
    /// is given rather than time, so it appears
    /// in the output logs
    if(call_monitors(eigenvalue, s, nsteps)) {
      // User signalled to quit
      
      output.write("Monitor signalled to quit. Returning\n");
      break;
    }
  }
  
  return 0;
}

void PowerSolver::outputVars(Options& output_options, bool save_repeat) {
  Timer time("io");
  // Include base class functionality
  this->Solver::outputVars(output_options, save_repeat);

  // Save the eigenvalue to the output
  output_options["eigenvalue"].assignRepeat(eigenvalue, "t", save_repeat, "Solver");
}

BoutReal PowerSolver::norm(Array<BoutReal> &state) {
  BoutReal total = 0.0, result;
  
  for(int i=0;i<nlocal;i++)
    total += state[i]*state[i];

  total /= static_cast<BoutReal>(nglobal);

  bout::globals::mpi->MPI_Allreduce(&total, &result, 1, MPI_DOUBLE, MPI_SUM,
                                    BoutComm::get());

  return sqrt(result);
}

void PowerSolver::divide(Array<BoutReal> &in, BoutReal value) {
  for(int i=0;i<nlocal;i++)
    in[i] /= value;
}
