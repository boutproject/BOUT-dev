#include "mpi.h"

#include "power.hxx"

#include <boutcomm.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include <output.hxx>

int PowerSolver::init(bool restarting, int nout, BoutReal tstep) {
  int msg_point = msg_stack.push("Initialising Power solver");

  /// Call the generic initialisation first
  if(Solver::init(restarting, nout, tstep))
    return 1;

  output << "\n\tPower eigenvalue solver\n";

  nsteps = nout; // Save number of output steps

  // Get options
  OPTION(options, curtime, 0.0);

  // Calculate number of variables
  nlocal = getLocalN();

  // Get total problem size
  if(MPI_Allreduce(&nlocal, &nglobal, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed in EulerSolver::init");
  }

  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
               n3Dvars(), n2Dvars(), nglobal, nlocal);

  // Allocate memory
  f0 = new BoutReal[nlocal];

  // Save the eigenvalue to the output
  dump.add(eigenvalue, "eigenvalue", 1);
  eigenvalue = 0.0;

  // Put starting values into f0
  save_vars(f0);

  msg_stack.pop(msg_point);

  return 0;
}

int PowerSolver::run() {
  int msg_point = msg_stack.push("PowerSolver::run()");

  // Make sure that f0 has a norm of 1
  divide(f0, norm(f0));

  for(int s=0;s<nsteps;s++) {

    load_vars(f0);
    run_rhs(curtime);
    save_derivs(f0);

    // Estimate eigenvalue
    eigenvalue = norm(f0);

    // Normalise
    divide(f0, eigenvalue);

    /// Write the restart file
    restart.write();

    if((archive_restart > 0) && (iteration % archive_restart == 0)) {
      restart.write("%s/BOUT.restart_%04d.%d.%s", restartdir.c_str(), iteration, MYPE, restartext.c_str());
    }

    /// Call the monitor function. The eigenvalue
    /// is given rather than time, so it appears
    /// in the output logs
    if(call_monitors(eigenvalue, s, nsteps)) {
      // User signalled to quit

      // Write restart to a different file
      restart.write("%s/BOUT.final.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());

      output.write("Monitor signalled to quit. Returning\n");
      break;
    }

    // Reset iteration and wall-time count
    rhs_ncalls = 0;
  }

  msg_stack.pop(msg_point);

  return 0;
}

BoutReal PowerSolver::norm(BoutReal *state) {
  BoutReal total = 0.0, result;

  for(int i=0;i<nlocal;i++)
    total += state[i]*state[i];

  total /= (BoutReal) nglobal;

  MPI_Allreduce(&total, &result, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get());

  return sqrt(result);
}

void PowerSolver::divide(BoutReal *in, BoutReal value) {
  for(int i=0;i<nlocal;i++)
    in[i] /= value;
}
