
#include "euler.hxx"

#include <boutcomm.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <bout/openmpwrap.hxx>

#include <cmath>

#include <output.hxx>

void EulerSolver::setMaxTimestep(BoutReal dt) {
  if(dt >= cfl_factor*timestep)
    return; // Already less than this
  
  timestep = dt*0.99 / cfl_factor; // Slightly below to avoid re-setting to same value over again
  timestep_reduced = true;
}

int EulerSolver::init(int nout, BoutReal tstep) {
  TRACE("Initialising Euler solver");
  
  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;
  
  output << "\n\tEuler solver\n";

  nsteps = nout; // Save number of output steps
  out_timestep = tstep;
  
  // Get options
  OPTION(options, timestep, tstep);
  OPTION(options, mxstep, 500); // Maximum number of steps between outputs
  OPTION(options, cfl_factor, 2.);

  // Calculate number of variables
  nlocal = getLocalN();
  
  // Get total problem size
  int neq;
  if(bout::globals::mpi->MPI_Allreduce(&nlocal, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed in EulerSolver::init");
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3Dvars(), n2Dvars(), neq, nlocal);
  
  // Allocate memory
  f0.reallocate(nlocal);
  f1.reallocate(nlocal);

  // Put starting values into f0
  save_vars(std::begin(f0));

  return 0;
}

int EulerSolver::run() {
  TRACE("EulerSolver::run()");
  
  for(int s=0;s<nsteps;s++) {
    BoutReal target = simtime + out_timestep;
    
    bool running = true;
    int internal_steps = 0;
    do {
      // Take a step
      BoutReal dt_limit = timestep; // Store the timestep
      
      if((simtime + timestep) >= target) {
        // Make sure the last timestep is on the output 
        timestep = target - simtime; 
        running = false;
      }
      
      BoutReal old_timestep = timestep;
      
      timestep_reduced = false;
      take_step(simtime, timestep, f0, f1);
      
      // Check with all processors if timestep was reduced
      
      BoutReal newdt_local = 10.*old_timestep; // Signal no change
      if(timestep_reduced)
        newdt_local = timestep;
      
      BoutReal newdt;
      if(bout::globals::mpi->MPI_Allreduce(&newdt_local, &newdt, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get())) {
        throw BoutException("MPI_Allreduce failed in EulerSolver::run");
      }

      // If timestep_reduced re-run
      if(newdt < old_timestep) { // At least one processor reduced the timestep
        timestep = newdt;
        take_step(simtime, timestep, f0, f1);
        dt_limit = timestep; // This becomes the new limit
        running = true; // Need another step
      }
      
      // Taken a step, swap buffers
      swap(f1, f0);
      simtime += timestep;

      internal_steps++;
      if(internal_steps > mxstep)
        throw BoutException("ERROR: MXSTEP exceeded. simtime=%e, timestep = %e\n", simtime, timestep);
      
      // Call timestep monitors
      call_timestep_monitors(simtime, timestep);
      
      timestep = dt_limit; // Change back to limiting timestep
    }while(running);

    load_vars(std::begin(f0)); // Put result into variables
    // Call rhs function to get extra variables at this time
    run_rhs(simtime);
    
    iteration++; // Advance iteration number
    
    /// Call the monitor function
    
    if(call_monitors(simtime, s, nsteps)) {
      // Stop simulation
      break;
    }
  }
  
  return 0;
}

void EulerSolver::take_step(BoutReal curtime, BoutReal dt, Array<BoutReal> &start,
                            Array<BoutReal> &result) {

  load_vars(std::begin(start));
  run_rhs(curtime);
  save_derivs(std::begin(result));

  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    result[i] = start[i] + dt*result[i];
}
