
#include "euler.hxx"

#include <boutcomm.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include <output.hxx>

EulerSolver::EulerSolver() : Solver() {
  
}

EulerSolver::~EulerSolver() {

}

void EulerSolver::setMaxTimestep(BoutReal dt) {
  if(dt >= cfl_factor*timestep)
    return; // Already less than this
  
  timestep = dt*0.99 / cfl_factor; // Slightly below to avoid re-setting to same value over again
  timestep_reduced = true;
}

int EulerSolver::init(bool restarting, int nout, BoutReal tstep) {
  int msg_point = msg_stack.push("Initialising Euler solver");
  
  /// Call the generic initialisation first
  if(Solver::init(restarting, nout, tstep))
    return 1;
  
  output << "\n\tEuler solver\n";

  nsteps = nout; // Save number of output steps
  out_timestep = tstep;
  
  // Get options
  Options *options = Options::getRoot();
  options = options->getSection("solver");
  OPTION(options, start_timestep, tstep);
  OPTION(options, mxstep, 500); // Maximum number of steps between outputs
  OPTION(options, cfl_factor, 2.);

  // Calculate number of variables
  nlocal = getLocalN();
  
  // Get total problem size
  int neq;
  if(MPI_Allreduce(&nlocal, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed in EulerSolver::init");
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3Dvars(), n2Dvars(), neq, nlocal);
  
  // Allocate memory
  f0 = new BoutReal[nlocal];
  f1 = new BoutReal[nlocal];
  
  // Put starting values into f0
  save_vars(f0);
  
  msg_stack.pop(msg_point);

  return 0;
}

int EulerSolver::run() {
  int msg_point = msg_stack.push("EulerSolver::run()");
  
  timestep = start_timestep;
  
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
      if(MPI_Allreduce(&newdt_local, &newdt, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get())) {
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
      
      timestep = dt_limit; // Change back to limiting timestep
    }while(running);
    
    iteration++; // Advance iteration number
    
    /// Write the restart file
    restart.write("%s/BOUT.restart.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());
    
    if((archive_restart > 0) && (iteration % archive_restart == 0)) {
      restart.write("%s/BOUT.restart_%04d.%d.%s", restartdir.c_str(), iteration, MYPE, restartext.c_str());
    }
    
    /// Call the monitor function
    
    if(call_monitors(simtime, s, nsteps)) {
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

void EulerSolver::take_step(BoutReal curtime, BoutReal dt, BoutReal *start, BoutReal *result) {
  
  load_vars(start);
  run_rhs(curtime);
  save_derivs(result);
  
  #pragma omp parallel for
  for(int i=0;i<nlocal;i++)
    result[i] = start[i] + dt*result[i];
}
