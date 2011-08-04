
#include "euler.hxx"

#include <utils.hxx>
#include <boutexception.hxx>

#include <cmath>

EulerSolver::EulerSolver() : Solver() {
  
}

EulerSolver::~EulerSolver() {

}

void EulerSolver::setMaxTimestep(BoutReal dt) {
  if(dt > timestep)
    return; // Already less than this
  
  timestep = dt;
  timestep_reduced = true;
}

int EulerSolver::init(rhsfunc f, int argc, char **argv, bool restarting, int nout, BoutReal tstep) {
#ifdef CHECK
  int msg_point = msg_stack.push("Initialising Euler solver");
#endif
  
  /// Call the generic initialisation first
  if(Solver::init(f, argc, argv, restarting, nout, tstep))
    return 1;
  
  output << "\n\tEuler solver\n";

  nsteps = nout; // Save number of output steps
  out_timestep = tstep;
  
  // Get options
  Options *options = Options::getRoot();
  options = options->getSection("solver");
  OPTION(options, start_timestep, tstep);
  OPTION(options, mxstep, 500); // Maximum number of steps between outputs
  
  // Calculate number of variables
  nlocal = getLocalN();
  
  // Get total problem size
  int neq;
  if(MPI_Allreduce(&nlocal, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    output.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3Dvars(), n2Dvars(), neq, nlocal);
  
  // Allocate memory
  f0 = new BoutReal[nlocal];
  f1 = new BoutReal[nlocal];
  
  // Put starting values into f0
  save_vars(f0);
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}

int EulerSolver::run(MonitorFunc monitor) {
#ifdef CHECK
  int msg_point = msg_stack.push("EulerSolver::run()");
#endif
  
  timestep = start_timestep;
  
  for(int s=0;s<nsteps;s++) {
    BoutReal target = simtime + out_timestep;
    
    bool running = true;
    int internal_steps = 0;
    do {
      // Take a step
      if((simtime + timestep) >= target) {
        // Make sure the last timestep is on the output 
        timestep = target - simtime; 
        running = false;
      }
      
      timestep_reduced = false;
      take_step(simtime, timestep, f0, f1);
      
      // If timestep_reduced re-run
      if(timestep_reduced)
        take_step(simtime, timestep, f0, f1);
      
      // Taken a step, swap buffers
      SWAP(f1, f0);
      simtime += timestep;

      internal_steps++;
      if(internal_steps > mxstep)
        throw BoutException("ERROR: MXSTEP exceeded. simtime=%e, timestep = %e\n", simtime, timestep);
      
    }while(running);
    
    iteration++; // Advance iteration number
    
    /// Write the restart file
    restart.write("%s/BOUT.restart.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());
    
    if((archive_restart > 0) && (iteration % archive_restart == 0)) {
      restart.write("%s/BOUT.restart_%04d.%d.%s", restartdir.c_str(), iteration, MYPE, restartext.c_str());
    }
    
    /// Call the monitor function
    
    if(monitor(simtime, s, nsteps)) {
      // User signalled to quit
      
      // Write restart to a different file
      restart.write("%s/BOUT.final.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());
      
      output.write("Monitor signalled to quit. Returning\n");
      break;
    }
    
    // Reset iteration and wall-time count
    rhs_ncalls = 0;
    rhs_wtime = 0.0;
  }
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
  
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
