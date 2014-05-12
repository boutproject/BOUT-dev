
#include "rk3-ssp.hxx"

#include <boutcomm.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include <output.hxx>

RK3SSP::RK3SSP(Options *opt) : Solver(opt), f(0) {
  
}

RK3SSP::~RK3SSP() {
  if(f != 0) {
    delete[] f;
    
    delete[] u1;
    delete[] u2;
    delete[] u3;
    delete[] L;
  }
}

void RK3SSP::setMaxTimestep(BoutReal dt) {
  if(dt > timestep)
    return; // Already less than this
  
  timestep = dt; // Won't be used this time, but next
}

int RK3SSP::init(bool restarting, int nout, BoutReal tstep) {

  int msg_point = msg_stack.push("Initialising RK3 SSP solver");
  
  /// Call the generic initialisation first
  if(Solver::init(restarting, nout, tstep))
    return 1;
  
  output << "\n\tRunge-Kutta 3rd-order SSP solver\n";

  nsteps = nout; // Save number of output steps
  out_timestep = tstep;
  max_dt = tstep;
  
  // Calculate number of variables
  nlocal = getLocalN();
  
  // Get total problem size
  int ntmp;
  if(MPI_Allreduce(&nlocal, &ntmp, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed!");
  }
  neq = ntmp;
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3Dvars(), n2Dvars(), neq, nlocal);
  
  // Allocate memory
  f = new BoutReal[nlocal];
  
  // memory for taking a single time step
  u1 = new BoutReal[nlocal];
  u2 = new BoutReal[nlocal];
  u3 = new BoutReal[nlocal];
  L = new BoutReal[nlocal];

  // Put starting values into f
  save_vars(f);
  
  // Get options
  OPTION(options, max_timestep, tstep); // Maximum timestep
  OPTION(options, timestep, max_timestep); // Starting timestep
  OPTION(options, mxstep, 500); // Maximum number of steps between outputs

  msg_stack.pop(msg_point);

  return 0;
}

int RK3SSP::run() {
  int msg_point = msg_stack.push("RK3SSP::run()");
  
  for(int s=0;s<nsteps;s++) {
    BoutReal target = simtime + out_timestep;
    
    BoutReal dt;
    bool running = true;
    int internal_steps = 0;
    do {
      // Take a single time step
      
      dt = timestep;
      running = true;
      if((simtime + dt) >= target) {
        dt = target - simtime; // Make sure the last timestep is on the output 
        running = false;
      }
        
      // No adaptive timestepping for now
      take_step(simtime, dt, f, f);
      
      simtime += dt;
      
      call_timestep_monitors(simtime, dt);
    }while(running);
    
    iteration++; // Advance iteration number
    
    /// Write the restart file
    restart.write();
    
    if((archive_restart > 0) && (iteration % archive_restart == 0)) {
      restart.write("%s/BOUT.restart_%04d.%s", restartdir.c_str(), iteration, restartext.c_str());
    }
    
    /// Call the monitor function
    
    if(call_monitors(simtime, s, nsteps)) {
      // User signalled to quit
      
      // Write restart to a different file
      restart.write("%s/BOUT.final.%s", restartdir.c_str(), restartext.c_str());
      
      output.write("Monitor signalled to quit. Returning\n");
      break;
    }
    
    // Reset iteration and wall-time count
    rhs_ncalls = 0;
  }
  
  msg_stack.pop(msg_point);
  
  return 0;
}

void RK3SSP::take_step(BoutReal curtime, BoutReal dt, BoutReal *start, BoutReal *result) {
  
  load_vars(start);
  run_rhs(curtime);
  save_derivs(L);
  
  #pragma omp parallel for
  for(int i=0;i<nlocal;i++)
    u1[i] = start[i] + dt*L[i];
  
  load_vars(u1);
  run_rhs(curtime + dt);
  save_derivs(L);
  
  #pragma omp parallel for 
  for(int i=0;i<nlocal;i++)
    u2[i] = 0.75*start[i] + 0.25*u1[i] + 0.25*dt*L[i];
  
  load_vars(u2);
  run_rhs(curtime + 0.5*dt);
  save_derivs(L);
 
  #pragma omp parallel for
  for(int i=0;i<nlocal;i++)
    result[i] = (1./3)*start[i] + (2./3.)*(u2[i] + dt*L[i]);
}
