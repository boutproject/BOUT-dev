
#include "rk4.hxx"

#include <boutcomm.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <bout/openmpwrap.hxx>

#include <cmath>

#include <output.hxx>

RK4Solver::RK4Solver(Options *options) : Solver(options) { canReset = true; }

void RK4Solver::setMaxTimestep(BoutReal dt) {
  if (dt > timestep)
    return; // Already less than this
  
  if (adaptive)
    timestep = dt; // Won't be used this time, but next
}

int RK4Solver::init(int nout, BoutReal tstep) {

  TRACE("Initialising RK4 solver");
  
  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;
  
  output << "\n\tRunge-Kutta 4th-order solver\n";

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
  
  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n",
	       n3Dvars(), n2Dvars(), neq, nlocal);
  
  // Allocate memory
  f0.reallocate(nlocal);
  f1.reallocate(nlocal);
  f2.reallocate(nlocal);

  // memory for taking a single time step
  k1.reallocate(nlocal);
  k2.reallocate(nlocal);
  k3.reallocate(nlocal);
  k4.reallocate(nlocal);
  k5.reallocate(nlocal);

  // Put starting values into f0
  save_vars(std::begin(f0));

  // Get options
  atol = (*options)["atol"].doc("Absolute tolerance").withDefault(1.e-5);
  rtol = (*options)["rtol"].doc("Relative tolerance").withDefault(1.e-3);
  max_timestep = (*options)["max_timestep"].doc("Maximum timestep").withDefault(tstep);
  timestep = (*options)["timestep"].doc("Starting timestep").withDefault(max_timestep);
  mxstep = (*options)["mxstep"].doc("Maximum number of steps between outputs").withDefault(500);
  adaptive = (*options)["adaptive"].doc("Adapt internal timestep using ATOL and RTOL.").withDefault(false);

  return 0;
}

int RK4Solver::run() {
  TRACE("RK4Solver::run()");
  
  for(int s=0;s<nsteps;s++) {
    BoutReal target = simtime + out_timestep;
    
    BoutReal dt;
    bool running = true;
    int internal_steps = 0;
    do {
      // Take a single time step
      
      do {
        dt = timestep;
        running = true;
        if((simtime + dt) >= target) {
          dt = target - simtime; // Make sure the last timestep is on the output 
          running = false;
        }
        if(adaptive) {
          // Take two half-steps
          take_step(simtime,          0.5*dt, f0, f1);
          take_step(simtime + 0.5*dt, 0.5*dt, f1, f2);
          
          // Take a full step
          take_step(simtime, dt, f0, f1);
          
          // Check accuracy
          BoutReal local_err = 0.;
          BOUT_OMP(parallel for reduction(+: local_err)   )
          for(int i=0;i<nlocal;i++) {
            local_err += fabs(f2[i] - f1[i]) / ( fabs(f1[i]) + fabs(f2[i]) + atol );
          }
        
          // Average over all processors
          BoutReal err;
          if(MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get())) {
            throw BoutException("MPI_Allreduce failed");
          }

          err /= static_cast<BoutReal>(neq);

          internal_steps++;
          if(internal_steps > mxstep)
            throw BoutException("ERROR: MXSTEP exceeded. timestep = {:e}, err={:e}\n",
                                timestep, err);

          if((err > rtol) || (err < 0.1*rtol)) {
            // Need to change timestep. Error ~ dt^5
            timestep /= pow(err / (0.5*rtol), 0.2);
            
            if((max_timestep > 0) && (timestep > max_timestep))
              timestep = max_timestep;
          }
          if(err < rtol) {
            break; // Acceptable accuracy
          }
        }else {
          // No adaptive timestepping
          take_step(simtime, dt, f0, f2);
          break;
        }
      }while(true);
      
      // Taken a step, swap buffers
      swap(f2, f0);
      simtime += dt;
      
      call_timestep_monitors(simtime, dt);
    }while(running);

    load_vars(std::begin(f0)); // Put result into variables
    // Call rhs function to get extra variables at this time
    run_rhs(simtime);
    
    iteration++; // Advance iteration number
    
    /// Call the monitor function
    
    if(call_monitors(simtime, s, nsteps)) {
      break; // Stop simulation
    }
  }
  
  return 0;
}

void RK4Solver::resetInternalFields(){
  //Zero out history
  for(int i=0;i<nlocal;i++){
    f1[i]=0; f2[i]=0;
  }
  
  //Copy fields into current step
  save_vars(std::begin(f0));
}

void RK4Solver::take_step(BoutReal curtime, BoutReal dt, Array<BoutReal> &start,
                          Array<BoutReal> &result) {

  load_vars(std::begin(start));
  run_rhs(curtime);
  save_derivs(std::begin(k1));

  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    k5[i] = start[i] + 0.5*dt*k1[i];

  load_vars(std::begin(k5));
  run_rhs(curtime + 0.5*dt);
  save_derivs(std::begin(k2));

  BOUT_OMP(parallel for )
  for(int i=0;i<nlocal;i++)
    k5[i] = start[i] + 0.5*dt*k2[i];

  load_vars(std::begin(k5));
  run_rhs(curtime + 0.5*dt);
  save_derivs(std::begin(k3));

  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    k5[i] = start[i] + dt*k3[i];

  load_vars(std::begin(k5));
  run_rhs(curtime + dt);
  save_derivs(std::begin(k4));

  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    result[i] = start[i] + (1./6.)*dt*(k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i]);
}
