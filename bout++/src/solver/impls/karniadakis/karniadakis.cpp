
#include "karniadakis.h"

KarniadakisSolver::KarniadakisSolver() : Solver()
{
  
}

KarniadakisSolver::~KarniadakisSolver()
{
  
}

int KarniadakisSolver::init(rhsfunc f, int argc, char **argv, bool restarting, int nout, BoutReal tstep)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Initialising Karniadakis solver");
#endif
  
  /// Call the generic initialisation first
  if(Solver::init(f, argc, argv, restarting, nout, tstep))
    return 1;
  
  nsteps = nout; // Save number of output steps
  out_timestep = tstep;
  
  // Choose timestep
  if(max_dt < 0.0) {
    output << "\tWARNING: Starting dt not set\n";
    max_dt = tstep;
  }
  
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}


void KarniadakisSolver::take_step(BoutReal dt)
{
  // S0 = S(f0)
  
  load_vars(f0);
  run_convective(time);
  save_derivs(S0);
  
  // f1 = (6./11.) * (3.*f0 - 1.5*fm1 + (1./3.)*fm2 + dt*(3.*S0 - 3.*Sm1 + Sm2))
  
  for(int i=0;i<nlocal;i++)
    f1[i] = (6./11.) * (3.*f0[i] - 1.5*fm1[i] + (1./3.)*fm2[i] + dt*(3.*S0[i] - 3.*Sm1[i] + Sm2[i]));
  
  // D0 = S(f0)
  load_vars(f0);
  run_diffusive(time);
  save_derivs(D0);
  
  // f1 = f1 + dt*D0
  for(int i=0;i<nlocal;i++)
    f1[i] += dt*D0[i];
}

int KarniadakisSolver::run(MonitorFunc monitor)
{
#ifdef CHECK
  int msg_point = msg_stack.push("KarniadakisSolver::run()");
#endif
  
  for(int i=0;i<nsteps;i++) {
    BoutReal target = time + out_timestep;
    
    // Run until reach target time
    BoutReal dt;
    bool running = true;
    do {
      dt = timestep;
      if((time + dt) >= target) {
	dt = target - time; // Make sure the last timestep is on the output 
	running = false;
      }
      take_step(dt);
      
      // Cycle buffers
      BoutReal *tmp = fm2;
      fm2 = fm1;
      fm1 = f0;
      f0 = f1;
      f1 = tmp;
      
      tmp = Sm2;
      Sm2 = Sm1;
      Sm1 = S0;
      S0 = tmp;
      
      time += dt;
    }while(running);
    
    /// Write the restart file
    restart.write("%s/BOUT.restart.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());
    
    if((archive_restart > 0) && (iteration % archive_restart == 0)) {
      restart.write("%s/BOUT.restart_%04d.%d.%s", restartdir.c_str(), iteration, MYPE, restartext.c_str());
    }
    
    /// Call the monitor function
    
    if(monitor(time, i, nsteps)) {
      // User signalled to quit
      
      // Write restart to a different file
      restart.write("%s/BOUT.final.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());
      
      output.write("Monitor signalled to quit. Returning\n");
      break;
    }
  }
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
  
  return 0;
}

