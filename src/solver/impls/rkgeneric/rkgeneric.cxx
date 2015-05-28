
#include "rkgeneric.hxx"
#include "rkschemefactory.hxx"
#include <bout/rkscheme.hxx>

#include <boutcomm.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include <output.hxx>

RKGenericSolver::RKGenericSolver(Options *options) : Solver(options) {
  f0 = 0; // Mark as uninitialised

  //Create scheme
  scheme=RKSchemeFactory::getInstance()->createRKScheme(options);
}

RKGenericSolver::~RKGenericSolver() {
  delete scheme;

  if(f0 != 0) {
    delete[] f0;
    delete[] f1;
    delete[] f2;
    delete[] tmpState;
  };
}

void RKGenericSolver::setMaxTimestep(BoutReal dt) {
  if(dt > timestep)
    return; // Already less than this
  
  if(adaptive)
    timestep = dt; // Won't be used this time, but next
}

int RKGenericSolver::init(bool restarting, int nout, BoutReal tstep) {

  int msg_point = msg_stack.push("Initialising RKGeneric solver");
  
  /// Call the generic initialisation first
  if(Solver::init(restarting, nout, tstep))
    return 1;

  output << "\n\tRunge-Kutta generic solver with scheme type "<<scheme->getType()<<"\n";

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
  f0 = new BoutReal[nlocal]; //Input
  f1 = new BoutReal[nlocal]; //Result--alternative order
  f2 = new BoutReal[nlocal]; //Result--follow order
  tmpState = new BoutReal[nlocal];

  // Put starting values into f0
  save_vars(f0);
  
  // Get options
  OPTION(options, atol, 1.e-5); // Absolute tolerance
  OPTION(options, rtol, 1.e-3); // Relative tolerance
  OPTION(options, max_timestep, tstep); // Maximum timestep
  OPTION(options, timestep, max_timestep); // Starting timestep
  OPTION(options, mxstep, 500); // Maximum number of steps between outputs
  OPTION(options, adaptive, true); // Prefer adaptive scheme

  msg_stack.pop(msg_point);

  //Initialise scheme
  scheme->init(nlocal);

  return 0;
}

int RKGenericSolver::run() {
  int msg_point = msg_stack.push("RKGenericSolver::run()");
  
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

	//Take a step
	take_step(simtime, dt, f0, f2, f1);

	//Calculate and check error if adaptive
        if(adaptive) {
	  BoutReal err;

	  //Get local part of relative error
	  BoutReal local_err = 0.;
	  for(int i=0;i<nlocal;i++) {
	    local_err += fabs(f2[i] - f1[i]) / ( fabs(f2[i]) + fabs(f1[i]) + atol );
	  }
	  //Reduce over procs
	  if(MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get())) {
	    throw BoutException("MPI_Allreduce failed");
	  }
	  //Normalise by number of values
          err /= (BoutReal) neq;

	  //Really the following should apply to both adaptive and non-adaptive
	  //approaches, but the non-adaptive can be determined without needing
	  //to do any solves so could perhaps be check during init instead.
          internal_steps++;
          if(internal_steps > mxstep)
            throw BoutException("ERROR: MXSTEP exceeded. timestep = %e, err=%e\n", timestep, err);

	  //Update the time step if required
          if((err > rtol) || (err < 0.1*rtol)) {

	    //Get new timestep
	    timestep=scheme->updateTimestep(timestep,rtol,err);

	    //Limit timestep to specified maximum
            if((max_timestep > 0) && (timestep > max_timestep))
              timestep = max_timestep;
          }

	  //If accuracy ok then break
          if(err < rtol) break;

        }else {
          // No adaptive timestepping so just accept step
          break;
        }
      }while(true);
      
      // Taken a step, swap buffers to put result into f0
      swap(f2, f0);
      simtime += dt;

      //Call the per internal timestep monitors
      call_timestep_monitors(simtime, dt);

    }while(running);
    
    load_vars(f0); // Put result into variables

    iteration++; // Advance iteration number
    
    /// Call the output step monitor function
    if(call_monitors(simtime, s, nsteps)) break; // Stop simulation
    
    // Reset iteration and wall-time count
    rhs_ncalls = 0;
  }
  
  msg_stack.pop(msg_point);
  
  return 0;
}

//Returns the evolved state vector along with an error estimate
void RKGenericSolver::take_step(const BoutReal timeIn, const BoutReal dt, const BoutReal *start, 
				BoutReal *resultFollow, BoutReal *resultAlt){

  //Calculate the intermediate stages
  for(int curStage=0;curStage<scheme->getStageCount();curStage++){
    //Use scheme to get this stage's time and state
    BoutReal curTime=scheme->setCurTime(timeIn,dt,curStage);
    scheme->setCurState(start, tmpState, nlocal, curStage, dt);

    //Get derivs for this stage
    load_vars(tmpState);
    run_rhs(curTime);
    save_derivs(scheme->steps[curStage]);
  };

  scheme->setOutputStates(start, resultFollow, resultAlt, nlocal, dt);

};
