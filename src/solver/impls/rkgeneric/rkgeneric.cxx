
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
  //Create scheme
  scheme=RKSchemeFactory::getInstance()->createRKScheme(options);
  canReset = true;
}

RKGenericSolver::~RKGenericSolver() {
  delete scheme;
}

void RKGenericSolver::setMaxTimestep(BoutReal dt) {
  if(dt > timestep)
    return; // Already less than this
  
  if(adaptive)
    timestep = dt; // Won't be used this time, but next
}

int RKGenericSolver::init(int nout, BoutReal tstep) {

  TRACE("Initialising RKGeneric solver");
  
  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;

  output << "\n\tRunge-Kutta generic solver with scheme type "<<scheme->getType()<<"\n";

  nsteps = nout; // Save number of output steps
  out_timestep = tstep;
  max_dt = tstep;
  
  // Calculate number of variables
  nlocal = getLocalN();
  
  // Get total problem size
  int ntmp;
  if(bout::globals::mpi->MPI_Allreduce(&nlocal, &ntmp, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed!");
  }
  neq = ntmp;
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3Dvars(), n2Dvars(), neq, nlocal);
  
  // Get options
  atol = (*options)["atol"].doc("Absolute tolerance").withDefault(1.e-5);
  rtol = (*options)["rtol"].doc("Relative tolerance").withDefault(1.e-3);
  max_timestep = (*options)["max_timestep"].doc("Maximum timestep").withDefault(tstep);
  timestep = (*options)["timestep"].doc("Starting timestep").withDefault(max_timestep);
  mxstep = (*options)["mxstep"].doc("Maximum number of steps between outputs").withDefault(500);
  adaptive = (*options)["adaptive"].doc("Adapt internal timestep using ATOL and RTOL.").withDefault(true);

  // Allocate memory
  f0.reallocate(nlocal); // Input
  f2.reallocate(nlocal); // Result--follow order
  tmpState.reallocate(nlocal);

  // Put starting values into f0
  save_vars(std::begin(f0));

  //Initialise scheme
  scheme->init(nlocal,neq,adaptive,atol,rtol,options);

  return 0;
}

void RKGenericSolver::resetInternalFields(){
  //Zero out history
  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++){
    tmpState[i]=0; f2[i]=0;
  }
  
  //Copy fields into current step
  save_vars(std::begin(f0));
}

int RKGenericSolver::run() {
  TRACE("RKGenericSolver::run()");
  
  for(int s=0;s<nsteps;s++) {
    BoutReal target = simtime + out_timestep;
    
    BoutReal dt;
    bool running = true;
    int internal_steps = 0;

    // Take a single output time step
    do {
      // Take a single internal time step
      do {
        dt = timestep;
        running = true;
        if((simtime + dt) >= target) {
          dt = target - simtime; // Make sure the last timestep is on the output 
          running = false;
        }

	BoutReal err;

	//Take a step
	err = take_step(simtime, dt, f0, f2);

	//Calculate and check error if adaptive
        if(adaptive) {
	  //Really the following should apply to both adaptive and non-adaptive
	  //approaches, but the non-adaptive can be determined without needing
	  //to do any solves so could perhaps be check during init instead.
          internal_steps++;
          if(internal_steps > mxstep)
            throw BoutException("ERROR: MXSTEP exceeded. timestep = %e, err=%e\n", timestep, err);

	  //Update the time step if required, note we ignore increases to the timestep
	  //when on the last internal step as here we may have an artificially small dt
	  if((err > rtol) || ((err < 0.1*rtol) && running)) {
	    
	    //Get new timestep
	    timestep=scheme->updateTimestep(dt,err);

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

    load_vars(std::begin(f0)); // Put result into variables

    run_rhs(simtime); //Ensure aux. variables are up to date

    iteration++; // Advance iteration number
    
    /// Call the output step monitor function
    if(call_monitors(simtime, s, nsteps)) break; // Stop simulation
  }
  
  return 0;
}

//Returns the evolved state vector along with an error estimate
BoutReal RKGenericSolver::take_step(const BoutReal timeIn, const BoutReal dt,
                                    const Array<BoutReal> &start,
                                    Array<BoutReal> &resultFollow) {

  //Calculate the intermediate stages
  for(int curStage=0;curStage<scheme->getStageCount();curStage++){
    //Use scheme to get this stage's time and state
    BoutReal curTime=scheme->setCurTime(timeIn,dt,curStage);
    scheme->setCurState(start, tmpState, curStage, dt);

    //Get derivs for this stage
    load_vars(std::begin(tmpState));
    run_rhs(curTime);
    save_derivs(&(scheme->steps(curStage,0)));
  }

  return scheme->setOutputStates(start, dt, resultFollow);
}
