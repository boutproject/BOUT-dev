
#ifdef BOUT_HAS_PETSC

#include "imex-bdf2.hxx"

#include <boutcomm.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include <output.hxx>

#include "petscsnes.h"

IMEXBDF2::IMEXBDF2(Options *opt) : Solver(opt), u(0) {
  
}

IMEXBDF2::~IMEXBDF2() {
  if(u) {
    delete[] u;
    delete[] u_1;
    delete[] u_2;
    
    delete[] f_1;
    delete[] f_2;
    
    delete[] rhs;
    
    VecDestroy(&snes_f);
    VecDestroy(&snes_x);
  }
}

/*
 * PETSc callback function, which evaluates the nonlinear
 * function to be solved by SNES.
 *
 * This function assumes the context void pointer is a pointer
 * to an IMEXBDF2 object.
 */ 
static PetscErrorCode FormFunction(SNES snes,Vec x, Vec f, void* ctx) {
  return static_cast<IMEXBDF2*>(ctx)->snes_function(x, f);
}

int IMEXBDF2::init(bool restarting, int nout, BoutReal tstep) {

  int msg_point = msg_stack.push("Initialising IMEX-BDF2 solver");
  
  /// Call the generic initialisation first
  if(Solver::init(restarting, nout, tstep))
    return 1;
  
  output << "\n\tIMEX-BDF2 time-integration solver\n";

  nsteps = nout; // Save number of output steps
  out_timestep = tstep;
  
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
  u = new BoutReal[nlocal];
  u_1 = new BoutReal[nlocal];
  u_2 = new BoutReal[nlocal];
  
  f_1 = new BoutReal[nlocal];
  f_2 = new BoutReal[nlocal];

  rhs = new BoutReal[nlocal];

  // Put starting values into u
  save_vars(u);
  
  // Get options
  OPTION(options, timestep, tstep); // Internal timestep
  OPTION(options, mxstep, 500); // Maximum number of steps between outputs
  
  ninternal = (int) (out_timestep / timestep);
  
  if((ninternal == 0) || (out_timestep / ninternal > timestep))
    ++ninternal;
  
  timestep = out_timestep / ninternal;
  output.write("\tUsing timestep = %e, %d internal steps per output\n", timestep, ninternal);

  // Initialise PETSc components
  int ierr;
  
  // Vectors
  ierr = VecCreate(BoutComm::get(), &snes_x);CHKERRQ(ierr);
  ierr = VecSetSizes(snes_x, nlocal, PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(snes_x);CHKERRQ(ierr);
  
  VecDuplicate(snes_x,&snes_f);
  
  // Nonlinear solver interface (SNES)
  SNESCreate(BoutComm::get(),&snes);
  
  // Set the callback function
  SNESSetFunction(snes,snes_f,FormFunction,this);
  
  // Set up the Jacobian
  //MatCreateSNESMF(snes,&Jmf);
  //SNESSetJacobian(snes,Jmf,Jmf,SNESComputeJacobianDefault,this);
  MatCreateAIJ(BoutComm::get(),
               nlocal,nlocal,  // Local sizes
               PETSC_DETERMINE, PETSC_DETERMINE, // Global sizes
               3,   // Number of nonzero entries in diagonal portion of local submatrix
               PETSC_NULL,
               0,   // Number of nonzeros per row in off-diagonal portion of local submatrix
               PETSC_NULL, 
               &Jmf);
  SNESSetJacobian(snes,Jmf,Jmf,SNESDefaultComputeJacobian,this);
  MatSetOption(Jmf,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);

  // Set tolerances
  BoutReal atol, rtol; // Tolerances for SNES solver
  options->get("atol", atol, 1e-16);
  options->get("rtol", rtol, 1e-10);
  SNESSetTolerances(snes,atol,rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

  // Predictor method
  options->get("predictor", predictor, 1);
  
  // Get runtime options
  SNESSetFromOptions(snes);
  
  msg_stack.pop(msg_point);

  return 0;
}

int IMEXBDF2::run() {
  int msg_point = msg_stack.push("IMEXBDF2::run()");
  
  // Multi-step scheme, so first steps are different
  bool starting = true;

  BoutReal dt = timestep;
  for(int s=0;s<nsteps;s++) {
    for(int i=0;i<ninternal;i++) {
      if(starting) {
        // Need to start multistep scheme using another method
        startup(simtime, dt);
        starting = false;
      }else {
        // Take a time step using IMEX-BDF2
        take_step(simtime, dt);
      }
      // No adaptive timestepping for now
      
      simtime += dt;
      
      call_timestep_monitors(simtime, dt);
    }
    
    load_vars(u); // Put result into variables
    run_rhs(simtime); // Run RHS to calculate auxilliary variables
    
    iteration++; // Advance iteration number
    
    /// Call the monitor function
    
    if(call_monitors(simtime, s, nsteps)) {
      // User signalled to quit
      break;
    }
    
    // Reset iteration and wall-time count
    rhs_ncalls = 0;
  }
  
  msg_stack.pop(msg_point);
  
  return 0;
}
/*!
 * Use forward-backward Euler to take a single step
 *
 * Inputs:
 * u   - Latest solution
 * 
 * Outputs:
 * u   - Latest solution
 * u_1 - Previous solution
 * f_1 - Time-derivative of previous solution
 * 
 */
void IMEXBDF2::startup(BoutReal curtime, BoutReal dt) {
  BoutReal *tmp;
  // Swap u and u_1
  tmp = u_1;
  u_1 = u;
  u = tmp;
  
  // Calculate time-derivative of u_1, put into f_1
  load_vars(u_1);
  run_convective(curtime);
  save_derivs(f_1);
  
  // Save to rhs vector
  for(int i=0;i<nlocal;i++)
    rhs[i] = u_1[i] + dt*f_1[i];
  
  switch(predictor) {
  case 1: {
    // Copy u_1 to u_2, since this will be used in predictor
    for(int i=0;i<nlocal;i++)
      u_2[i] = u_1[i];
    break;
  }
  case 2: {
    // Copy u_1 to u_2 and u, since these will be used in predictor
    for(int i=0;i<nlocal;i++)
      u[i] = u_2[i] = u_1[i];
    break;
  }
  }
  

  // Now need to solve u - dt*G(u) = rhs
  // Using run_diffusive as G
  solve_implicit(curtime+dt, dt);
}

/*!
 * Take a full IMEX-BDF2 step. Note that this assumes
 * that two time points are already available (in u and u_1). 
 * This therefore requires a startup step first
 * 
 * Inputs:
 * u   - Latest solution
 * u_1 - Previous solution
 * f_1 - Time-derivative of previous solution
 *
 * Outputs:
 * u   - Latest Solution
 * u_1 - Previous solution
 * f_1 - Time-derivative of previous solution
 * u_2 - Solution before last
 * f_2 - Time-derivative of u_2
 */
void IMEXBDF2::take_step(BoutReal curtime, BoutReal dt) {
  
  BoutReal *tmp;
  // Move f_1 to f_2, then f_1 will be overwritten with new time-derivative
  tmp = f_1;
  f_1 = f_2;
  f_2 = tmp;
  // Rotate u -> u_1, u_1 -> u_2, u_2 -> u . U later overwritten
  tmp = u_2;
  u_2 = u_1;
  u_1 = u;
  u = tmp;

  // Calculate time-derivative of u_1, put into f_1
  load_vars(u_1);
  run_convective(curtime);
  save_derivs(f_1);
  
  // Save to rhs vector
  for(int i=0;i<nlocal;i++)
    rhs[i] = (4./3)*u_1[i] - (1./3)*u_2[i] + (4./3)*dt*f_1[i] - (2./3)*dt*f_2[i];

  // Now need to solve u - (2./3)*dt*G(u) = rhs
  solve_implicit(curtime+dt, (2./3)*dt);
}

/*
 * Solves u - gamma*G(u) = rhs
 * 
 * where u is the result, G(u) is the stiff part of the rhs (run_diffusive)
 * and gamma is a factor depending on the time-step and method
 *
 * Inputs:
 * rhsvec  
 * 
 * 
 */
PetscErrorCode IMEXBDF2::solve_implicit(BoutReal curtime, BoutReal gamma) {
  implicit_curtime = curtime;
  implicit_gamma = gamma;
  
  // Set initial guess at the solution
  BoutReal *xdata;
  int ierr;
  ierr = VecGetArray(snes_x,&xdata);CHKERRQ(ierr);

  switch(predictor) {
  case 0: {
    // Constant, so next step is same as last step
    for(int i=0;i<nlocal;i++) {
      xdata[i] = u_1[i];     // Use previous solution
    }
    break;
  }
  case 1: {
    // Linear extrapolation from last two steps
    for(int i=0;i<nlocal;i++) {
      xdata[i] = 2.*u_1[i] - u_2[i];
    }
    break;
  }
  case 2: {
    // Quadratic extrapolation. Uses the fact that u has not yet been overwritten
    // and still contains u_3
    for(int i=0;i<nlocal;i++) {
      xdata[i] = 3.*u_1[i] - 3.*u_2[i] + u[i];
    }
  }
  default: {
    // Assume that there is no non-linear solve, so G = 0
    for(int i=0;i<nlocal;i++) {
      xdata[i] = rhs[i];   // If G = 0
    }
  }
  }
  ierr = VecRestoreArray(snes_x,&xdata);CHKERRQ(ierr);
  
  /*
  output << "Computing Jacobian\n";
  MatStructure  flag;
  implicit_curtime = curtime;
  implicit_gamma = gamma;
  SNESComputeFunction(snes, snes_x, snes_f);
  SNESComputeJacobian(snes,snes_x,&Jmf,&Jmf,&flag);
  MatView(Jmf, 	PETSC_VIEWER_STDOUT_SELF);
  */
  
  SNESSolve(snes,NULL,snes_x);
  
  // Find out if converged
  SNESConvergedReason reason;
  SNESGetConvergedReason(snes,&reason);
  if(reason < 0) {
    // Diverged
    throw BoutException("SNES failed to converge. Reason: %d\n", reason);
  }
  
  int its;
  SNESGetIterationNumber(snes,&its);
  
  //output << "Number of SNES iterations: " << its << endl;
  
  // Put the result into u
  ierr = VecGetArray(snes_x,&xdata);CHKERRQ(ierr);
  for(int i=0;i<nlocal;i++)
    u[i] = xdata[i];
  ierr = VecRestoreArray(snes_x,&xdata);CHKERRQ(ierr);
}

// f = (x - gamma*G(x)) - rhs
PetscErrorCode IMEXBDF2::snes_function(Vec x, Vec f) {
  
  BoutReal *xdata, *fdata;
  int ierr;
  
  // Get data from PETSc into BOUT++ fields
  ierr = VecGetArray(x,&xdata);CHKERRQ(ierr);
  
  load_vars(xdata);
  
  // Call RHS function
  run_diffusive(implicit_curtime);
  
  // Copy derivatives back
  ierr = VecGetArray(f,&fdata);CHKERRQ(ierr);
  save_derivs(fdata);
  
  // G(x) now in fdata
  for(int i=0;i<nlocal;i++) {
    //output.write("\n%d, %e, %e, %e ", i, xdata[i], fdata[i], rhs[i]);
    fdata[i] = xdata[i] - implicit_gamma * fdata[i] - rhs[i];
    //output.write("-> %e\n", fdata[i]);
  }
  
  // Restore data arrays to PETSc
  ierr = VecRestoreArray(f,&fdata);CHKERRQ(ierr);
  ierr = VecRestoreArray(x,&xdata);CHKERRQ(ierr);
  
  return 0;
}

#endif // BOUT_HAS_PETSC
