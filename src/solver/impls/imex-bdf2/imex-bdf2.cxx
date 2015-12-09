
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

/*!
 * PETSc callback function, which evaluates the nonlinear
 * function to be solved by SNES.
 *
 * This function assumes the context void pointer is a pointer
 * to an IMEXBDF2 object.
 */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
static PetscErrorCode FormFunction(SNES snes,Vec x, Vec f, void* ctx) {
  return static_cast<IMEXBDF2*>(ctx)->snes_function(x, f, false);
}

/*!
 * PETSc callback function for forming Jacobian
 *
 * This function can be a linearised form of FormFunction
 */
#undef __FUNCT__
#define __FUNCT__ "FormFunctionForDifferencing"
static PetscErrorCode FormFunctionForDifferencing(void* ctx, Vec x, Vec f) {
  return static_cast<IMEXBDF2*>(ctx)->snes_function(x, f, true);
}

#undef __FUNCT__
#define __FUNCT__ "imexbdf2PCapply"
static PetscErrorCode imexbdf2PCapply(PC pc,Vec x,Vec y) {
  int ierr;

  // Get the context
  IMEXBDF2 *s;
  ierr = PCShellGetContext(pc,(void**)&s);CHKERRQ(ierr);

  PetscFunctionReturn(s->precon(x, y));
}


/*!
 * Initialisation routine. Called once before solve.
 *
 */
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
  saveVars(u);

  // Get options
  OPTION(options, timestep, tstep); // Internal timestep

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

  /////////////////////////////////////////////////////
  // Set up the Jacobian

  bool matrix_free;
  OPTION(options, matrix_free, true); // Default is matrix free
  if(matrix_free) {
    /*!
      PETSc SNES matrix free Jacobian, using a different
      operator for differencing.

      See PETSc examples
      http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tests/ex7.c.html
      and this thread:
      http://lists.mcs.anl.gov/pipermail/petsc-users/2014-January/020075.html

     */
    MatCreateSNESMF(snes,&Jmf);

    // Set a function to be called for differencing
    // This can be a linearised form of the SNES function
    MatMFFDSetFunction(Jmf,FormFunctionForDifferencing,this);

    // Calculate Jacobian matrix free using FormFunctionForDifferencing
    SNESSetJacobian(snes,Jmf,Jmf,MatMFFDComputeJacobian,this);
  }else {
    /*!
     * Calculate Jacobian using finite differences.
     * NOTE: Slow!
     */
    MatCreateAIJ(BoutComm::get(),
                 nlocal,nlocal,  // Local sizes
                 PETSC_DETERMINE, PETSC_DETERMINE, // Global sizes
                 3,   // Number of nonzero entries in diagonal portion of local submatrix
                 PETSC_NULL,
                 0,   // Number of nonzeros per row in off-diagonal portion of local submatrix
                 PETSC_NULL,
                 &Jmf);

#ifdef BOUT_HAS_PETSC_3_3
    // Before 3.4
    SNESSetJacobian(snes,Jmf,Jmf,SNESDefaultComputeJacobian,this);
#else
    SNESSetJacobian(snes,Jmf,Jmf,SNESComputeJacobianDefault,this);
#endif

    MatSetOption(Jmf,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
  }

  /////////////////////////////////////////////////////
  // Set tolerances
  BoutReal atol, rtol; // Tolerances for SNES solver
  options->get("atol", atol, 1e-16);
  options->get("rtol", rtol, 1e-10);
  SNESSetTolerances(snes,atol,rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

  /////////////////////////////////////////////////////
  // Predictor method
  options->get("predictor", predictor, 1);

  /////////////////////////////////////////////////////
  // Preconditioner

  bool use_precon;
  OPTION(options, use_precon,   false);
  if(use_precon && have_user_precon()) {
    output.write("\tUsing user-supplied preconditioner\n");

    // Get KSP context from SNES
    KSP ksp;
    SNESGetKSP(snes, &ksp);

    // Get PC context from KSP
    PC pc;
    KSPGetPC(ksp,&pc);

    // Set a Shell (matrix-free) preconditioner type
    PCSetType(pc, PCSHELL);

    // Specify the preconditioner function
    PCShellSetApply(pc,imexbdf2PCapply);
    // Context used to supply object pointer
    PCShellSetContext(pc,this);
  }

  /////////////////////////////////////////////////////
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

    loadVars(u);// Put result into variables
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
  loadVars(u_1);
  run_convective(curtime);
  saveDerivs(f_1);

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

  loadVars(u_1);
  run_convective(curtime);
  saveDerivs(f_1);

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
    break;
  }
  default: {
    // Assume that there is no non-linear solve, so G = 0
    for(int i=0;i<nlocal;i++) {
      xdata[i] = rhs[i];   // If G = 0
    }
  }
  }
  //output.write("\nIMEX: Solving, %e, %e, %e, (%e)\n", u[0], u_2[0], u_1[0], xdata[0]);

  ierr = VecRestoreArray(snes_x,&xdata);CHKERRQ(ierr);

  /*
  output << "Computing Jacobian\n";
  MatStructure  flag;
  implicit_curtime = curtime;
  implicit_gamma = gamma;
  SNESComputeFunction(snes, snes_x, snes_f);
  SNESComputeJacobian(snes,snes_x,&Jmf,&Jmf,&flag);
  MatView(Jmf,  PETSC_VIEWER_STDOUT_SELF);
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
  //output.write("\nIMEX: Done -> %e\n", xdata[0]);

  for(int i=0;i<nlocal;i++)
    u[i] = xdata[i];
  ierr = VecRestoreArray(snes_x,&xdata);CHKERRQ(ierr);
}

// f = (x - gamma*G(x)) - rhs
PetscErrorCode IMEXBDF2::snes_function(Vec x, Vec f, bool linear) {
  BoutReal *xdata, *fdata;
  int ierr;

  // Get data from PETSc into BOUT++ fields
  ierr = VecGetArray(x,&xdata);CHKERRQ(ierr);

  loadVars(xdata);

  // Call RHS function
  run_diffusive(implicit_curtime, linear);

  // Copy derivatives back
  ierr = VecGetArray(f,&fdata);CHKERRQ(ierr);
  saveDerivs(fdata);

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

/*
 * Preconditioner function
 */
PetscErrorCode IMEXBDF2::precon(Vec x, Vec f) {
  if(!have_user_precon()) {
    // No user preconditioner
    throw BoutException("No user preconditioner");
  }

  int ierr;

  // Get data from PETSc into BOUT++ fields
  Vec solution;
  SNESGetSolution(snes, &solution);
  BoutReal *soldata;
  ierr = VecGetArray(x,&soldata);CHKERRQ(ierr);
  load_vars(soldata);
  ierr = VecRestoreArray(solution,&soldata);CHKERRQ(ierr);

  // Load vector to be inverted into ddt() variables
  BoutReal *xdata;
  ierr = VecGetArray(x,&xdata);CHKERRQ(ierr);
  load_derivs(xdata);
  ierr = VecRestoreArray(x,&xdata);CHKERRQ(ierr);

  // Run the preconditioner
  run_precon(implicit_curtime, implicit_gamma, 0.0);

  // Save the solution from F_vars
  BoutReal *fdata;
  ierr = VecGetArray(f,&fdata);CHKERRQ(ierr);
  save_derivs(fdata);
  ierr = VecRestoreArray(f,&fdata);CHKERRQ(ierr);

  return 0;
}

/*!
 * Loop over arrays, using template parameter
 * to specify the operation to be performed at each point
 *
 */
template< class Op >
void IMEXBDF2::loopVars(BoutReal *u) {
  // Loop over 2D variables
  for(vector< VarStr<Field2D> >::const_iterator it = f2d.begin(); it != f2d.end(); ++it) {
    Op op(it->var, it->F_var); // Initialise the operator

    if(it->evolve_bndry) {
      // Include boundary regions

      // Inner X
      if(mesh->firstX() && !mesh->periodicX) {
        for(int jx=0;jx<mesh->xstart;++jx)
          for(int jy=mesh->ystart;jy<=mesh->yend;++jy) {
            op.run(jx, jy, u); ++u;
          }
      }

      // Outer X
      if(mesh->lastX() && !mesh->periodicX) {
        for(int jx=mesh->xend+1;jx<mesh->ngx;++jx)
          for(int jy=mesh->ystart;jy<=mesh->yend;++jy) {
            op.run(jx, jy, u); ++u;
          }
      }
      // Lower Y
      for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); ++xi) {
        for(int jy=0;jy<mesh->ystart;++jy) {
          op.run(*xi, jy, u); ++u;
        }
      }

      // Upper Y
      for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); ++xi) {
        for(int jy=mesh->yend+1;jy<mesh->ngy;++jy) {
          op.run(*xi, jy, u); ++u;
        }
      }
    }

    // Bulk of points
    for(int jx=mesh->xstart; jx <= mesh->xend; ++jx)
      for(int jy=mesh->ystart; jy <= mesh->yend; ++jy) {
        op.run(jx, jy, u); ++u;
      }
  }

  // Loop over 3D variables
  for(vector< VarStr<Field3D> >::const_iterator it = f3d.begin(); it != f3d.end(); ++it) {
    Op op(it->var, it->F_var); // Initialise the operator
    if(it->evolve_bndry) {
      // Include boundary regions

      // Inner X
      if(mesh->firstX() && !mesh->periodicX) {
        for(int jx=0;jx<mesh->xstart;++jx)
          for(int jy=mesh->ystart;jy<=mesh->yend;++jy)
            for(int jz=0; jz < mesh->ngz-1; ++jz) {
              op.run(jx, jy, jz, u); ++u;
            }
      }

      // Outer X
      if(mesh->lastX() && !mesh->periodicX) {
        for(int jx=mesh->xend+1;jx<mesh->ngx;++jx)
          for(int jy=mesh->ystart;jy<=mesh->yend;++jy)
            for(int jz=0; jz < mesh->ngz-1; ++jz) {
              op.run(jx, jy, jz, u); ++u;
            }
      }
      // Lower Y
      for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); ++xi) {
        for(int jy=0;jy<mesh->ystart;++jy)
          for(int jz=0; jz < mesh->ngz-1; ++jz) {
            op.run(*xi, jy, jz, u); ++u;
          }
      }

      // Upper Y
      for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); ++xi) {
        for(int jy=mesh->yend+1;jy<mesh->ngy;++jy)
          for(int jz=0; jz < mesh->ngz-1; ++jz) {
            op.run(*xi, jy, jz, u); ++u;
          }
      }
    }

    // Bulk of points
    for(int jx=mesh->xstart; jx <= mesh->xend; ++jx)
      for(int jy=mesh->ystart; jy <= mesh->yend; ++jy)
        for(int jz=0; jz < mesh->ngz-1; ++jz) {
          op.run(jx, jy, jz, u); ++u;
        }
  }
}

///////////////////////////////////////////////////////////////////

class SaveVarOp {
public:
  // Initialise with a Field2D iterator
  SaveVarOp(Field2D *var, Field2D *F_var) : var2D(var) {}
  // Initialise with a Field3D iterator
  SaveVarOp(Field3D *var, Field3D *F_var) : var3D(var) {}

  // Perform operation on 2D field
  inline void run(int jx, int jy, BoutReal *u) {
    *u = (*var2D)(jx,jy);
  }

  // Perform operation on 3D field
  inline void run(int jx, int jy, int jz, BoutReal *u) {
    *u = (*var3D)(jx,jy,jz);
  }
private:
  Field2D *var2D;
  Field3D *var3D;
};

/*!
 * Copy data from fields into array
 */
void IMEXBDF2::saveVars(BoutReal *u) {
  //loopVars<SaveVarOp>(u);
  save_vars(u);
}

///////////////////////////////////////////////////////////////////

class LoadVarOp {
public:
  // Initialise with a Field2D iterator
  LoadVarOp(Field2D *var, Field2D *F_var) : var2D(var) {}
  // Initialise with a Field3D iterator
  LoadVarOp(Field3D *var, Field3D *F_var) : var3D(var) {}

  // Perform operation on 2D field
  inline void run(int jx, int jy, BoutReal *u) {
    (*var2D)(jx,jy) = *u;
  }

  // Perform operation on 3D field
  inline void run(int jx, int jy, int jz, BoutReal *u) {
    (*var3D)(jx,jy,jz) = *u;
  }
private:
  Field2D *var2D;
  Field3D *var3D;
};

/*!
 * Copy data from array into fields
 */
void IMEXBDF2::loadVars(BoutReal *u) {
  //loopVars<LoadVarOp>(u);
  load_vars(u);
}

///////////////////////////////////////////////////////////////////

class SaveDerivsOp {
public:
  // Initialise with a Field2D iterator
  SaveDerivsOp(Field2D *var, Field2D *F_var) : F_var2D(F_var) {}
  // Initialise with a Field3D iterator
  SaveDerivsOp(Field3D *var, Field3D *F_var) : F_var3D(F_var) {}

  // Perform operation on 2D field
  inline void run(int jx, int jy, BoutReal *u) {
    *u = (*F_var2D)(jx,jy);
  }

  // Perform operation on 3D field
  inline void run(int jx, int jy, int jz, BoutReal *u) {
    *u = (*F_var3D)(jx,jy,jz);
  }
private:
  Field2D *F_var2D;
  Field3D *F_var3D;
};

/*!
 * Copy time derivatives from fields into array
 */
void IMEXBDF2::saveDerivs(BoutReal *u) {
  //loopVars<SaveDerivsOp>(u);
  save_derivs(u);
}

#endif // BOUT_HAS_PETSC
