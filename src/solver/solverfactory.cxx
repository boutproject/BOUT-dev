#include "bout/solverfactory.hxx"


#include "impls/cvode/cvode.hxx"
#include "impls/arkode/arkode.hxx"
#include "impls/petsc-3.1/petsc-3.1.hxx"
#include "impls/petsc-3.2/petsc-3.2.hxx"
#include "impls/petsc-3.3/petsc-3.3.hxx"
#include "impls/petsc-3.4/petsc-3.4.hxx"
#include "impls/petsc-3.5/petsc-3.5.hxx"
#include "impls/petsc/petsc.hxx"
#include "impls/slepc-3.4/slepc-3.4.hxx"
#include "impls/ida/ida.hxx"
#include "impls/pvode/pvode.hxx"
#include "impls/karniadakis/karniadakis.hxx"
#include "impls/rk4/rk4.hxx"
#include "impls/euler/euler.hxx"
#include "impls/rk3-ssp/rk3-ssp.hxx"
#include "impls/power/power.hxx"
#include "impls/imex-bdf2/imex-bdf2.hxx"
#include "impls/snes/snes.hxx"
#include "impls/rkgeneric/rkgeneric.hxx"

#include <boutexception.hxx>

SolverFactory* SolverFactory::instance = NULL;

SolverFactory* SolverFactory::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new SolverFactory();
  }
  return instance;
}

inline SolverType SolverFactory::getDefaultSolverType() {
  SolverType type = NULL;

  #if defined BOUT_HAS_CVODE
    type = SOLVERCVODE;
  #elif defined BOUT_HAS_IDA
    type = SOLVERIDA;
    //#elif defined BOUT_HAS_PETSC
    //type = SOLVERPETSC;
  #else
    type = SOLVERPVODE;
  #endif

  return type;
}

Solver* SolverFactory::createSolver(Options *options) {
  SolverType type = getDefaultSolverType();

  if(options == NULL) 
    options = Options::getRoot()->getSection("solver");
  
  string solver_option;
  options->get("type", solver_option, "");

  if(!solver_option.empty()) type = solver_option.c_str();

  return createSolver(type, options);
}

Solver* SolverFactory::createSolver(SolverType &type, Options *options) {
  if(options == NULL)
    options = Options::getRoot()->getSection("solver");
  
  if(!strcasecmp(type, SOLVERPVODE)) {
    return new PvodeSolver(options);
  } else if(!strcasecmp(type, SOLVERCVODE)) {
    return new CvodeSolver(options);
  } else if(!strcasecmp(type, SOLVERIDA)) {
    return new IdaSolver(options);
  } else if(!strcasecmp(type, SOLVERPETSC)) {
    return new PetscSolver(options); 
  } else if(!strcasecmp(type, SOLVERSLEPC)) {
    return new SlepcSolver(options);
  } else if(!strcasecmp(type, SOLVERKARNIADAKIS)) {
    return new KarniadakisSolver(options);
  } else if(!strcasecmp(type, SOLVERRK4)) {
    return new RK4Solver(options);
  } else if(!strcasecmp(type, SOLVEREULER)) {
    return new EulerSolver(options);
  } else if(!strcasecmp(type, SOLVERRK3SSP)) {
    return new RK3SSP(options);
  } else if(!strcasecmp(type, SOLVERPOWER)) {
    return new PowerSolver;
  } else if(!strcasecmp(type, SOLVERARKODE)){
    return new ArkodeSolver(options);
  } else if(!strcasecmp(type, SOLVERIMEXBDF2)) {
    return new IMEXBDF2(options);
  } else if(!strcasecmp(type, SOLVERSNES)) {
    return new SNESSolver(options);
  } else if(!strcasecmp(type, SOLVERRKGENERIC)){
    return new RKGenericSolver(options);
  }

  // Need to throw an error saying 'Supplied option "type"' was not found
  throw BoutException("No such solver exists in this build, type: %s", type);
}
