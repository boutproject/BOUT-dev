#include "solverfactory.hxx"

#include "impls/cvode/cvode.hxx"
#include "impls/petsc-3.1/petsc-3.1.hxx"
#include "impls/petsc-3.2/petsc-3.2.hxx"
#include "impls/petsc-3.3/petsc-3.3.hxx"
#include "impls/petsc/petsc.hxx"
#include "impls/ida/ida.hxx"
#include "impls/pvode/pvode.hxx"
#include "impls/karniadakis/karniadakis.hxx"
#include "impls/rk4/rk4.hxx"
#include "impls/euler/euler.hxx"

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
    return new CvodeSolver;
  } else if(!strcasecmp(type, SOLVERIDA)) {
    return new IdaSolver;
  } else if(!strcasecmp(type, SOLVERPETSC)) {
    return new PetscSolver;
  } else if(!strcasecmp(type, SOLVERKARNIADAKIS)) {
    return new KarniadakisSolver;
  } else if(!strcasecmp(type, SOLVERRK4)) {
    return new RK4Solver;
  } else if(!strcasecmp(type, SOLVEREULER)) {
    return new EulerSolver;
  }

  // Need to throw an error saying 'Supplied option "type"' was not found
  throw BoutException("No such solver exists in this build, type: %s", type);
}
