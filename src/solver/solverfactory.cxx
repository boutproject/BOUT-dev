#include "solverfactory.hxx"

#include "impls/cvode/cvode.hxx"
#include "impls/petsc-3.1/petsc-3.1.hxx"
#include "impls/petsc/petsc.hxx"
#include "impls/ida/ida.hxx"
#include "impls/pvode/pvode.hxx"
#include "impls/karniadakis/karniadakis.hxx"
#include "impls/rk4/rk4.hxx"

#include <boutexception.hxx>

SolverFactory* SolverFactory::instance = NULL;

SolverFactory* SolverFactory::getInstance()
{
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

Solver* SolverFactory::createSolver() {
  SolverType type = getDefaultSolverType();

  Options *options = Options::getRoot();
  options = options->getSection("solver");
  string solver_option;
/*  options.get("solver_type", solver_option, type);
  string solver_option;*/
  options->get("type", solver_option, "");

  if(!solver_option.empty()) type = solver_option.c_str();

  return createSolver(type);
}

Solver* SolverFactory::createSolver(SolverType &type)
{
  if(!strcasecmp(type, SOLVERPVODE)) {
    return new PvodeSolver;
  } else if(!strcasecmp(type, SOLVERCVODE)) {
    return new CvodeSolver;
  } else if(!strcasecmp(type, SOLVERIDA)) {
    return new IdaSolver;
  } else if(!strcasecmp(type, SOLVERPETSC31)) {
    return new Petsc31Solver;
  } else if(!strcasecmp(type, SOLVERPETSC)) {
    return new PetscSolver;
  } else if(!strcasecmp(type, SOLVERKARNIADAKIS)) {
    return new KarniadakisSolver;
  } else if(!strcasecmp(type, SOLVERRK4)) {
    return new RK4Solver;
  }

  // Need to throw an error saying 'Supplied option "type"' was not found
  throw BoutException("No such solver exists in this build, type: %s", type);
}
