#include "solverfactory.h"

#include "impls/cvode/cvode.h"
#include "impls/petsc/petsc.h"
#include "impls/ida/ida.h"
#include "impls/pvode/pvode.h"

#include <stdexcept>
#include "boutexception.h"

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
  #elif defined BOUT_HAS_PETSC
    type = SOLVERPETSC;
  #else
    type = SOLVERPVODE;
  #endif
  
  return type;
}

Solver* SolverFactory::createSolver() {
  SolverType type = getDefaultSolverType();
  
  options.setSection(NULL);
  const char* solver_option = options.getString("solver_type");
  
  if(solver_option) type = solver_option;
  
  
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
    } else if(!strcasecmp(type, SOLVERPETSC)) {
      return new PetscSolver;
    }  

  // Need to throw an error saying 'Supplied option "type"' was not found
    throw BoutException("No such solver exists in this build, type: %s", type);
}
