#include "bout/solverfactory.hxx"

SolverFactory* SolverFactory::instance = nullptr;

SolverFactory *SolverFactory::getInstance() {
  if (instance == nullptr) {
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

Solver *SolverFactory::createSolver(Options *options) {
  SolverType type = getDefaultSolverType();

  if (options == nullptr) {
    options = Options::getRoot()->getSection("solver");
  }

  std::string solver_option;
  options->get("type", solver_option, "");

  if (!solver_option.empty()) {
    type = solver_option.c_str();
  }

  return createSolver(type, options);
}
