#include "bout/solverfactory.hxx"

#include "impls/arkode/arkode.hxx"
#include "impls/cvode/cvode.hxx"
#include "impls/euler/euler.hxx"
#include "impls/ida/ida.hxx"
#include "impls/imex-bdf2/imex-bdf2.hxx"
#include "impls/karniadakis/karniadakis.hxx"
#include "impls/petsc-3.1/petsc-3.1.hxx"
#include "impls/petsc-3.2/petsc-3.2.hxx"
#include "impls/petsc-3.3/petsc-3.3.hxx"
#include "impls/petsc-3.4/petsc-3.4.hxx"
#include "impls/petsc-3.5/petsc-3.5.hxx"
#include "impls/petsc/petsc.hxx"
#include "impls/power/power.hxx"
#include "impls/pvode/pvode.hxx"
#include "impls/rk3-ssp/rk3-ssp.hxx"
#include "impls/rk4/rk4.hxx"
#include "impls/rkgeneric/rkgeneric.hxx"
#include "impls/slepc-3.4/slepc-3.4.hxx"
#include "impls/snes/snes.hxx"

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
