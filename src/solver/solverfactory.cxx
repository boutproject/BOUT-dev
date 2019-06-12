#include "bout/solverfactory.hxx"

#include "impls/arkode/arkode.hxx"
#include "impls/cvode/cvode.hxx"
#include "impls/euler/euler.hxx"
#include "impls/ida/ida.hxx"
#include "impls/imex-bdf2/imex-bdf2.hxx"
#include "impls/karniadakis/karniadakis.hxx"
#include "impls/petsc/petsc.hxx"
#include "impls/power/power.hxx"
#include "impls/pvode/pvode.hxx"
#include "impls/rk3-ssp/rk3-ssp.hxx"
#include "impls/rk4/rk4.hxx"
#include "impls/rkgeneric/rkgeneric.hxx"
#include "impls/slepc/slepc.hxx"
#include "impls/snes/snes.hxx"
#include "impls/split-rk/split-rk.hxx"

SolverFactory* SolverFactory::instance = nullptr;

SolverFactory* SolverFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new SolverFactory();
  }
  return instance;
}

SolverType SolverFactory::getDefaultSolverType() {
  return
#if defined BOUT_HAS_CVODE
      SOLVERCVODE;
#elif defined BOUT_HAS_IDA
      SOLVERIDA;
#else
      SOLVERPVODE;
#endif
}

Solver* SolverFactory::createSolver(Options* options) {
  if (options == nullptr) {
    options = Options::getRoot()->getSection("solver");
  }

  auto type = (*options)["type"].withDefault(getDefaultSolverType());

  return createSolver(type, options);
}
