
class SolverFactory;

#ifndef __SOLVER_FACTORY_H__
#define __SOLVER_FACTORY_H__

#include "bout/generic_factory.hxx"
#include "bout/solver.hxx"

#include <functional>
#include <string>
#include <iostream>

class Options;

/// Factory specialisation for Solvers
///
/// TODO: Inherits from Factory<> in order to keep the same interface during
/// transition period. After next major version, replace whole class with:
///
///     using SolverFactory = Factory<Solver, SolverType,
///                                   std::function<Solver*(Options*)>>;
class SolverFactory : public Factory<Solver, std::function<Solver*(Options*)>> {
 public:
  SolverType getDefaultSolverType();
  
  static SolverFactory *getInstance();

  Solver* createSolver(Options *options = nullptr);
  Solver* createSolver(SolverType &name) {
    return create(name, Options::getRoot()->getSection("solver"));
  }
  Solver* createSolver(SolverType &name, Options *options) {
    return create(name, options);
  }
private:
  SolverFactory() {}
  static SolverFactory* instance;
};

/// Specialisation of Factory registration helper class
template<typename DerivedType>
class RegisterInFactory<Solver, DerivedType> {
public:
  RegisterInFactory(const std::string &name) {
    SolverFactory::getInstance()->add(
        name, [](Options *options) -> Solver * { return new DerivedType(options); });
  }
};

/// Simpler name for Factory registration helper class
///
/// Usage:
///
///     #include <bout/solverfactory.hxx>
///     namespace {
///     RegisterSolver<MySolver> registersolvermine("mysolver");
///     }
template<typename DerivedType>
using RegisterSolver = RegisterInFactory<Solver, DerivedType>;

#endif // __SOLVER_FACTORY_H__

