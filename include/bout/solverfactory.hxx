#ifndef __SOLVER_FACTORY_H__
#define __SOLVER_FACTORY_H__

#include "bout/generic_factory.hxx"
#include "bout/solver.hxx"

#include <string>

template<>
struct StandardFactoryTraits<Solver> {
  static constexpr auto type_name = "Solver";
  static constexpr auto section_name = "solver";
  static std::string getDefaultType();
};

using SolverFactory = StandardFactory<Solver>;

/// Simpler name for Factory registration helper class
///
/// Usage:
///
///     #include <bout/solverfactory.hxx>
///     namespace {
///     RegisterSolver<MySolver> registersolvermine("mysolver");
///     }
template <typename DerivedType>
using RegisterSolver = RegisterInStandardFactory<Solver, DerivedType>;

#endif // __SOLVER_FACTORY_H__
