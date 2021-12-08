#ifndef FAKESOLVER_H
#define FAKESOLVER_H

#include "gtest/gtest.h"

#include "bout/solver.hxx"

#include <algorithm>
#include <string>
#include <vector>

class FakeSolver : public Solver {
public:
  FakeSolver(Options* options) : Solver(options) { has_constraints = true; }
  ~FakeSolver() = default;

  int run() override {
    run_called = true;
    if ((*options)["throw_run"].withDefault(false)) {
      throw BoutException("Deliberate exception in FakeSolver::run");
    }
    return (*options)["fail_run"].withDefault(0);
  }
  bool run_called{false};

  int init(int nout, BoutReal tstep) override {
    init_called = true;
    if (Solver::init(nout, tstep)) {
      return 1;
    }
    return (*options)["fail_init"].withDefault(0);
  }
  bool init_called{false};

  void changeHasConstraints(bool new_value) { has_constraints = new_value; }

  auto listField2DNames() -> std::vector<std::string> {
    std::vector<std::string> result{};
    std::transform(begin(f2d), end(f2d), std::back_inserter(result),
                   [](const VarStr<Field2D>& f) { return f.name; });
    return result;
  }

  auto listField3DNames() -> std::vector<std::string> {
    std::vector<std::string> result{};
    std::transform(begin(f3d), end(f3d), std::back_inserter(result),
                   [](const VarStr<Field3D>& f) { return f.name; });
    return result;
  }

  auto listVector2DNames() -> std::vector<std::string> {
    std::vector<std::string> result{};
    std::transform(begin(v2d), end(v2d), std::back_inserter(result),
                   [](const VarStr<Vector2D>& f) { return f.name; });
    return result;
  }

  auto listVector3DNames() -> std::vector<std::string> {
    std::vector<std::string> result{};
    std::transform(begin(v3d), end(v3d), std::back_inserter(result),
                   [](const VarStr<Vector3D>& f) { return f.name; });
    return result;
  }

  // Shims for protected functions
  auto getMaxTimestepShim() const -> BoutReal { return max_dt; }
  using Solver::getLocalN;
  using Solver::hasPreconditioner;
  using Solver::runPreconditioner;
  using Solver::globalIndex;
  using Solver::getMonitors;
  using Solver::call_monitors;
  using Solver::call_timestep_monitors;
  using Solver::hasJacobian;
  using Solver::runJacobian;
};

#endif // FAKESOLVER_H
