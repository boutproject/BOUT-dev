#ifndef FAKESOLVER_H
#define FAKESOLVER_H

#include "gtest/gtest.h"

#include "bout/solver.hxx"
#include "bout/solverfactory.hxx"

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
  auto getLocalNShim() -> int { return getLocalN(); }
  auto haveUserPreconShim() -> bool { return hasPreconditioner(); }
  auto runPreconShim(BoutReal t, BoutReal gamma, BoutReal delta) -> int {
    return runPreconditioner(t, gamma, delta);
  }
  auto globalIndexShim(int local_start) -> Field3D { return globalIndex(local_start); }
  auto getMonitorsShim() const -> const std::list<Monitor*>& { return getMonitors(); }
  auto callMonitorsShim(BoutReal simtime, int iter, int NOUT) -> int {
    return call_monitors(simtime, iter, NOUT);
  }
  auto callTimestepMonitorsShim(BoutReal simtime, BoutReal lastdt) -> int {
    return call_timestep_monitors(simtime, lastdt);
  }
  using Solver::hasJacobian;
  using Solver::runJacobian;
};

#endif // FAKESOLVER_H
