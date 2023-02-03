/**************************************************************************
 * Generic Adams Bashforth multistep scheme
 *
 * Always available, since doesn't depend on external library
 *
 **************************************************************************
 * Written by D Dickinson 2019
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

class AdamsBashforthSolver;

#ifndef __ADAMSBASHFORTH_SOLVER_H__
#define __ADAMSBASHFORTH_SOLVER_H__

#include <bout/bout_types.hxx>
#include <bout/solver.hxx>

#include <deque>

namespace {
RegisterSolver<AdamsBashforthSolver> registersolveradamsbashforth("adams-bashforth");
}

class AdamsBashforthSolver : public Solver {
public:
  explicit AdamsBashforthSolver(Options* options = nullptr);
  ~AdamsBashforthSolver() = default;

  void resetInternalFields() override;

  // Utilities only used by the CTU bracket approach
  void setMaxTimestep(BoutReal dt) override;
  BoutReal getCurrentTimestep() override { return timestep; }

  // Setup solver and scheme
  int init() override;

  // Actually evolve
  int run() override;

private:
  // Take a single timestep of specified order. If adaptive also calculates
  // and returns an error estimate.
  BoutReal take_step(BoutReal timeIn, BoutReal dt, int order, Array<BoutReal>& current,
                     Array<BoutReal>& result);

  // Holds the current/next state
  Array<BoutReal> state, nextState;

  // State history - we use deque's to make it easy to add/remove from
  // either end.  Whilst this looks like it might be expensive for
  // states (rather than say using std::rotate with a std::vector to
  // just move things around) we're relying on the Array store making
  // it cheap to get a "new" array.
  std::deque<Array<BoutReal>> history; // History of d state/dt values
  std::deque<BoutReal> times;          // Times at which above states calculated

  // Inputs
  /// Tolerances for adaptive timestepping
  BoutReal atol, rtol;
  /// Maximum number of internal steps between outputs
  int mxstep;
  /// Adapt timestep?
  bool adaptive;
  /// Adapt order?
  bool adaptive_order;
  /// If true and adaptive the solution used is the more accurate one.
  bool followHighOrder;
  /// Factor we scale timestep estimate by when adapting.
  BoutReal dtFac;
  /// The maximum order scheme to use.
  int maximum_order;
  /// Maximum timestep
  BoutReal max_timestep;
  /// The internal timestep
  BoutReal timestep;

  // Internal vars
  int current_order; // The current order of the scheme
  int nlocal, neq;   // Number of variables on local processor and in total
};

#endif // __ADAMSBASHFORTH_SOLVER_H__
