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

#include "mpi.h"

#include <bout/solver.hxx>
#include <bout/solverfactory.hxx>
#include <bout_types.hxx>

#include <deque>

namespace {
RegisterSolver<AdamsBashforthSolver> registersolveradamsbashforth("adams-bashforth");
}

struct AdamsBashforthHelper {
  BoutReal lagrange_at_position_denominator(const std::deque<BoutReal>& grid,
                                            const int position, const int order) const {
    AUTO_TRACE();
    ASSERT2(position < order);
    ASSERT2(order <= grid.size());

    const auto xj = grid[position];

    BoutReal result = 1.0;
    for (int i = 0; i < order; i++) {
      result /= (i != position) ? (xj - grid[i]) : 1.0;
    }
    return result;
  };

  BoutReal lagrange_at_position_numerator(const BoutReal varX,
                                          const std::deque<BoutReal>& grid,
                                          const int position, const int order) const {
    AUTO_TRACE();
    ASSERT2(position < order);
    ASSERT2(order <= grid.size());
    BoutReal result = 1.0;
    for (int i = 0; i < order; i++) {
      result *= (i != position) ? (varX - grid[i]) : 1.0;
    }
    return result;

    // // Above could be rewritten as following but not sure this is more readable and
    // // the floating comparison that seems to be required is less nice. Possibly that
    // // we could use std::iota(grid.size()) as the iterator args and then use that value
    // // to index grid (which we'd have to capture) instead but again that seems more
    // complex.
    // const auto tmp = grid[position];
    // return std::accumulate(std::begin(grid), std::end(grid), 1.0,
    //                        [varX, tmp](BoutReal current, BoutReal gridVal) {
    //                          return current * ((gridVal != tmp)? (varX - gridVal) :
    //                          1.0);
    //                        });
  };

  // Integrate using newton-cotes 9 rule
  BoutReal integrate_lagrange_curve_nc9(const BoutReal theStart, const BoutReal theEnd,
                                        const std::deque<BoutReal>& points,
                                        const int position, const int order) const {
    AUTO_TRACE();
    constexpr int size = 9;
    constexpr BoutReal fac = 4.0 / 14175.0;
    constexpr std::array<BoutReal, size> facs{989.0 * fac,   5888.0 * fac,  -928.0 * fac,
                                              10496.0 * fac, -4540.0 * fac, 10496.0 * fac,
                                              -928.0 * fac,  5888.0 * fac,  989.0 * fac};
    constexpr BoutReal stepFac = 1.0 / (size - 1.0);
    const BoutReal stepSize = (theEnd - theStart) * stepFac;

    BoutReal result{0.0};
    for (int i = 0; i < size; i++) {
      result += facs[i] * lagrange_at_position_numerator(theStart + i * stepSize, points,
                                                         position, order);
    }
    return stepSize * result * lagrange_at_position_denominator(points, position, order);
  };

  // Integrate using newton-cotes 8 rule
  BoutReal integrate_lagrange_curve_nc8(const BoutReal theStart, const BoutReal theEnd,
                                        const std::deque<BoutReal>& points,
                                        const int position, const int order) const {
    AUTO_TRACE();
    constexpr int size = 8;
    constexpr BoutReal fac = 7.0 / 17280.0;
    constexpr std::array<BoutReal, size> facs{751.0 * fac,  3577.0 * fac, 1323.0 * fac,
                                              2989.0 * fac, 2989.0 * fac, 1323.0 * fac,
                                              3577.0 * fac, 751.0 * fac};
    constexpr BoutReal stepFac = 1.0 / (size - 1.0);
    const BoutReal stepSize = (theEnd - theStart) * stepFac;

    BoutReal result{0.0};
    for (int i = 0; i < size; i++) {
      result += facs[i] * lagrange_at_position_numerator(theStart + i * stepSize, points,
                                                         position, order);
    }
    return stepSize * result * lagrange_at_position_denominator(points, position, order);
  };

  // Integrate using newton-cotes 7 rule
  BoutReal integrate_lagrange_curve_nc7(const BoutReal theStart, const BoutReal theEnd,
                                        const std::deque<BoutReal>& points,
                                        const int position, const int order) const {
    AUTO_TRACE();
    constexpr int size = 7;
    constexpr BoutReal fac = 1.0 / 140.0;
    constexpr std::array<BoutReal, size> facs{41.0 * fac,  216.0 * fac, 27.0 * fac,
                                              272.0 * fac, 27.0 * fac,  216.0 * fac,
                                              41.0 * fac};
    constexpr BoutReal stepFac = 1.0 / (size - 1.0);
    const BoutReal stepSize = (theEnd - theStart) * stepFac;

    BoutReal result{0.0};
    for (int i = 0; i < size; i++) {
      result += facs[i] * lagrange_at_position_numerator(theStart + i * stepSize, points,
                                                         position, order);
    }
    return stepSize * result * lagrange_at_position_denominator(points, position, order);
  };

  // Integrate using newton-cotes 6 rule
  BoutReal integrate_lagrange_curve_nc6(const BoutReal theStart, const BoutReal theEnd,
                                        const std::deque<BoutReal>& points,
                                        const int position, const int order) const {
    AUTO_TRACE();
    constexpr int size = 6;
    constexpr BoutReal fac = 5.0 / 288.0;
    constexpr std::array<BoutReal, size> facs{19.0 * fac, 75.0 * fac, 50.0 * fac,
                                              50.0 * fac, 75.0 * fac, 19.0 * fac};
    constexpr BoutReal stepFac = 1.0 / (size - 1.0);
    const BoutReal stepSize = (theEnd - theStart) * stepFac;

    BoutReal result{0};
    for (int i = 0; i < size; i++) {
      result += facs[i] * lagrange_at_position_numerator(theStart + i * stepSize, points,
                                                         position, order);
    }
    return stepSize * result * lagrange_at_position_denominator(points, position, order);
  };

  // Integrate using newton-cotes 5 rule (Boole)
  BoutReal integrate_lagrange_curve_nc5(const BoutReal theStart, const BoutReal theEnd,
                                        const std::deque<BoutReal>& points,
                                        const int position, const int order) const {
    AUTO_TRACE();
    constexpr int size = 5;
    constexpr BoutReal fac = 2.0 / 45.0;
    constexpr std::array<BoutReal, size> facs{7.0 * fac, 32.0 * fac, 12.0 * fac,
                                              32.0 * fac, 7.0 * fac};
    constexpr BoutReal stepFac = 1.0 / (size - 1.0);
    const BoutReal stepSize = (theEnd - theStart) * stepFac;

    BoutReal result{0.0};
    for (int i = 0; i < size; i++) {
      result += facs[i] * lagrange_at_position_numerator(theStart + i * stepSize, points,
                                                         position, order);
    }
    return stepSize * result * lagrange_at_position_denominator(points, position, order);
  };

  // Integrate using newton-cotes 4 rule (Simpson 3/8)
  BoutReal integrate_lagrange_curve_nc4(const BoutReal theStart, const BoutReal theEnd,
                                        const std::deque<BoutReal>& points,
                                        const int position, const int order) const {
    AUTO_TRACE();
    constexpr int size = 4;
    constexpr BoutReal fac = 3.0 / 8.0;
    constexpr std::array<BoutReal, size> facs{1.0 * fac, 3.0 * fac, 3.0 * fac, 1.0 * fac};
    constexpr BoutReal stepFac = 1.0 / (size - 1.0);
    const BoutReal stepSize = (theEnd - theStart) * stepFac;

    BoutReal result{0.0};
    for (int i = 0; i < size; i++) {
      result += facs[i] * lagrange_at_position_numerator(theStart + i * stepSize, points,
                                                         position, order);
    }
    return stepSize * result * lagrange_at_position_denominator(points, position, order);
  };

  // Integrate using newton-cotes 3 rule (Simpson)
  BoutReal integrate_lagrange_curve_nc3(const BoutReal theStart, const BoutReal theEnd,
                                        const std::deque<BoutReal>& points,
                                        const int position, const int order) const {
    AUTO_TRACE();
    constexpr int size = 3;
    constexpr BoutReal fac = 1.0 / 3.0;
    constexpr std::array<BoutReal, size> facs{1.0 * fac, 4.0 * fac, 1.0 * fac};
    constexpr BoutReal stepFac = 1.0 / (size - 1.0);
    const BoutReal stepSize = (theEnd - theStart) * stepFac;

    BoutReal result{0.0};
    for (int i = 0; i < size; i++) {
      result += facs[i] * lagrange_at_position_numerator(theStart + i * stepSize, points,
                                                         position, order);
    }
    return stepSize * result * lagrange_at_position_denominator(points, position, order);
  };

  // Integrate using newton-cotes 2 rule (Trap)
  BoutReal integrate_lagrange_curve_nc2(const BoutReal theStart, const BoutReal theEnd,
                                        const std::deque<BoutReal>& points,
                                        const int position, const int order) const {
    AUTO_TRACE();
    constexpr int size = 2;
    constexpr BoutReal fac = 1.0 / 2.0;
    constexpr std::array<BoutReal, size> facs{1.0 * fac, 1.0 * fac};
    constexpr BoutReal stepFac = 1.0 / (size - 1.0);
    const BoutReal stepSize = (theEnd - theStart) * stepFac;

    BoutReal result{0.0};
    for (int i = 0; i < size; i++) {
      result += facs[i] * lagrange_at_position_numerator(theStart + i * stepSize, points,
                                                         position, order);
    }
    return stepSize * result * lagrange_at_position_denominator(points, position, order);
  };

  // Integrate lagrange polynomial to find the coefficienst of the requested order (note
  // don't currently
  // request an order just try to work it out from number of points).
  BoutReal integrate_lagrange_curve(const BoutReal theStart, const BoutReal theEnd,
                                    const std::deque<BoutReal>& points,
                                    const int position, const int order) const {
    AUTO_TRACE();
    ASSERT2(order <= points.size());

    switch (order) {
    case 1:
      return integrate_lagrange_curve_nc2(theStart, theEnd, points, position, order);
    case 2:
      return integrate_lagrange_curve_nc3(theStart, theEnd, points, position, order);
    case 3:
      return integrate_lagrange_curve_nc4(theStart, theEnd, points, position, order);
    case 4:
      return integrate_lagrange_curve_nc5(theStart, theEnd, points, position, order);
    case 5:
      return integrate_lagrange_curve_nc6(theStart, theEnd, points, position, order);
    case 6:
      return integrate_lagrange_curve_nc7(theStart, theEnd, points, position, order);
    case 7:
      return integrate_lagrange_curve_nc8(theStart, theEnd, points, position, order);
    default:
      return integrate_lagrange_curve_nc9(theStart, theEnd, points, position, order);
    }
  };

  // Calculate the set of Adams-Bashforth coefficients required to get from t = points[0]
  // to t = nextPoint
  // at the requested order.
  std::vector<BoutReal>
  get_adams_bashforth_coefficients(const BoutReal nextPoint,
                                   const std::deque<BoutReal>& points,
                                   const int order) const {
    AUTO_TRACE();
    ASSERT2(order <= points.size());

    std::vector<BoutReal> result;

    for (int i = 0; i < order; i++) {
      result.emplace_back(
          integrate_lagrange_curve(points[0], nextPoint, points, i, order));
    };

    return result;
  };
};

class AdamsBashforthSolver : public Solver {
public:
  AdamsBashforthSolver(Options* options = nullptr);
  ~AdamsBashforthSolver() = default;

  void resetInternalFields() override;

  // Utilities only used by the CTU bracket approach
  void setMaxTimestep(BoutReal dt) override;
  BoutReal getCurrentTimestep() override { return timestep; }

  // Setup solver and scheme
  int init(int nout, BoutReal tstep) override;

  // Actually evolve
  int run() override;

private:
  // Take a single timestep of specified order. If adaptive also calculates
  // and returns an error estimate.
  BoutReal take_step(const BoutReal timeIn, const BoutReal dt, const int order,
                     Array<BoutReal>& current, Array<BoutReal>& result);

  // Finds the maximum absolute error, i.e. Max(Abs(stateApprox - stateAccurate))
  // over all processors.
  BoutReal get_error(const Array<BoutReal>& stateApprox,
                     const Array<BoutReal>& stateAccurate) const {
    AUTO_TRACE();
    BoutReal local_result = 0.0;
    BoutReal err = 0.0;

    for (int i = 0; i < nlocal; i++) {
      local_result = std::max(std::abs(stateAccurate[i] - stateApprox[i]), local_result);

      // The below is the more typical error calculation used in other solvers.
      // We prefer the above definition as it provides a way to get a reasonable
      // estimate of the limiting timestep.
      // local_result = std::max(std::abs(stateAccurate[i] -
      //                        stateApprox[i]) / (std::abs(stateAccurate[i]) +
      //                        std::abs(stateApprox[i]) + atol), local_result);
    }

    // Reduce over procs
    if (MPI_Allreduce(&local_result, &err, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get())) {
      throw BoutException("MPI_Allreduce failed");
    }
    return err;
  };

  // Coefficient calculator
  AdamsBashforthHelper coefficients_calculator;

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
  BoutReal atol, rtol;   // Tolerances for adaptive timestepping
  BoutReal max_timestep; // Maximum timestep
  int mxstep;            // Maximum number of internal steps between outputs
  bool adaptive;         // Adapt timestep?
  bool
      followHighOrder; // If true and adaptive the solution used is the more accurate one.
  BoutReal dtFac;      // Factor we scale timestep estimate by when adapting.
  int maximum_order;   // The maximum order scheme to use.
  BoutReal timestep;   // The internal timestep

  // Internal vars
  BoutReal out_timestep; // The output timestep
  int current_order;     // The current order of the scheme
  int nsteps;            // Number of output steps
  int nlocal, neq;       // Number of variables on local processor and in total
};

#endif // __ADAMSBASHFORTH_SOLVER_H__
