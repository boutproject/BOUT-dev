/**************************************************************************
 * A maximum step size for CVODE learned from its Krylov linear solves
 *
 * With right preconditioning (or FGMRES) CVODE's error estimate no longer
 * includes the preconditioner, so it can take steps for which a weak
 * preconditioner leaves GMRES unable to converge within maxl iterations.
 * CVODE reacts to each such failure by cutting the step, then grows it back
 * to the same size and fails again. StepLimiter keeps a memory instead:
 *
 *  - after a failed linear solve it cuts the limit by the residual reduction
 *    GMRES reached, extrapolated to the iterations it would have needed;
 *  - it holds the limit when a step's solves come close to maxl;
 *  - otherwise it relaxes the limit, but keeps it below a margin of the
 *    smallest step that recently failed (the "wall"), which itself grows
 *    slowly so the limit can follow an easier problem;
 *  - it releases the limit once it is far above the steps being taken.
 *
 * Everything here is pure arithmetic on the recorded solves, so it can be
 * unit tested without CVODE.
 **************************************************************************/

#ifndef BOUT_CVODE_STEP_LIMIT_H
#define BOUT_CVODE_STEP_LIMIT_H

#include "bout/bout_types.hxx"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace bout {
namespace cvode {

/// One Krylov linear solve, as seen by the step control
struct LinearSolveRecord {
  BoutReal h{0.0};     ///< Step size being attempted
  BoutReal gamma{0.0}; ///< gamma in (I - gamma J) x = b
  BoutReal tol{0.0};   ///< Tolerance passed to the linear solver
  BoutReal r0{-1.0};   ///< Initial residual in the solver's stop norm; < 0 if unknown
  BoutReal rend{0.0};  ///< Final residual in the same norm
  int iters{0};        ///< Iterations taken
  bool failed{
      false}; ///< Not converged (including stopped at maxl with a reduced residual)
};

/// Tuning parameters (the options cvode_lsc_*)
struct StepLimitParams {
  BoutReal hold{0.8};         ///< Hold when iterations exceed hold * budget
  BoutReal relax{1.3};        ///< Growth of the limit per quiet step
  BoutReal safety{0.7};       ///< Target fraction of the budget after a cut
  BoutReal wall_margin{0.85}; ///< Limit <= wall_margin * wall; 0 turns the wall off
  BoutReal wall_relax{1.01};  ///< Growth of the wall per step without a failure
};

enum class StepLimitAction { none, cut, hold, relax, release };

/// Least-squares slope p of log(iterations) against log(gamma) over recent
/// converged solves, modelling iterations ~ gamma^p. Returns 1 when there are
/// too few points or too little spread in gamma, and clamps to [0.5, 2].
class IterationFit {
public:
  static constexpr std::size_t capacity = 50;
  static constexpr int min_iterations = 5;

  void add(BoutReal gamma, int iterations) {
    if (gamma <= 0.0 or iterations < min_iterations) {
      return;
    }
    const std::pair<BoutReal, BoutReal> point{
        std::log(gamma), std::log(static_cast<BoutReal>(iterations))};
    if (points.size() < capacity) {
      points.push_back(point);
    } else {
      points[next] = point;
    }
    next = (next + 1) % capacity;
  }

  BoutReal slope() const {
    if (points.size() < 5) {
      return 1.0;
    }
    BoutReal mean_x = 0.0;
    BoutReal mean_y = 0.0;
    for (const auto& point : points) {
      mean_x += point.first;
      mean_y += point.second;
    }
    const auto n = static_cast<BoutReal>(points.size());
    mean_x /= n;
    mean_y /= n;
    BoutReal sxx = 0.0;
    BoutReal sxy = 0.0;
    for (const auto& point : points) {
      sxx += (point.first - mean_x) * (point.first - mean_x);
      sxy += (point.first - mean_x) * (point.second - mean_y);
    }
    if (sxx / n < 0.01) { // gamma spread under about 10 %
      return 1.0;
    }
    return std::clamp(sxy / sxx, 0.5, 2.0);
  }

private:
  std::vector<std::pair<BoutReal, BoutReal>> points;
  std::size_t next{0};
};

/// Factor by which to cut the step after a failed solve.
///
/// If GMRES reduced the residual from r0 to rend in `iters` iterations and
/// needed `tol`, extrapolating linearly in log(residual) it would have needed
///   k = iters * ln(tol / r0) / ln(rend / r0)
/// iterations. With iterations ~ gamma^p, the step that fits
/// safety * budget iterations is h * (safety * budget / k)^(1/p).
/// The factor is clamped to [0.2, 0.8]. Without a usable residual history it
/// is 0.2, or 0.5 if r0 is unknown (left preconditioning).
inline BoutReal cutFactor(const LinearSolveRecord& solve, BoutReal budget, BoutReal p,
                          BoutReal safety) {
  BoutReal factor = 0.2;
  if (solve.r0 > 0.0 and solve.rend > 0.0 and solve.rend < solve.r0
      and solve.tol < solve.r0 and solve.iters > 0) {
    const BoutReal needed =
        solve.iters * std::log(solve.tol / solve.r0) / std::log(solve.rend / solve.r0);
    factor = std::pow(safety * budget / needed, 1.0 / p);
  } else if (solve.r0 < 0.0) {
    factor = 0.5;
  }
  return std::clamp(factor, 0.2, 0.8);
}

/// The learned step limit
class StepLimiter {
public:
  explicit StepLimiter(StepLimitParams params = {}) : params(params) {}

  /// Update after an internal step from the solves made during it.
  /// @param solves  Every linear solve since the previous step, including
  ///                those of failed attempts
  /// @param h_last  Size of the step just taken
  /// @param budget  Iterations GMRES may take (maxl)
  StepLimitAction update(const std::vector<LinearSolveRecord>& solves, BoutReal h_last,
                         BoutReal budget) {
    constexpr BoutReal inf = std::numeric_limits<BoutReal>::infinity();
    int max_iters = 0;
    bool any_failed = false;
    for (const auto& solve : solves) {
      if (solve.failed) {
        any_failed = true;
      } else {
        fit.add(solve.gamma, solve.iters);
      }
      max_iters = std::max(max_iters, solve.iters);
    }
    p = fit.slope();
    h_last = std::abs(h_last);
    BoutReal hlim = (limit > 0.0) ? limit : inf;
    StepLimitAction action = StepLimitAction::none;

    if (any_failed) {
      for (const auto& solve : solves) {
        if (not solve.failed) {
          continue;
        }
        const BoutReal h_fail = std::abs(solve.h);
        hlim = std::min(hlim, cutFactor(solve, budget, p, params.safety) * h_fail);
        if (params.wall_margin > 0.0) {
          wall = (wall > 0.0) ? std::min(wall, h_fail) : h_fail;
        }
      }
      action = StepLimitAction::cut;
    } else {
      if (wall > 0.0) {
        wall *= params.wall_relax;
      }
      if (max_iters > params.hold * budget) {
        hlim = std::min(hlim, h_last);
        action = StepLimitAction::hold;
      } else if (hlim < inf) {
        hlim *= params.relax;
        if (wall > 0.0) {
          hlim = std::min(hlim, params.wall_margin * wall);
        }
        action = StepLimitAction::relax;
        if (hlim > 10.0 * h_last) {
          hlim = inf;
          wall = 0.0;
          action = StepLimitAction::release;
        }
      }
    }
    if (action != StepLimitAction::none) {
      limit = (hlim < inf) ? hlim : 0.0;
    }
    return action;
  }

  /// Current limit; 0 means none
  BoutReal getLimit() const { return limit; }
  /// Smallest recently failing step; 0 if none
  BoutReal getWall() const { return wall; }
  /// Last fitted exponent of iterations against gamma
  BoutReal getExponent() const { return p; }

private:
  StepLimitParams params;
  IterationFit fit;
  BoutReal limit{0.0};
  BoutReal wall{0.0};
  BoutReal p{1.0};
};

} // namespace cvode
} // namespace bout

#endif // BOUT_CVODE_STEP_LIMIT_H
