#include "gtest/gtest.h"

#include "../../../src/solver/impls/cvode/cvode_step_limit.hxx"

#include <cmath>
#include <vector>

using bout::cvode::cutFactor;
using bout::cvode::IterationFit;
using bout::cvode::LinearSolveRecord;
using bout::cvode::StepLimitAction;
using bout::cvode::StepLimiter;
using bout::cvode::StepLimitParams;

namespace {
LinearSolveRecord converged(BoutReal h, BoutReal gamma, int iters) {
  LinearSolveRecord solve;
  solve.h = h;
  solve.gamma = gamma;
  solve.tol = 1e-3;
  solve.r0 = 1.0;
  solve.rend = 1e-4;
  solve.iters = iters;
  solve.failed = false;
  return solve;
}

LinearSolveRecord failedSolve(BoutReal h, BoutReal r0, BoutReal rend, BoutReal tol,
                              int iters) {
  LinearSolveRecord solve;
  solve.h = h;
  solve.gamma = h;
  solve.tol = tol;
  solve.r0 = r0;
  solve.rend = rend;
  solve.iters = iters;
  solve.failed = true;
  return solve;
}
} // namespace

TEST(CvodeStepLimitTest, CutFactorFromResidualReduction) {
  // 50 iterations took the residual from 1 to 1e-2; 1e-4 needs 100.
  // With p = 1, safety 0.7, budget 50: factor = 0.7 * 50 / 100 = 0.35
  const auto solve = failedSolve(100.0, 1.0, 1e-2, 1e-4, 50);
  EXPECT_NEAR(cutFactor(solve, 50.0, 1.0, 0.7), 0.35, 1e-12);
  // With p = 2 the step needs to shrink less: sqrt(0.35)
  EXPECT_NEAR(cutFactor(solve, 50.0, 2.0, 0.7), std::sqrt(0.35), 1e-12);
}

TEST(CvodeStepLimitTest, CutFactorClamped) {
  // Nearly converged, with iterations growing as gamma^2: the estimate
  // sqrt(0.7 * 50 / 50.5) = 0.83 exceeds 0.8, so it is clamped
  EXPECT_DOUBLE_EQ(cutFactor(failedSolve(1.0, 1.0, 1.1e-4, 1e-4, 50), 50.0, 2.0, 0.7),
                   0.8);
  // Barely reduced: the estimate is tiny, so it is clamped at 0.2
  EXPECT_DOUBLE_EQ(cutFactor(failedSolve(1.0, 1.0, 0.99, 1e-4, 50), 50.0, 1.0, 0.7), 0.2);
}

TEST(CvodeStepLimitTest, CutFactorWithoutResidualHistory) {
  // No reduction at all
  EXPECT_DOUBLE_EQ(cutFactor(failedSolve(1.0, 1.0, 1.0, 1e-4, 50), 50.0, 1.0, 0.7), 0.2);
  // Initial residual unknown (left preconditioning)
  EXPECT_DOUBLE_EQ(cutFactor(failedSolve(1.0, -1.0, 0.5, 1e-4, 50), 50.0, 1.0, 0.7), 0.5);
}

TEST(CvodeStepLimitTest, IterationFitDefaultsAndSlope) {
  IterationFit fit;
  EXPECT_DOUBLE_EQ(fit.slope(), 1.0); // no points

  // Too little spread in gamma
  for (int i = 0; i < 10; ++i) {
    fit.add(1.0, 10);
  }
  EXPECT_DOUBLE_EQ(fit.slope(), 1.0);

  // iterations = 2 * gamma^1.5 over a decade of gamma
  IterationFit fit15;
  for (int i = 0; i < 20; ++i) {
    const BoutReal gamma = 10.0 * std::pow(10.0, i / 19.0);
    fit15.add(gamma, static_cast<int>(std::lround(2.0 * std::pow(gamma, 1.5))));
  }
  EXPECT_NEAR(fit15.slope(), 1.5, 0.01);

  // Solves under 5 iterations are ignored
  IterationFit small;
  for (int i = 0; i < 20; ++i) {
    small.add(1.0 + i, 3);
  }
  EXPECT_DOUBLE_EQ(small.slope(), 1.0);
}

TEST(CvodeStepLimitTest, IterationFitClamped) {
  IterationFit steep;
  for (int i = 0; i < 20; ++i) {
    const BoutReal gamma = 1.0 + 0.1 * i;
    steep.add(gamma, static_cast<int>(std::lround(5.0 * std::pow(gamma, 4.0))));
  }
  EXPECT_DOUBLE_EQ(steep.slope(), 2.0);
}

TEST(CvodeStepLimitTest, NoLimitWithoutFailures) {
  StepLimiter limiter;
  const auto action = limiter.update({converged(10.0, 10.0, 10)}, 10.0, 50.0);
  EXPECT_EQ(action, StepLimitAction::none);
  EXPECT_DOUBLE_EQ(limiter.getLimit(), 0.0);
}

TEST(CvodeStepLimitTest, CutSetsLimitAndWall) {
  StepLimiter limiter;
  // Failed at h = 100, then CVODE succeeded at 25
  const std::vector<LinearSolveRecord> solves{failedSolve(100.0, 1.0, 1e-2, 1e-4, 50),
                                              converged(25.0, 25.0, 20)};
  EXPECT_EQ(limiter.update(solves, 25.0, 50.0), StepLimitAction::cut);
  EXPECT_NEAR(limiter.getLimit(), 35.0, 1e-9); // 0.35 * 100
  EXPECT_DOUBLE_EQ(limiter.getWall(), 100.0);
}

TEST(CvodeStepLimitTest, HoldNearIterationBudget) {
  StepLimiter limiter;
  limiter.update({failedSolve(100.0, 1.0, 1e-2, 1e-4, 50)}, 25.0, 50.0);
  // 45 iterations > 0.8 * 50: hold at the step just taken
  EXPECT_EQ(limiter.update({converged(30.0, 30.0, 45)}, 30.0, 50.0),
            StepLimitAction::hold);
  EXPECT_DOUBLE_EQ(limiter.getLimit(), 30.0);
}

TEST(CvodeStepLimitTest, RelaxCappedByWallWhichGrows) {
  StepLimitParams params;
  StepLimiter limiter(params);
  limiter.update({failedSolve(100.0, 1.0, 1e-2, 1e-4, 50)}, 25.0, 50.0); // limit 35
  BoutReal wall = 100.0;
  BoutReal limit = 35.0;
  for (int step = 0; step < 6; ++step) {
    EXPECT_EQ(limiter.update({converged(limit, limit, 10)}, limit, 50.0),
              StepLimitAction::relax);
    wall *= params.wall_relax;
    limit = std::min(limit * params.relax, params.wall_margin * wall);
    EXPECT_NEAR(limiter.getWall(), wall, 1e-9);
    EXPECT_NEAR(limiter.getLimit(), limit, 1e-9);
  }
  // The limit has reached the wall's margin and only grows with the wall
  EXPECT_NEAR(limiter.getLimit(), params.wall_margin * limiter.getWall(), 1e-9);
}

TEST(CvodeStepLimitTest, ReleaseFarAboveSteps) {
  StepLimitParams params;
  params.wall_margin = 0.0; // no wall, so the limit can relax freely
  StepLimiter limiter(params);
  limiter.update({failedSolve(100.0, 1.0, 1e-2, 1e-4, 50)}, 25.0, 50.0); // limit 35
  // Steps stay at 1: the limit relaxes until it exceeds 10 * h, then is dropped
  StepLimitAction action = StepLimitAction::relax;
  int steps = 0;
  while (action == StepLimitAction::relax and steps < 100) {
    action = limiter.update({converged(1.0, 1.0, 10)}, 1.0, 50.0);
    ++steps;
  }
  EXPECT_EQ(action, StepLimitAction::release);
  EXPECT_DOUBLE_EQ(limiter.getLimit(), 0.0);
  EXPECT_DOUBLE_EQ(limiter.getWall(), 0.0);
}
