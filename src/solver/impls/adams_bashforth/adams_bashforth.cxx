#include "adams_bashforth.hxx"

#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <utils.hxx>

#include <output.hxx>

namespace {
BoutReal lagrange_at_position_denominator(const std::deque<BoutReal>& grid,
                                          const int position, const int order) {
  AUTO_TRACE();

  const auto xj = grid[position];

  BoutReal result = 1.0;
  for (int i = 0; i < order; i++) {
    result /= (i != position) ? (xj - grid[i]) : 1.0;
  }
  return result;
}

BoutReal lagrange_at_position_numerator(const BoutReal varX,
                                        const std::deque<BoutReal>& grid,
                                        const int position, const int order) {
  AUTO_TRACE();
  BoutReal result = 1.0;
  for (int i = 0; i < order; i++) {
    result *= (i != position) ? (varX - grid[i]) : 1.0;
  }
  return result;
}

template <std::size_t N>
BoutReal lagrange_interpolate(BoutReal start, BoutReal end,
                              const std::deque<BoutReal>& points, const int position,
                              const std::array<BoutReal, N>& facs) {
  const BoutReal stepSize = (end - start) / (N - 1.0);

  BoutReal result{0.0};
  for (std::size_t i = 0; i < N; i++) {
    result +=
        facs[i]
        * lagrange_at_position_numerator(start + i * stepSize, points, position, N - 1);
  }
  return stepSize * result * lagrange_at_position_denominator(points, position, N - 1);
}

// Integrate using newton-cotes 9 rule
BoutReal integrate_lagrange_curve_nc9(const BoutReal start, const BoutReal end,
                                      const std::deque<BoutReal>& points,
                                      const int position) {
  AUTO_TRACE();
  constexpr std::size_t size = 9;
  constexpr BoutReal fac = 4.0 / 14175.0;
  constexpr std::array<BoutReal, size> facs{989.0 * fac,   5888.0 * fac,  -928.0 * fac,
                                            10496.0 * fac, -4540.0 * fac, 10496.0 * fac,
                                            -928.0 * fac,  5888.0 * fac,  989.0 * fac};
  return lagrange_interpolate(start, end, points, position, facs);
}

// Integrate using newton-cotes 8 rule
BoutReal integrate_lagrange_curve_nc8(const BoutReal start, const BoutReal end,
                                      const std::deque<BoutReal>& points,
                                      const int position) {
  AUTO_TRACE();
  constexpr std::size_t size = 8;
  constexpr BoutReal fac = 7.0 / 17280.0;
  constexpr std::array<BoutReal, size> facs{751.0 * fac,  3577.0 * fac, 1323.0 * fac,
                                            2989.0 * fac, 2989.0 * fac, 1323.0 * fac,
                                            3577.0 * fac, 751.0 * fac};
  return lagrange_interpolate(start, end, points, position, facs);
}

// Integrate using newton-cotes 7 rule
BoutReal integrate_lagrange_curve_nc7(const BoutReal start, const BoutReal end,
                                      const std::deque<BoutReal>& points,
                                      const int position) {
  AUTO_TRACE();
  constexpr std::size_t size = 7;
  constexpr BoutReal fac = 1.0 / 140.0;
  constexpr std::array<BoutReal, size> facs{41.0 * fac,  216.0 * fac, 27.0 * fac,
                                            272.0 * fac, 27.0 * fac,  216.0 * fac,
                                            41.0 * fac};
  return lagrange_interpolate(start, end, points, position, facs);
}

// Integrate using newton-cotes 6 rule
BoutReal integrate_lagrange_curve_nc6(const BoutReal start, const BoutReal end,
                                      const std::deque<BoutReal>& points,
                                      const int position) {
  AUTO_TRACE();
  constexpr std::size_t size = 6;
  constexpr BoutReal fac = 5.0 / 288.0;
  constexpr std::array<BoutReal, size> facs{19.0 * fac, 75.0 * fac, 50.0 * fac,
                                            50.0 * fac, 75.0 * fac, 19.0 * fac};
  return lagrange_interpolate(start, end, points, position, facs);
};

// Integrate using newton-cotes 5 rule (Boole)
BoutReal integrate_lagrange_curve_nc5(const BoutReal start, const BoutReal end,
                                      const std::deque<BoutReal>& points,
                                      const int position) {
  AUTO_TRACE();
  constexpr std::size_t size = 5;
  constexpr BoutReal fac = 2.0 / 45.0;
  constexpr std::array<BoutReal, size> facs{7.0 * fac, 32.0 * fac, 12.0 * fac, 32.0 * fac,
                                            7.0 * fac};
  return lagrange_interpolate(start, end, points, position, facs);
}

// Integrate using newton-cotes 4 rule (Simpson 3/8)
BoutReal integrate_lagrange_curve_nc4(const BoutReal start, const BoutReal end,
                                      const std::deque<BoutReal>& points,
                                      const int position) {
  AUTO_TRACE();
  constexpr std::size_t size = 4;
  constexpr BoutReal fac = 3.0 / 8.0;
  constexpr std::array<BoutReal, size> facs{1.0 * fac, 3.0 * fac, 3.0 * fac, 1.0 * fac};
  return lagrange_interpolate(start, end, points, position, facs);
}

// Integrate using newton-cotes 3 rule (Simpson)
BoutReal integrate_lagrange_curve_nc3(const BoutReal start, const BoutReal end,
                                      const std::deque<BoutReal>& points,
                                      const int position) {
  AUTO_TRACE();
  constexpr std::size_t size = 3;
  constexpr BoutReal fac = 1.0 / 3.0;
  constexpr std::array<BoutReal, size> facs{1.0 * fac, 4.0 * fac, 1.0 * fac};
  return lagrange_interpolate(start, end, points, position, facs);
}

// Integrate using newton-cotes 2 rule (Trap)
BoutReal integrate_lagrange_curve_nc2(const BoutReal start, const BoutReal end,
                                      const std::deque<BoutReal>& points,
                                      const int position) {
  AUTO_TRACE();
  constexpr std::size_t size = 2;
  constexpr BoutReal fac = 1.0 / 2.0;
  constexpr std::array<BoutReal, size> facs{1.0 * fac, 1.0 * fac};
  return lagrange_interpolate(start, end, points, position, facs);
}

// Integrate lagrange polynomial to find the coefficienst of the requested order
BoutReal integrate_lagrange_curve(const BoutReal start, const BoutReal end,
                                  const std::deque<BoutReal>& points, const int position,
                                  const int order) {
  AUTO_TRACE();

  switch (order) {
  case 1:
    return integrate_lagrange_curve_nc2(start, end, points, position);
  case 2:
    return integrate_lagrange_curve_nc3(start, end, points, position);
  case 3:
    return integrate_lagrange_curve_nc4(start, end, points, position);
  case 4:
    return integrate_lagrange_curve_nc5(start, end, points, position);
  case 5:
    return integrate_lagrange_curve_nc6(start, end, points, position);
  case 6:
    return integrate_lagrange_curve_nc7(start, end, points, position);
  case 7:
    return integrate_lagrange_curve_nc8(start, end, points, position);
  default:
    return integrate_lagrange_curve_nc9(start, end, points, position);
  }
}

// Calculate the set of Adams-Bashforth coefficients required to get from t = points[0]
// to t = nextPoint
// at the requested order.
std::vector<BoutReal> get_adams_bashforth_coefficients(const BoutReal nextPoint,
                                                       const std::deque<BoutReal>& points,
                                                       const int order) {
  AUTO_TRACE();
  ASSERT2(order <= points.size());

  std::vector<BoutReal> result;
  result.reserve(order);

  for (int i = 0; i < order; i++) {
    result.emplace_back(integrate_lagrange_curve(points[0], nextPoint, points, i, order));
  }

  return result;
}

// In-place Adams-Bashforth integration
void AB_integrate_update(Array<BoutReal>& update, BoutReal timestep,
                         const std::deque<BoutReal>& times,
                         const std::deque<Array<BoutReal>>& history, int order) {

  const auto AB_coefficients = get_adams_bashforth_coefficients(timestep, times, order);

  for (std::size_t j = 0; j < static_cast<std::size_t>(order); ++j) {
    const BoutReal factor = AB_coefficients[j];
    BOUT_OMP(parallel for);
    for (std::size_t i = 0; i < static_cast<std::size_t>(update.size()); ++i) {
      update[i] += history[j][i] * factor;
    }
  }
}

// Integrate \p history with Adams-Bashforth of order \p order
Array<BoutReal> AB_integrate(int nlocal, BoutReal timestep,
                             const std::deque<BoutReal>& times,
                             const std::deque<Array<BoutReal>>& history, int order) {
  Array<BoutReal> update(nlocal);

  // Zero-initialise to ensure we can operate on the contiguous
  // history arrays in order
  std::fill(std::begin(update), std::end(update), 0.0);

  AB_integrate_update(update, timestep, times, history, order);
  return update;
}

// Free function to return an estimate of the factor by which a
// timestep giving aerror = error should be scaled to give aerror =
// tolerance when using a scheme of order = order, where aerror =
// abs(soln_accurate - soln_approx)
BoutReal get_timestep_limit(const BoutReal error, const BoutReal tolerance,
                            const int order) {
  return std::exp(-std::log(error / tolerance) / order);
};

/// Finds the maximum absolute error, i.e. Max(Abs(stateApprox - stateAccurate))
/// over all processors.
BoutReal get_error(const Array<BoutReal>& stateApprox,
                   const Array<BoutReal>& stateAccurate) {
  AUTO_TRACE();
  BoutReal local_result = 0.0;
  BoutReal err = 0.0;

  const auto nlocal = stateAccurate.size();
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
  if (MPI_Allreduce(&local_result, &err, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get()) != 0) {
    throw BoutException("MPI_Allreduce failed");
  }
  return err;
}
} // namespace

AdamsBashforthSolver::AdamsBashforthSolver(Options* options) : Solver(options) {
  AUTO_TRACE();
  canReset = true;
}

void AdamsBashforthSolver::setMaxTimestep(BoutReal dt) {
  AUTO_TRACE();
  if (dt > timestep)
    return; // Already less than this

  if (adaptive) // Should we throw if we're not adaptive as we've tried to set a timestep
                // limit but couldn't?
    timestep = dt; // Won't be used this time, but next
}

int AdamsBashforthSolver::init(int nout, BoutReal tstep) {

  TRACE("Initialising AdamsBashforth solver");

  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;

  output << "\n\tAdams-Bashforth (explicit) multistep solver\n";

  nsteps = nout; // Save number of output steps
  out_timestep = tstep;
  max_dt = tstep;

  // Calculate number of variables
  nlocal = getLocalN();

  // Get total problem size
  int ntmp;
  if (MPI_Allreduce(&nlocal, &ntmp, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed!");
  }
  neq = ntmp;

  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n", n3Dvars(),
               n2Dvars(), neq, nlocal);

  // Get options
  atol = (*options)["atol"]
             .doc("Absolute tolerance")
             .withDefault(1.e-5); // Not used, just here for parity
  rtol = (*options)["rtol"].doc("Relative tolerance").withDefault(1.e-5);
  dtFac = (*options)["dtFac"]
              .doc("Factor by which we scale timestep estimate when adapating")
              .withDefault(0.75);
  max_timestep = (*options)["max_timestep"].doc("Maximum timestep").withDefault(tstep);
  timestep = (*options)["timestep"].doc("Starting timestep").withDefault(max_timestep);
  mxstep = (*options)["mxstep"]
               .doc("Maximum number of steps taken between outputs")
               .withDefault(50000);
  adaptive =
      (*options)["adaptive"].doc("Adapt internal timestep using rtol.").withDefault(true);
  adaptive_order = (*options)["adaptive_order"]
                       .doc("Adapt algorithm order using rtol.")
                       .withDefault(true);

  maximum_order =
      (*options)["order"].doc("The requested maximum order of the scheme").withDefault(5);
  followHighOrder =
      (*options)["followHighOrder"]
          .doc("If true and adaptive then use the more accurate solution as result.")
          .withDefault(true);

  // Check if the requested timestep in the non-adaptive case would lead to us
  // effectively violating the MXSTEP specified.
  if (not adaptive and (out_timestep / timestep > mxstep)) {
    throw BoutException("ERROR: Requested timestep would lead to MXSTEP being exceeded. "
                        "timestep = {:e}, MXSTEP={:d}\n",
                        timestep, mxstep);
  }

  // Put starting values into states
  state.reallocate(nlocal);
  nextState.reallocate(nlocal);
  std::fill(std::begin(nextState), std::end(nextState), 0.0);
  save_vars(std::begin(state));

  // Set the starting order
  current_order = 1;

  return 0;
}

void AdamsBashforthSolver::resetInternalFields() {
  AUTO_TRACE();

  // History and times
  history.clear();
  times.clear();

  // Order
  current_order = 1;

  // States
  std::fill(std::begin(nextState), std::end(nextState), 0.0);
  save_vars(std::begin(state));
}

int AdamsBashforthSolver::run() {
  AUTO_TRACE();

  // Just for developer diagnostics
  int nwasted = 0;
  int nwasted_following_fail = 0;

  for (int s = 0; s < nsteps; s++) {
    BoutReal target = simtime + out_timestep;

    bool running = true;
    int internal_steps = 0;

    // Take a single output time step
    while (running) {
      // Here's the derivative calculation at the current time
      // Find d state/dt and store in history -- this doesn't
      // need repeating whilst adapting timestep
      run_rhs(simtime);
      history.emplace_front(nlocal);
      save_derivs(std::begin(history[0]));
      times.emplace_front(simtime);

      // Just for developer diagnostics - set to true when the previous
      // attempt at a time step failed.
      bool previous_fail = false;

      // Flag to indicate if we want to use a lower order method
      bool use_lower = false;

      BoutReal dt;

      // Take a single internal time step
      while (true) {
        // Limit the timestep to the specified maximum
        timestep = std::min(timestep, max_timestep);

        // We actually use dt to reflect the timestep actually used for an advance
        // as we may modify it.
        dt = timestep;

        // If running is true then we haven't yet finished this output step
        running = true;

        // Check if we're going to reach our target time and adjust
        // the timestep to ensure we don't go past it. Note this means
        // that even when non-adaptive we may end up changing the
        // timestep occassionally. This is ok here as the timestep
        // code is completely general, but could potentially be an
        // issue for other solvers.
        if ((simtime + dt) >= target) {
          dt = target - simtime;
          running = false;
        }

        // Take a step and get the error if adaptive
        const BoutReal err = take_step(simtime, dt, current_order, state, nextState);

        // Calculate and check error if adaptive
        if (not adaptive) {
          break;
        }

        // Really the following should apply to both adaptive and non-adaptive
        // approaches, but the non-adaptive can be determined without needing
        // to do any solves so we check during init instead.
        ++internal_steps;
        if (internal_steps > mxstep) {
          throw BoutException("ERROR: MXSTEP exceeded. timestep = {:e}, err={:e}\n",
                              timestep, err);
        }

        // Estimate the limiting timestep and update.  This is
        // really the estimate of the timestep for this current step
        // that would just satisfy the tolerance. In cases where we
        // move to the next step we actually end up using this new
        // timestep for the next step.
        BoutReal dt_lim = dt * get_timestep_limit(err, rtol, current_order);

        if (err < rtol) { // Successful step

          // Now we can consider what result we would get at
          // lower/higher order Our timestep limit gets smaller as
          // the order increases for fixed error, hence we really
          // want to use the lowest order that satisfies the
          // tolerance. Or in other words we want to use the order
          // that gives us the biggest timestep. For now we just see
          // what the error is when using one order lower.
          //
          // For now we only do this when we've had a successful
          // step, in general we might want to do this for failing
          // steps as well, but as the error drops quicker with
          // higher orders we might hope higher order is better when
          // the error condition is not met.
          if (adaptive_order and current_order > 1) {
            Array<BoutReal> lowerNextState(nlocal);
            // Currently we just reuse the existing code to take a
            // step but just do it with lower order
            // coefficients. This means we have to do another rhs
            // call. We might be able to get away with reusing the
            // half point derivatives from the higher order method
            // here instead, which would save the rhs call. This may
            // mean we don't trust the error as much and hence have
            // to scale the timestep more conservatively but this
            // may be worth it.
            //
            // Actually currently we do skip the second rhs call
            // and instead try to reuse the existing data.
            const BoutReal lowerErr =
                take_step(simtime, dt, current_order - 1, state, lowerNextState);

            const BoutReal lower_dt_lim =
                dt * get_timestep_limit(lowerErr, rtol, current_order - 1);

            // Decide if we want to use the lower order method based
            // on which gives us the biggest timestep.
            use_lower = lower_dt_lim > dt_lim;

            // If we decide the lower order is better then swap/set
            // the associated values to use the lower order result.
            if (use_lower) {
              dt_lim = lower_dt_lim;
              swap(nextState, lowerNextState);
              current_order = current_order - 1;
            }
          }

          // Try to limit increases in the timestep to no more than 10%.
          // We could/should make these numbers runtime to give more
          // control to the users, just wary of option overload.
          timestep = std::min(timestep * 1.1, dt_lim);

          // For developers
          previous_fail = false;

          break;
        }
        // Be more conservative if we've failed;
        timestep = 0.9 * dt_lim;

        // For developers
        if (previous_fail) {
          nwasted_following_fail++;
        }
        previous_fail = true;
        nwasted++;
      }

      // Ditch last history point if we have enough
      if (times.size() == maximum_order)
        times.pop_back();
      if (history.size() == maximum_order)
        history.pop_back();

      if (current_order < maximum_order) {
        // Don't increase the order if we wanted to use the lower order.
        if (not use_lower)
          current_order++;
      }

      // Taken an internal step, update times
      simtime += dt;

      // Put the new state into state.
      swap(state, nextState);

      // Put the state into the fields
      load_vars(std::begin(state));

      // Call the per internal timestep monitors
      call_timestep_monitors(simtime, dt);
    };

    // Put result into variables
    load_vars(std::begin(state));

    // Ensure aux. variables are up to date. In the future it would be nice to
    // provide a calc_aux(simtime) method on PhysicsModel (that could default to
    // calling rhs) which ensures the aux. variables are up to date in order to
    // avoid any additional unrequired work associated with run_rhs.
    run_rhs(simtime);

    // Advance iteration number
    iteration++;

    // Call the output step monitor function
    if (call_monitors(simtime, s, nsteps))
      break; // Stop simulation
  }

#if CHECK > 4
  output.write("\nNumber of wasted steps = {} and following a fail {}\n\n", nwasted,
               nwasted_following_fail);
#endif
  return 0;
}

// Updates the internal state (?) along with an error estimate?
// Should probably just try taking a step of given size with given
// order, leaving the error calculation for calling code
BoutReal AdamsBashforthSolver::take_step(const BoutReal timeIn, const BoutReal dt,
                                         const int order, Array<BoutReal>& current,
                                         Array<BoutReal>& result) {
  AUTO_TRACE();

  // The initial error is 0.0
  BoutReal err = 0.0;

  Array<BoutReal> full_update = AB_integrate(nlocal, timeIn + dt, times, history, order);

  // Calculate the new state given the history and current state.
  // Could possibly skip the following calculation if adaptive and following the high
  // order method.
  // Possible to write this using algorithms, but until c++ 17 probably prefer the
  // explicit loop as
  // clearer, compatible with OMP and empirically slightly faster.
  // std::transform(std::begin(current), std::end(current), std::begin(full_update),
  //                std::begin(result), std::plus<BoutReal>{});
  if (not(adaptive and followHighOrder)) {
    BOUT_OMP(parallel for);
    for (int i = 0; i < nlocal; i++) {
      result[i] = current[i] + full_update[i];
    }
  }

  if (adaptive) {

    // Use this variable to say how big the first small timestep should be as a fraction
    // of the large timestep, dt. Here fixed to 0.5 to take two equally sized half steps
    // but left here to enable developer experimentation.
    constexpr BoutReal firstPart = 0.5;

    // Take a small time step - note we don't need to call the rhs again just yet
    Array<BoutReal> half_update =
        AB_integrate(nlocal, timeIn + (dt * firstPart), times, history, order);

    // -------------------------------------------
    // Now do the second small timestep -- note we need to call rhs again
    // -------------------------------------------

    // Add storage to history and the current time to times
    history.emplace_front(nlocal);
    times.emplace_front(timeIn + dt * firstPart);

    // Put intermediate result into variables, call rhs and save the derivatives
    // Try to cheat for now with this HACK. If the order /=
    // current_order then call must be part of the adapative_order code
    // so don't recalculate just reuse stored derivatives.
    if (order == current_order) {
      Array<BoutReal> result2(nlocal);

      // Now we have to calculate the state after the first small step as we will need to
      // use this to calculate the derivatives at this point.
      // std::transform(std::begin(current), std::end(current), std::begin(half_update),
      //                std::begin(result2), std::plus<BoutReal>{});
      BOUT_OMP(parallel for);
      for (int i = 0; i < nlocal; i++) {
        result2[i] = current[i] + half_update[i];
      };

      load_vars(std::begin(result2));
      // This is typically the most expensive part of this routine.
      //
      run_rhs(timeIn + firstPart * dt);

      // Restore fields to the original state
      load_vars(std::begin(current));
    }
    save_derivs(std::begin(history[0]));

    // Finish the time step
    AB_integrate_update(half_update, timeIn + dt, times, history, order);

    // Drop the temporary history information
    history.pop_front();
    times.pop_front();

    // Here we calculate the error by comparing the updates rather than output states
    // this is to avoid issues where we have large fields but small derivatives (i.e. to
    // avoid possible numerical issues at looking at the difference between two large
    // numbers).
    err = get_error(full_update, half_update);

    // Note here we don't add a small change onto result, we recalculate using the
    // "full" two half step half_update. Rather than using result2 we just replace
    // result here as we want to use this smaller step result
    if (followHighOrder) {
      BOUT_OMP(parallel for);
      for (int i = 0; i < nlocal; i++) {
        result[i] = current[i] + half_update[i];
      };
    }
  }

  return err;
}

