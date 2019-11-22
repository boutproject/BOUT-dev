#include "adams_bashforth.hxx"

#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <utils.hxx>

#include <output.hxx>

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

  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n", n3Dvars(),
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
                        "timestep = %e, MXSTEP=%i\n",
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
        if (adaptive) {
          // Really the following should apply to both adaptive and non-adaptive
          // approaches, but the non-adaptive can be determined without needing
          // to do any solves so we check during init instead.
          internal_steps++;
          if (internal_steps > mxstep)
            throw BoutException("ERROR: MXSTEP exceeded. timestep = %e, err=%e\n",
                                timestep, err);

          // Estimate the limiting timestep and update.  This is
          // really the estimate of the timestep for this current step
          // that would just satisfy the tolerance. In cases where we
          // move to the next step we actually end up using this new
          // timestep for the next step.
          const BoutReal dt_lim = dt * exp(-log(err / rtol) / current_order);

          if (err < rtol) {
            // Try to limit increases in the timestep to no more than 5%.
            timestep = std::min(timestep * 1.05, dt_lim * 0.75);
            previous_fail = false;
            break;
          } else {
            // Be more conservative if we've failed;
            timestep = 0.75 * dt_lim;
            if (previous_fail) {
              nwasted_following_fail++;
            }
            previous_fail = true;
            nwasted++;
          }
        }
      }

      // Taken an internal step, update times
      simtime += dt;

      if (current_order == maximum_order) {
        // Ditch last history point
        times.pop_back();
        history.pop_back();
      } else {
        // Here we unconditionally increase the order if we've yet to
        // reach the maximum. In general it is probably better to
        // consider an adaptive order scheme to try to maximise the
        // timestep we can take. This is something to explore in the
        // future.
        current_order++;
      }

      // Call the per internal timestep monitors
      call_timestep_monitors(simtime, dt);

      // Put the new state into state.
      swap(state, nextState);
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
  output << "\nNumber of wasted steps = " << nwasted << " and following a fail "
         << nwasted_following_fail << "\n\n"
         << endl;
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

  // Calculate the coefficients for a single step of size dt
  const auto coefs =
      coefficients_calculator.get_adams_bashforth_coefficients(timeIn + dt, times, order);

  // Create some storage for the update to the state (i.e. state(timeIn + dt) = current +
  // full_update).
  Array<BoutReal> full_update(nlocal);

  // Note we split the work here into initialisation with std::fill
  // and a separate double loop to calculate the update. This is
  // to ensure we can operate on the contiguous arrays in history
  // in order.
  std::fill(std::begin(full_update), std::end(full_update), 0.0);

  for (int j = 0; j < order; j++) {
    const BoutReal factor = coefs[j];

    BOUT_OMP(parallel for);
    for (int i = 0; i < nlocal; i++) {
      full_update[i] += history[j][i] * factor;
    }
  }

  // Calculate the new state given the history and current state.
  // Could possibly skip the following calculation if adaptive and following the high
  // order method.
  // Possible to write this using algorithms, but until c++ 17 probably prefer the
  // explicit loop as
  // clearer, compatible with OMP and empirically slightly faster.
  // std::transform(std::begin(current), std::end(current), std::begin(full_update),
  //                std::begin(result), std::plus<BoutReal>{});
  BOUT_OMP(parallel for);
  for (int i = 0; i < nlocal; i++) {
    result[i] = current[i] + full_update[i];
  };

  if (adaptive) {

    // Create some storage for the small step update and corresponding resulting state
    Array<BoutReal> result2(nlocal);
    Array<BoutReal> half_update(nlocal);

    // Use this variable to say how big the first small timestep should be as a fraction
    // of the large timestep, dt. Here fixed to 0.5 to take two equally sized half steps
    // but left here to enable developer experimentation.
    constexpr BoutReal firstPart = 0.5;

    // -------------------------------------------
    // Take a small time step - note we don't need to call the rhs again just yet
    // -------------------------------------------

    // Calculate the coefficients to get to timeIn + dt * firstPart
    const auto coefsFirstStep = coefficients_calculator.get_adams_bashforth_coefficients(
        timeIn + dt * firstPart, times, order);

    // Initialise the update array to 0.
    std::fill(std::begin(half_update), std::end(half_update), 0.0);

    for (int j = 0; j < order; j++) {
      const BoutReal factor = coefsFirstStep[j];

      BOUT_OMP(parallel for);
      for (int i = 0; i < nlocal; i++) {
        half_update[i] += history[j][i] * factor;
      }
    }

    // Now we have to calculate the state after the first small step as we will need to
    // use this to calculate the derivatives at this point.
    // std::transform(std::begin(current), std::end(current), std::begin(half_update),
    //                std::begin(result2), std::plus<BoutReal>{});
    BOUT_OMP(parallel for);
    for (int i = 0; i < nlocal; i++) {
      result2[i] = current[i] + half_update[i];
    };

    // -------------------------------------------
    // Now do the second small timestep -- note we need to call rhs again
    // -------------------------------------------

    // Add storage to history and the current time to times
    history.emplace_front(nlocal);
    times.emplace_front(timeIn + dt * firstPart);

    // Put intermediate result into variables, call rhs and save the derivatives
    load_vars(std::begin(result2));
    run_rhs(timeIn
            + firstPart
                  * dt); // This is typically the most expensive part of this routine
    save_derivs(std::begin(history[0]));

    // Calculate the coefficients to get to timeIn + dt
    const auto coefsSecondStep = coefficients_calculator.get_adams_bashforth_coefficients(
        timeIn + dt, times, order);

    for (int j = 0; j < order; j++) {
      const BoutReal factor = coefsSecondStep[j];

      BOUT_OMP(parallel for);
      for (int i = 0; i < nlocal; i++) {
        half_update[i] += history[j][i] * factor;
      }
    }

    // Restore fields to the original state
    load_vars(std::begin(current));

    // Drop the temporary history information
    history.pop_front();
    times.pop_front();

    // Note here we don't add a small change onto result2, we recalculate using the
    // "full" two half step half_update
    // std::transform(std::begin(current), std::end(current), std::begin(half_update),
    //                std::begin(result2), std::plus<BoutReal>{});
    BOUT_OMP(parallel for);
    for (int i = 0; i < nlocal; i++) {
      result2[i] = current[i] + half_update[i];
    };

    // Here we calculate the error by comparing the updates rather than output states
    // this is to avoid issues where we have large fields but small derivatives (i.e. to
    // avoid possible numerical issues at looking at the difference between two large
    // numbers).
    err = get_error(full_update, half_update);

    // Swap the result to use the more accurate value if requested
    if (followHighOrder) {
      swap(result, result2);
    }
  }

  return err;
}
