#include "split-rk.hxx"

int SplitRK::init(int nout, BoutReal tstep) {
  AUTO_TRACE();

  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;

  output.write(_("\n\tSplit Runge-Kutta-Legendre and SSP-RK3 solver\n"));

  nsteps = nout; // Save number of output steps
  out_timestep = tstep;

  // Calculate number of variables
  nlocal = getLocalN();
  
  // Get total problem size
  if(MPI_Allreduce(&nlocal, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed!");
  }
  
  // Allocate memory
  state.reallocate(nlocal);

  // memory for taking a single time step
  u1.reallocate(nlocal);
  u2.reallocate(nlocal);
  u3.reallocate(nlocal);
  dydt.reallocate(nlocal);
  
  // Put starting values into f
  save_vars(std::begin(state));

  // Get options. Default values for many of these are set in constructor.
  auto &opt = *options;
  timestep = opt["timestep"]
                 .doc("Internal timestep. This may be rounded down.")
                 .withDefault(out_timestep);

  adaptive = opt["adaptive"].doc("Use accuracy tolerances to adapt timestep?").withDefault(adaptive);

  atol = opt["atol"].doc("Absolute tolerance").withDefault(atol);
  rtol = opt["rtol"].doc("Relative tolerance").withDefault(rtol);

  max_timestep = opt["max_timestep"].doc("Maximum timestep. Negative means no limit.").withDefault(out_timestep);

  max_timestep_change = opt["max_timestep_change"]
                            .doc("Maximum factor by which the timestep should be changed. Must be >1")
                            .withDefault(max_timestep_change);
  ASSERT0(max_timestep_change > 1.0);
  
  mxstep = opt["mxstep"]
               .doc("Maximum number of internal steps between outputs")
               .withDefault(mxstep);
  ASSERT0(mxstep > 0);

  adapt_period = opt["adapt_period"]
          .doc("Number of steps between tolerance checks. 1 means check every step.")
          .withDefault(adapt_period);

  if (adaptive) {
    // Need additional storage to compare results after time steps
    state1.reallocate(nlocal);
    state2.reallocate(nlocal);
  }
  
  ASSERT0(adapt_period > 0);
  
  int ninternal_steps = static_cast<int>(std::ceil(out_timestep / timestep));
  ASSERT0(ninternal_steps > 0);
  
  timestep = out_timestep / ninternal_steps;
  output.write(_("\tUsing a timestep %e\n"), timestep);
  
  nstages = opt["nstages"].doc("Number of stages in RKL step. Must be > 1").withDefault(10);
  ASSERT0(nstages > 1);

  diagnose = opt["diagnose"].doc("Print diagnostic information?").withDefault(diagnose);

  return 0;
}

int SplitRK::run() {
  AUTO_TRACE();

  for (int step = 0; step < nsteps; step++) {
    // Take an output step

    BoutReal target = simtime + out_timestep;

    BoutReal dt;  // The next timestep to take
    bool running = true;  // Changed to false to break out of inner loop
    int internal_steps = 0;  // Quit if this exceeds mxstep
    
    do {
      // Take a single time step

      if (adaptive and (internal_steps % adapt_period == 0)) {
        do {
          // Keep adapting the timestep until the error is within tolerances
          
          dt = timestep;
          running = true; // Reset after maybe adapting timestep
          if ((simtime + dt) >= target) {
            dt = target - simtime; // Make sure the last timestep is on the output 
            running = false; // Fall out of this inner loop after this step
          }
          
          // Take two half-steps
          take_step(simtime,          0.5*dt, state, state1);
          take_step(simtime + 0.5*dt, 0.5*dt, state1, state2);
          
          // Take a full step
          take_step(simtime, dt, state, state1);

          // Check accuracy
          BoutReal local_err = 0.;
          BOUT_OMP(parallel for reduction(+: local_err)   )
          for (int i = 0; i < nlocal; i++) {
            local_err += fabs(state2[i] - state1[i]) / (fabs(state1[i]) + fabs(state2[i]) + atol);
          }
          
          // Average over all processors
          BoutReal err;
          if (MPI_Allreduce(&local_err, &err, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get())) {
            throw BoutException("MPI_Allreduce failed");
          }

          err /= static_cast<BoutReal>(neq);
          
          internal_steps++;
          if (internal_steps > mxstep) {
            throw BoutException("ERROR: MXSTEP exceeded. timestep = {:e}, err={:e}\n",
                                timestep, err);
          }

          if (diagnose) {
            output.write("\nError: %e. atol=%e, rtol=%e\n", err, atol, rtol);
          }

          if ((err > rtol) || (err < 0.1 * rtol)) {
            // Need to change timestep. Error ~ dt^2

            BoutReal factor = pow((0.5 * rtol) / err, 1./3);

            if (factor > max_timestep_change) {
              factor = max_timestep_change;
            } else if (factor < 1. / max_timestep_change) {
              factor = 1. / max_timestep_change;
            }

            timestep *= factor;

            if ((max_timestep > 0) && (timestep > max_timestep)) {
              timestep = max_timestep;
            }

            if (diagnose) {
              output.write("\tAdapting. timestep %e (factor %e). Max=%e\n", timestep, factor, max_timestep);
            }
          }
          if (err < rtol) {
            swap(state, state2); // Put result in state
            break; // Acceptable accuracy
          }
        } while (true);
      } else {
        // No adaptive timestepping. Take a single step

        dt = timestep;
        running = true; // Reset after maybe adapting timestep
        if ((simtime + dt) >= target) {
          dt = target - simtime; // Make sure the last timestep is on the output 
          running = false; // Fall out of this inner loop after this step
        }
        
        take_step(simtime, timestep, state, state);
        internal_steps++;
      }
      
      simtime += dt;
      call_timestep_monitors(simtime, timestep);
      
    } while (running);
      
    load_vars(std::begin(state)); // Put result into variables
    // Call rhs function to get extra variables at this time
    run_rhs(simtime);
    
    iteration++; // Advance iteration number
    
    /// Call the monitor function
    
    if(call_monitors(simtime, step, nsteps)) {
      // User signalled to quit
      break;
    }
  }
  return 0;
}

void SplitRK::take_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start,
                        Array<BoutReal>& result) {
  // Half step
  take_diffusion_step(curtime, 0.5*dt, start, result);
  
  // Full step
  take_advection_step(curtime, dt, result, result);

  // Half step
  take_diffusion_step(curtime + 0.5*dt, 0.5*dt, result, result);
}

void SplitRK::take_diffusion_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start,
                                  Array<BoutReal>& result) {

  const BoutReal weight = dt * 4./(SQ(nstages) + nstages - 2);
  
  load_vars(std::begin(start));
  run_diffusive(curtime);
  save_derivs(std::begin(dydt));   // dydt = f(y0)

  // Stage j = 1
  // y_m2 = y0 + weight/3.0 * f(y0)  -> u2

  BOUT_OMP(parallel for)
  for (int i = 0; i < dydt.size(); i++) {
    u2[i] = start[i] + (weight/3.0) * dydt[i];
  }
  
  // Stage j = 2
  // mu = 1.5, nu terms cancel
  load_vars(std::begin(u2));
  run_diffusive(curtime + (weight/3.0) * dt);
  save_derivs(std::begin(u3)); // f(y_m2) -> u3
  
  BOUT_OMP(parallel for)
  for (int i = 0; i < u3.size(); i++) {
    u1[i] = 1.5 * (u2[i] + weight * u3[i]) - 0.5 * start[i] - weight * dydt[i];
  }
  
  BoutReal b_jm2 = 1. / 3; // b_{j - 2}
  BoutReal b_jm1 = 1. / 3; // b_{j - 1}
  
  for (int j = 3; j <= nstages; j++) {
    
    BoutReal b_j = (SQ(j) + j - 2.0) / (2.*j * (j + 1.));
    
    BoutReal mu = (2.*j - 1.)/j * b_j / b_jm1;
    BoutReal nu = -(j - 1.)/j * b_j / b_jm2;
    BoutReal a_jm1 = 1. - b_jm1;

    load_vars(std::begin(u1));
    run_diffusive(curtime);
    save_derivs(std::begin(u3)); // f(y_m1) -> u3
    
    BOUT_OMP(parallel for)
    for (int i = 0; i < u3.size(); i++) {
      // Next stage result in u3
      u3[i] = mu * (u1[i] + weight * (u3[i] - a_jm1 * dydt[i])) + nu * u2[i]
              + (1. - mu - nu) * start[i];
    }

    // Cycle values
    b_jm2 = b_jm1;
    b_jm1 = b_j;

    // Cycle u2 <- u1 <- u3 <- u2
    // so that no new memory is allocated, and no arrays point to the same data
    swap(u1, u2);
    swap(u1, u3); 

    // Most recent now in u1, then u2, then u3
  }
  swap(u1, result);
}

void SplitRK::take_advection_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start,
                                  Array<BoutReal>& result) {
  const int nlocal = getLocalN();
  
  load_vars(std::begin(start));
  run_convective(curtime);
  save_derivs(std::begin(dydt));

  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    u1[i] = start[i] + dt*dydt[i];

  load_vars(std::begin(u1));
  run_convective(curtime + dt);
  save_derivs(std::begin(dydt));

  BOUT_OMP(parallel for )
  for(int i=0;i<nlocal;i++)
    u2[i] = 0.75*start[i] + 0.25*u1[i] + 0.25*dt*dydt[i];

  load_vars(std::begin(u2));
  run_convective(curtime + 0.5*dt);
  save_derivs(std::begin(dydt));

  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    result[i] = (1./3)*start[i] + (2./3.)*(u2[i] + dt*dydt[i]);
}
