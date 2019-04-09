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
  const int nlocal = getLocalN();
  
  // Allocate memory
  state.reallocate(nlocal);

  // memory for taking a single time step
  u1.reallocate(nlocal);
  u2.reallocate(nlocal);
  u3.reallocate(nlocal);
  L.reallocate(nlocal);
  
  // Put starting values into f
  save_vars(std::begin(state));

  // Get options
  auto &opt = *options;
  timestep = opt["timestep"]
                 .doc("Internal timestep. This may be rounded down.")
                 .withDefault(out_timestep);

  ninternal_steps = static_cast<int>(std::ceil(out_timestep / timestep));

  ASSERT0(ninternal_steps > 0);

  timestep = out_timestep / ninternal_steps;
  output.write(_("\tUsing a timestep %e"), timestep);
  
  nstages = opt["nstages"].doc("Number of stages in RKL step. Must be > 1").withDefault(10);
  ASSERT0(nstages > 1);
  
  return 0;
}

int SplitRK::run() {
  AUTO_TRACE();

  for (int step = 0; step < nsteps; step++) {
    // Take an output step

    for (int internal_step = 0; internal_step < ninternal_steps; internal_step++) {
      // Take a single step
      state = take_step(simtime, timestep, state);
      
      simtime += timestep;
      
      call_timestep_monitors(simtime, timestep);
    }

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

Array<BoutReal> SplitRK::take_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start) {
  return take_diffusion_step(curtime + 0.5*dt, 0.5*dt,  // Half step
                             take_advection_step(curtime, dt,  // Full step
                                                 take_diffusion_step(curtime, 0.5*dt, start))); // Half step
}

Array<BoutReal> SplitRK::take_diffusion_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start) {

  const BoutReal weight = dt * 4./(SQ(nstages) + nstages - 2);
  
  load_vars(std::begin(start));
  run_diffusive(curtime);
  save_derivs(std::begin(L));   // L = f(y0)

  // Stage j = 1
  // y_m2 = y0 + weight/3.0 * f(y0)  -> u2

  BOUT_OMP(parallel for)
  for (int i = 0; i < L.size(); i++) {
    u2[i] = start[i] + (weight/3.0) * L[i];
  }
  
  // Stage j = 2
  // mu = 0.75, nu terms cancel
  load_vars(std::begin(u2));
  run_diffusive(curtime);
  save_derivs(std::begin(u3)); // f(y_m2) -> u3
  
  BOUT_OMP(parallel for)
  for (int i = 0; i < u3.size(); i++) {
    u1[i] = 0.75 * u2[i] + 0.25 * start[i] + 0.75 * weight * u3[i] - 0.5 * weight * L[i];
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
      u3[i] = mu * u1[i] + nu * u2[i] + (1. - mu - nu) * start[i]
        + mu * weight * u3[i] - a_jm1 * mu * weight * L[i];
    }

    // Cycle values
    b_jm2 = b_jm1;
    b_jm1 = b_j;

    // Cycle u2 <- u1 <- u3 <- u2
    // so that no new memory is allocated, and no arrays point to the same data
    Array<BoutReal> tmp = u2;
    u2 = u1;
    u1 = u3;
    u3 = tmp;

    // Most recent now in u1, then u2, then u3
  }
  return u1;
}

Array<BoutReal> SplitRK::take_advection_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start) {
  const int nlocal = getLocalN();
  
  load_vars(std::begin(start));
  run_convective(curtime);
  save_derivs(std::begin(L));

  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    u1[i] = start[i] + dt*L[i];

  load_vars(std::begin(u1));
  run_convective(curtime + dt);
  save_derivs(std::begin(L));

  BOUT_OMP(parallel for )
  for(int i=0;i<nlocal;i++)
    u2[i] = 0.75*start[i] + 0.25*u1[i] + 0.25*dt*L[i];

  load_vars(std::begin(u2));
  run_convective(curtime + 0.5*dt);
  save_derivs(std::begin(L));

  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    result[i] = (1./3)*start[i] + (2./3.)*(u2[i] + dt*L[i]);
}
