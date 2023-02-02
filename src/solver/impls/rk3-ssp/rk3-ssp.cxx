
#include "rk3-ssp.hxx"

#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/msg_stack.hxx>
#include <bout/openmpwrap.hxx>
#include <bout/utils.hxx>
#include <cmath>

#include <bout/output.hxx>

RK3SSP::RK3SSP(Options* opt)
    : Solver(opt), max_timestep((*options)["max_timestep"]
                                    .doc("Maximum timestep")
                                    .withDefault(getOutputTimestep())),
      timestep((*options)["timestep"].doc("Starting timestep").withDefault(max_timestep)),
      mxstep((*options)["mxstep"]
                 .doc("Maximum number of steps between outputs")
                 .withDefault(500)) {}

void RK3SSP::setMaxTimestep(BoutReal dt) {
  if (dt > timestep) {
    return; // Already less than this
  }

  timestep = dt; // Won't be used this time, but next
}

int RK3SSP::init() {
  TRACE("Initialising RK3 SSP solver");

  Solver::init();
  output << "\n\tRunge-Kutta 3rd-order SSP solver\n";

  // Calculate number of variables
  nlocal = getLocalN();

  // Get total problem size
  int ntmp;
  if (bout::globals::mpi->MPI_Allreduce(&nlocal, &ntmp, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed!");
  }
  neq = ntmp;

  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n", n3Dvars(),
               n2Dvars(), neq, nlocal);

  // Allocate memory
  f.reallocate(nlocal);

  // memory for taking a single time step
  u1.reallocate(nlocal);
  u2.reallocate(nlocal);
  u3.reallocate(nlocal);
  L.reallocate(nlocal);

  // Put starting values into f
  save_vars(std::begin(f));

  return 0;
}

int RK3SSP::run() {
  TRACE("RK3SSP::run()");

  for (int s = 0; s < getNumberOutputSteps(); s++) {
    BoutReal target = simtime + getOutputTimestep();

    BoutReal dt;
    bool running = true;
    do {
      // Take a single time step

      dt = timestep;
      running = true;
      if ((simtime + dt) >= target) {
        dt = target - simtime; // Make sure the last timestep is on the output
        running = false;
      }
      output.write("t = {:e}, dt = {:e}\n", simtime, dt);
      // No adaptive timestepping for now
      take_step(simtime, dt, f, f);

      simtime += dt;

      call_timestep_monitors(simtime, dt);
    } while (running);

    load_vars(std::begin(f)); // Put result into variables
    // Call rhs function to get extra variables at this time
    run_rhs(simtime);

    if (call_monitors(simtime, s, getNumberOutputSteps())) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

void RK3SSP::take_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start,
                       Array<BoutReal>& result) {

  load_vars(std::begin(start));
  run_rhs(curtime);
  save_derivs(std::begin(L));

  BOUT_OMP(parallel for)
  for (int i = 0; i < nlocal; i++) {
    u1[i] = start[i] + dt * L[i];
  }

  load_vars(std::begin(u1));
  run_rhs(curtime + dt);
  save_derivs(std::begin(L));

  BOUT_OMP(parallel for )
  for (int i = 0; i < nlocal; i++) {
    u2[i] = 0.75 * start[i] + 0.25 * u1[i] + 0.25 * dt * L[i];
  }

  load_vars(std::begin(u2));
  run_rhs(curtime + 0.5 * dt);
  save_derivs(std::begin(L));

  BOUT_OMP(parallel for)
  for (int i = 0; i < nlocal; i++) {
    result[i] = (1. / 3) * start[i] + (2. / 3.) * (u2[i] + dt * L[i]);
  }
}
