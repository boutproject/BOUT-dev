/**************************************************************************
 * Karniadakis split-operator solver
 * 
 * Formulation from:
 * "GEM - An Energy Conserving Electromagnetic Gyrofluid Model"
 *  by Bruce D Scott. arXiv:physics/0501124v1 23 Jan 2005 
 *
 * Original paper:
 *   J. Comput. Phys. 97 (1991) p414-443
 * 
 * Always available, since doesn't depend on external library
 * 
 * Solves a system df/dt = S(f) + D(f)
 * 
 * where S is the RHS of each equation, and D is the diffusion terms
 * 
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
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

#include "karniadakis.hxx"

#include <boutcomm.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <bout/openmpwrap.hxx>

KarniadakisSolver::KarniadakisSolver(Options *options) : Solver(options) {
  canReset = true;  
}

int KarniadakisSolver::init(int nout, BoutReal tstep) {
  TRACE("Initialising Karniadakis solver");
  
  output_error << "\nWARNING:\n"
    "        The Karniadakis solver is now deprecated and will be removed in BOUT++ 5.0!\n"
    "        Try the \"splitrk\", \"imexbdf2\" (requires PETSc) or \"arkode\" (requires SUNDIALS)\n"
    "        solvers for other split-schemes\n\n";

  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;
  
  output << "\n\tKarniadakis solver\n";
 
  nsteps = nout; // Save number of output steps
  out_timestep = tstep;
  
  // Calculate number of variables
  nlocal = getLocalN();
  
  // Get total problem size
  int neq;
  if(MPI_Allreduce(&nlocal, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    output_error.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3Dvars(), n2Dvars(), neq, nlocal);
  
  // Allocate memory

  f1.reallocate(nlocal);
  f0.reallocate(nlocal);
  fm1.reallocate(nlocal);
  fm2.reallocate(nlocal);

  S0.reallocate(nlocal);
  Sm1.reallocate(nlocal);
  Sm2.reallocate(nlocal);

  D0.reallocate(nlocal);

  first_time = true;

  // Put starting values into f0
  save_vars(std::begin(f0));

  // Get options
  OPTION(options, timestep, tstep);
  
  // Make sure timestep divides into tstep
  
  // Number of sub-steps, rounded up
  nsubsteps = static_cast<int>(std::round(tstep / timestep));

  output.write("\tNumber of substeps: %e / %e -> %d\n", tstep, timestep, nsubsteps);

  timestep = tstep / static_cast<BoutReal>(nsubsteps);

  return 0;
}

int KarniadakisSolver::run() {
  TRACE("KarniadakisSolver::run()");
  
  for(int i=0;i<nsteps;i++) {
    // Run through a fixed number of steps
    for(int j=0; j<nsubsteps; j++) {
      // Advance f0 -> f1
      take_step(timestep);
      
      // Cycle buffers
      auto tmp = fm2;
      fm2 = fm1;
      fm1 = f0;
      f0 = f1;
      f1 = tmp;
      
      tmp = Sm2;
      Sm2 = Sm1;
      Sm1 = S0;
      S0 = tmp;
      
      simtime += timestep;
      
      call_timestep_monitors(simtime, timestep);
    }
    iteration++;
    
    // Call RHS to communicate and get auxilliary variables
    load_vars(std::begin(f0));
    run_rhs(simtime);
    
    /// Call the monitor function
    
    if(call_monitors(simtime, i, nsteps)) {
      // User signalled to quit
      break;
    }
  }
  
  return 0;
}

void KarniadakisSolver::resetInternalFields(){
  //Make sure all other fields get reset
  first_time=true;

  //Copy fields into vector
  save_vars(std::begin(f0));
}

void KarniadakisSolver::take_step(BoutReal dt) {
  // S0 = S(f0)

  load_vars(std::begin(f0));
  run_convective(simtime);
  save_derivs(std::begin(S0));

  if(first_time) {
    // Initialise values
    BOUT_OMP(parallel for)
    for(int i=0;i<nlocal;i++) {
    //fm1[i] = fm2[i] = f0[i];
      fm1[i] = f0[i] - dt*S0[i];
      fm2[i] = fm1[i] - dt*S0[i];
      Sm1[i] = Sm2[i] = S0[i];
    }
    first_time = false;
  }

  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    f1[i] = (6./11.) * (3.*f0[i] - 1.5*fm1[i] + (1./3.)*fm2[i] + dt*(3.*S0[i] - 3.*Sm1[i] + Sm2[i]));
  
  // D0 = S(f0)
  load_vars(std::begin(f0));
  run_diffusive(simtime);
  save_derivs(std::begin(D0));

  // f1 = f1 + dt*D0
  BOUT_OMP(parallel for)
  for(int i=0;i<nlocal;i++)
    f1[i] += (6./11.) * dt*D0[i];
}
