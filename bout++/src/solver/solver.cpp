/**************************************************************************
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

#include "globals.h"
#include "solver.h"
#include <string.h>

#include "initialprofiles.h"
#include "interpolation.h"
#include "solverfactory.h"

/**************************************************************************
 * Constructor
 **************************************************************************/

Solver::Solver() {
  // Set flags to defaults
  has_constraints = false;
  initialised = false;

  // Zero timing
  rhs_wtime = 0.0;
  rhs_ncalls = 0;

  // Restart directory
  restartdir = string("data");
}

/**************************************************************************
 * Add fields
 **************************************************************************/

void Solver::add(Field2D &v, Field2D &F_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Adding 2D field: Solver::add(%s)", name);
#endif

  if(initialised) {
    bout_error("Error: Cannot add to solver after initialisation\n");
  }
  
  // Set boundary conditions
  v.setBoundary(name);

  VarStr<Field2D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &F_v;
  d.name = string(name);

  f2d.push_back(d);

#ifdef TRACK
  var.name = name;
#endif

  /// Generate initial perturbations.
  /// NOTE: This could be done in init, but this would prevent the user
  ///       from modifying the initial perturbation (e.g. to prevent unphysical situations)
  ///       before it's loaded into the solver. If restarting, this perturbation
  ///       will be over-written anyway
  initial_profile(name, v);
  v.applyBoundary();

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::add(Field3D &v, Field3D &F_v, const char* name)
{ 
#ifdef CHECK
  int msg_point = msg_stack.push("Adding 3D field: Solver::add(%s)", name);
#endif

  if(initialised) {
    bout_error("Error: Cannot add to solver after initialisation\n");
  }

  // Set boundary conditions
  v.setBoundary(name);

  if(mesh->StaggerGrids && (v.getLocation() != CELL_CENTRE)) {
    output.write("\tVariable %s shifted to %s\n", name, strLocation(v.getLocation()));
    F_v.setLocation(v.getLocation()); // Make sure both at the same location
  }

  VarStr<Field3D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &F_v;
  d.location = v.getLocation();
  d.name = string(name);
  
  f3d.push_back(d);

#ifdef TRACK
  var.name = name;
#endif

  initial_profile(name, v);
  v.applyBoundary(); // Make sure initial profile obeys boundary conditions

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::add(Vector2D &v, Vector2D &F_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Adding 2D vector: Solver::add(%s)", name);
#endif

  if(initialised) {
    bout_error("Error: Cannot add to solver after initialisation\n");
  }

  // Set boundary conditions
  v.setBoundary(name);
  
  VarStr<Vector2D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &F_v;
  d.covariant = v.covariant;
  d.name = string(name);

  v2d.push_back(d);

  /// NOTE: No initial_profile call, because this will be done for each
  ///       component individually.
  
  /// Add suffix, depending on co- /contravariance
  if(v.covariant) {
    add(v.x, F_v.x, (d.name+"_x").c_str());
    add(v.y, F_v.y, (d.name+"_y").c_str());
    add(v.z, F_v.z, (d.name+"_z").c_str());
  }else {
    add(v.x, F_v.x, (d.name+"x").c_str());
    add(v.y, F_v.y, (d.name+"y").c_str());
    add(v.z, F_v.z, (d.name+"z").c_str());
  }
  
  /// Make sure initial profile obeys boundary conditions
  v.applyBoundary();

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::add(Vector3D &v, Vector3D &F_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Adding 3D vector: Solver::add(%s)", name);
#endif

  if(initialised) {
    bout_error("Error: Cannot add to solver after initialisation\n");
  }
  
  // Set boundary conditions
  v.setBoundary(name);

  VarStr<Vector3D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &F_v;
  d.covariant = v.covariant;
  d.name = string(name);
  
  v3d.push_back(d);

  // Add suffix, depending on co- /contravariance
  if(v.covariant) {
    add(v.x, F_v.x, (d.name+"_x").c_str());
    add(v.y, F_v.y, (d.name+"_y").c_str());
    add(v.z, F_v.z, (d.name+"_z").c_str());
  }else {
    add(v.x, F_v.x, (d.name+"x").c_str());
    add(v.y, F_v.y, (d.name+"y").c_str());
    add(v.z, F_v.z, (d.name+"z").c_str());
  }

  v.applyBoundary();

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * Constraints
 **************************************************************************/

void Solver::constraint(Field2D &v, Field2D &C_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 2D scalar: Solver::constraint(%s)", name);
#endif

  if(!has_constraints)
    bout_error("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    bout_error("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    bout_error("WARNING: Constraint requested for variable with NULL name\n");
  
  VarStr<Field2D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.name = string(name);

  f2d.push_back(d);

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::constraint(Field3D &v, Field3D &C_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 3D scalar: Solver::constraint(%s)", name);
#endif

  if(!has_constraints)
    bout_error("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    bout_error("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    bout_error("WARNING: Constraint requested for variable with NULL name\n");

  VarStr<Field3D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.location = v.getLocation();
  d.name = string(name);
  
  f3d.push_back(d);

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::constraint(Vector2D &v, Vector2D &C_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 2D vector: Solver::constraint(%s)", name);
#endif

  if(!has_constraints)
    bout_error("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    bout_error("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    bout_error("WARNING: Constraint requested for variable with NULL name\n");
    
  VarStr<Vector2D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.covariant = v.covariant;
  d.name = string(name);
  
  v2d.push_back(d);

  // Add suffix, depending on co- /contravariance
  if(v.covariant) {
    constraint(v.x, C_v.x, (d.name+"_x").c_str());
    constraint(v.y, C_v.y, (d.name+"_x").c_str());
    constraint(v.z, C_v.z, (d.name+"_x").c_str());
  }else {
    constraint(v.x, C_v.x, (d.name+"x").c_str());
    constraint(v.y, C_v.y, (d.name+"x").c_str());
    constraint(v.z, C_v.z, (d.name+"x").c_str());
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::constraint(Vector3D &v, Vector3D &C_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 3D vector: Solver::constraint(%s)", name);
#endif

  if(!has_constraints)
    bout_error("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    bout_error("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    bout_error("WARNING: Constraint requested for variable with NULL name\n");

  VarStr<Vector3D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.covariant = v.covariant;
  d.name = string(name);
  
  v3d.push_back(d);

  // Add suffix, depending on co- /contravariance
  if(v.covariant) {
    constraint(v.x, C_v.x, (d.name+"_x").c_str());
    constraint(v.y, C_v.y, (d.name+"_x").c_str());
    constraint(v.z, C_v.z, (d.name+"_x").c_str());
  }else {
    constraint(v.x, C_v.x, (d.name+"x").c_str());
    constraint(v.y, C_v.y, (d.name+"x").c_str());
    constraint(v.z, C_v.z, (d.name+"x").c_str());
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * Initialisation
 **************************************************************************/

int Solver::init(rhsfunc f, int argc, char **argv, bool restarting, int nout, BoutReal tstep)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Solver::init()");
#endif

  if(initialised)
    bout_error("ERROR: Solver is already initialised\n");

  phys_run = f;

  output.write("Initialising solver\n");
  
  /// GET GLOBAL OPTIONS
  options.setSection("");

  options.get("archive", archive_restart, -1);
/*    archive_restart = -1; // Not archiving restart files*/

  if(archive_restart > 0) {
    output.write("Archiving restart files every %d iterations\n",
        archive_restart);
  }

  /// Get restart file extension
  string dump_ext, restart_ext;

  options.get<string>("dump_format", dump_ext, DEFAULT_FILE_EXT);
  
  options.get<string>("restart_format", restart_ext, dump_ext);

  /// Set the restart file format
  restart.setFormat(data_format(restart_ext.c_str()));
  restartext = string(restart_ext);

  /// Add basic variables to the restart file
  restart.add(simtime,  "tt",    0);
  restart.add(iteration, "hist_hi", 0);
  
  MPI_Comm_size(MPI_COMM_WORLD, &NPES);
  MPI_Comm_rank(MPI_COMM_WORLD, &MYPE);
  
  restart.add(NPES, "NPES", 0);
  restart.add(mesh->NXPE, "NXPE", 0);

  /// Add variables to the restart and dump files.
  /// NOTE: Since vector components are already in the field arrays,
  ///       only loop over scalars, not vectors
  for(vector< VarStr<Field2D> >::iterator it = f2d.begin(); it != f2d.end(); it++) {
    // Add to restart file (not appending)
    restart.add(*(it->var), it->name.c_str(), 0);
    
    // Add to dump file (appending)
    dump.add(*(it->var), it->name.c_str(), 1);
    
    /// NOTE: Initial perturbations have already been set in add()
    
    /// Make sure boundary condition is satisfied
    it->var->applyBoundary();
  }  
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    // Add to restart file (not appending)
    restart.add(*(it->var), it->name.c_str(), 0);
    
    // Add to dump file (appending)
    dump.add(*(it->var), it->name.c_str(), 1);
    
    /// Make sure boundary condition is satisfied
    it->var->applyBoundary();
  }

  if(restarting) {
    /// Load state from the restart file
    
    // Copy processor numbers for comparison after. Very useful for checking
    // that the restart file is for the correct number of processors etc.
    int tmp_NP = NPES;
    int tmp_NX = mesh->NXPE;
    
#ifdef CHECK
    int msg_pt2 = msg_stack.push("Loading restart file");
#endif
    
    /// Load restart file
    if(restart.read("%s/BOUT.restart.%d.%s", restartdir.c_str(), MYPE, restartext.c_str()) != 0) {
      output.write("Error: Could not read restart file\n");
      return(2);
    }

    if(NPES == 0) {
      // Old restart file
      output.write("WARNING: Cannot verify processor numbers\n");
      NPES = tmp_NP;
      mesh->NXPE = tmp_NX;
    }else {
      // Check the processor numbers match
      if(NPES != tmp_NP) {
	output.write("ERROR: Number of processors (%d) doesn't match restart file number (%d)\n",
		     tmp_NP, NPES);
	return(1);
      }
      if(mesh->NXPE != tmp_NX) {
	output.write("ERROR: Number of X processors (%d) doesn't match restart file number (%d)\n",
		     tmp_NX, mesh->NXPE);
	return(1);
      }

      output.write("Restarting at iteration %d, simulation time %e\n", iteration, simtime);
    }

#ifdef CHECK
    msg_stack.pop(msg_pt2);
#endif
    
  }else {
    // Not restarting
    simtime = 0.0; iteration = 0;
  }
  
  /// Mark as initialised. No more variables can be added
  initialised = true;

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}

void Solver::setRestartDir(const string &dir)
{
  restartdir = dir;
}

/**************************************************************************
 * Useful routines (protected)
 **************************************************************************/

int Solver::getLocalN()
{
  int n2d = n2Dvars();
  int n3d = n3Dvars();
  
  int ncz = mesh->ngz-1;
  int MYSUB = mesh->yend - mesh->ystart + 1;

  int local_N = (mesh->xend - mesh->xstart + 1) *
    (mesh->yend - mesh->ystart + 1)*(n2d + ncz*n3d); // NOTE: Not including extra toroidal point

  //////////// Find boundary regions ////////////
  
  // Y up
  RangeIter *xi = mesh->iterateBndryUpperY();
  for(xi->first(); !xi->isDone(); xi->next()) {
    local_N +=  (mesh->ngy - mesh->yend - 1) * (n2d + ncz * n3d);
  }
  delete xi;
  
  // Y down
  xi = mesh->iterateBndryLowerY();
  for(xi->first(); !xi->isDone(); xi->next()) {
    local_N +=  mesh->ystart * (n2d + ncz * n3d);
  }
  delete xi;
  
  // X inner
  if(mesh->firstX()) {
    local_N += mesh->xstart * MYSUB * (n2d + ncz * n3d);
    output.write("\tBoundary region inner X\n");
  }

  // X outer
  if(mesh->lastX()) {
    local_N += (mesh->ngx - mesh->xend - 1) * MYSUB * (n2d + ncz * n3d);
    output.write("\tBoundary region outer X\n");
  }
  
  return local_N;
}

Solver* Solver::Create()
{  
  return SolverFactory::getInstance()->createSolver();
}

int Solver::run_rhs(BoutReal t)
{
  int status = (*phys_run)(t);

  // Make sure vectors in correct basis
  for(int i=0;i<v2d.size();i++) {
    if(v2d[i].covariant) {
      v2d[i].F_var->toCovariant();
    }else
      v2d[i].F_var->toContravariant();
  }
  for(int i=0;i<v3d.size();i++) {
    if(v3d[i].covariant) {
      v3d[i].F_var->toCovariant();
    }else
      v3d[i].F_var->toContravariant();
  }

  // Make sure 3D fields are at the correct cell location
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    if((*it).location != ((*it).F_var)->getLocation()) {
      //output.write("SOLVER: Interpolating\n");
      *((*it).F_var) = interp_to(*((*it).F_var), (*it).location);
    }
  }

  // Apply boundary conditions to the time-derivatives
  for(vector< VarStr<Field2D> >::iterator it = f2d.begin(); it != f2d.end(); it++) {
    it->var->applyTDerivBoundary();
  }
  
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    it->var->applyTDerivBoundary();
  }

#ifdef CHECK
  msg_stack.push("Solver checking time derivatives");
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++)
    it->F_var->checkData();
  msg_stack.pop();
#endif
  
  return status;
}
