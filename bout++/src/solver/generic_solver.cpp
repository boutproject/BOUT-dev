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
#include "generic_solver.h"

#include "initialprofiles.h"
#include "boundary.h"
#include "interpolation.h"

/**************************************************************************
 * Constructor
 **************************************************************************/

GenericSolver::GenericSolver() {
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

void GenericSolver::add(Field2D &v, Field2D &F_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Adding 2D field: GenericSolver::add(%s)", name);
#endif

  if(initialised) {
    bout_error("Error: Cannot add to solver after initialisation\n");
  }
  
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

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void GenericSolver::add(Field3D &v, Field3D &F_v, const char* name)
{ 
#ifdef CHECK
  int msg_point = msg_stack.push("Adding 3D field: GenericSolver::add(%s)", name);
#endif

  if(initialised) {
    bout_error("Error: Cannot add to solver after initialisation\n");
  }

  if(mesh->StaggerGrids && (v.getLocation() != CELL_CENTRE)) {
    output.write("\tVariable %s shifted to %s\n", name, strLocation(v.getLocation()));
    F_v.setLocation(v.getLocation()); // Make sure both at the same location
  }
  
  // Print the boundary conditions
  print_boundary(name);

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

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void GenericSolver::add(Vector2D &v, Vector2D &F_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Adding 2D vector: GenericSolver::add(%s)", name);
#endif

  if(initialised) {
    bout_error("Error: Cannot add to solver after initialisation\n");
  }

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

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void GenericSolver::add(Vector3D &v, Vector3D &F_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Adding 3D vector: GenericSolver::add(%s)", name);
#endif

  if(initialised) {
    bout_error("Error: Cannot add to solver after initialisation\n");
  }

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

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * Constraints
 **************************************************************************/

void GenericSolver::constraint(Field2D &v, Field2D &C_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 2D scalar: GenericSolver::constraint(%s)", name);
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

void GenericSolver::constraint(Field3D &v, Field3D &C_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 3D scalar: GenericSolver::constraint(%s)", name);
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

void GenericSolver::constraint(Vector2D &v, Vector2D &C_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 2D vector: GenericSolver::constraint(%s)", name);
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

void GenericSolver::constraint(Vector3D &v, Vector3D &C_v, const char* name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 3D vector: GenericSolver::constraint(%s)", name);
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

int GenericSolver::init(rhsfunc f, int argc, char **argv, bool restarting, int nout, real tstep)
{
#ifdef CHECK
  int msg_point = msg_stack.push("GenericSolver::init()");
#endif
  
  if(initialised)
    bout_error("ERROR: Solver is already initialised\n");

  output.write("Initialising solver\n");
  
  /// GET GLOBAL OPTIONS
  options.setSection(NULL);

  if(options.getInt("archive", archive_restart)) {
    archive_restart = -1; // Not archiving restart files
  }else if(archive_restart > 0) {
    output.write("Archiving restart files every %d iterations\n",
        archive_restart);
  }

  /// Get restart file extension
  const char *dump_ext, *restart_ext;
  if((dump_ext = options.getString("dump_format")) == NULL) {
    // Set default extension
    dump_ext = DEFAULT_FILE_EXT;
  }
  
  if((restart_ext = options.getString("restart_format")) == NULL) 
    restart_ext = dump_ext;

  /// Set the restart file format
  restart.setFormat(data_format(restart_ext));
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
  }  
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    // Add to restart file (not appending)
    restart.add(*(it->var), it->name.c_str(), 0);
    
    // Add to dump file (appending)
    dump.add(*(it->var), it->name.c_str(), 1);
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

void GenericSolver::setRestartDir(const string &dir)
{
  restartdir = dir;
}

/**************************************************************************
 * Useful routines (protected)
 **************************************************************************/

int GenericSolver::getLocalN()
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
  if(mesh->first_x()) {
    local_N += mesh->xstart * MYSUB * (n2d + ncz * n3d);
    output.write("\tBoundary region inner X\n");
  }

  // X outer
  if(mesh->last_x()) {
    local_N += (mesh->ngx - mesh->xend - 1) * MYSUB * (n2d + ncz * n3d);
    output.write("\tBoundary region outer X\n");
  }
  
  return local_N;
}
