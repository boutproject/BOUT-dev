/**************************************************************************
 * Implementation of the Mesh class, handling input files compatible with
 * BOUT / BOUT-06.
 *
 * Changelog
 * ---------
 *
 * 2010-05 Ben Dudson <bd512@york.ac.uk>
 *      * Initial version, adapted from grid.cpp and topology.cpp
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

#include "boutmesh.h"

#include "globals.h"
#include "utils.h"
#include "fft.h"

#include "dcomplex.h"

#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

BoutMesh::~BoutMesh()
{
  // Delete the communication handles
  clear_handles();
  
  // Delete the boundary regions
  for(vector<BoundaryRegion*>::iterator it = boundary.begin(); it != boundary.end(); it++)
    delete (*it);
}

int BoutMesh::load()
{
#ifdef CHECK
  int msg = msg_stack.push("BoutMesh::load()");
#endif
  
  output << "Loading mesh" << endl;

  //////////////
  // Number of processors
  
  MPI_Comm_size(MPI_COMM_WORLD, &NPES);
  MPI_Comm_rank(MPI_COMM_WORLD, &MYPE);

  //////////////
  // Grid sizes
  
  if(get(nx, "nx"))
    return 1;
  
  if(get(ny, "ny"))
    return 1;
  
  output << "\tGrid size: " << nx << " by " << ny << endl;

  options.setSection("");
  options.get("MXG", MXG, 2);
  options.get("MYG", MYG, 2);
  
  options.get("NXPE", NXPE, 1); // Decomposition in the radial direction
  if((NPES % NXPE) != 0) {
    output.write("Error: Number of processors (%d) not divisible by NPs in x direction (%d). Aborting\n",
		 NPES, NXPE);
    return(1);
  }

  NYPE = NPES / NXPE;
  
  /// Get X and Y processor indices
  PE_YIND = MYPE / NXPE;
  PE_XIND = MYPE % NXPE;
  
  // Work out other grid size quantities

  /// MXG at each end needed for edge boundary regions
  MX = nx - 2*MXG;
  
  /// Split MX points between NXPE processors
  MXSUB = MX / NXPE;
  if((MX % NXPE) != 0) {
    output.write("\tERROR: Cannot split %d X points equally between %d processors\n",
		 MX, NXPE);
    return false;
  }

  /// NOTE: No grid data reserved for Y boundary cells - copy from neighbours
  MY = ny;
  MYSUB = MY / NYPE;
  if((MY % NYPE) != 0) {
    output.write("\tERROR: Cannot split %d Y points equally between %d processors\n",
		 MY, NYPE);
    return false;
  }
  
  /// Get mesh options
  options.setSection(""); // Global options
  int MZ;
  OPTION(MZ,           65);
  if(!is_pow2(MZ-1)) {
    if(is_pow2(MZ)) {
      MZ++;
      output.write("WARNING: Number of toroidal points increased to %d\n", MZ);
    }else {
      output.write("Error: Number of toroidal points must be 2^n + 1");
      return 1;
    }
  }
  OPTION(TwistShift,   false);
  OPTION(TwistOrder,   0);
  OPTION(ShiftOrder,   0);
  OPTION(ShiftXderivs, false);
  OPTION(IncIntShear,  false);
  OPTION(BoundaryOnCell, false); // Determine location of boundary
  OPTION(StaggerGrids,   false); // Stagger grids
  
  OPTION(async_send, false); // Whether to use asyncronous sends

  if(ShiftXderivs) {
    output.write("Using shifted X derivatives. Interpolation: ");
    if(mesh->ShiftOrder == 0) {
      output.write("FFT\n");
    }else
      output.write("%d-point\n", mesh->ShiftOrder);
  }

  options.get("zperiod",   zperiod,      1);
  if(zperiod == 1) {
    options.get("ZMIN",         ZMIN,         0.0);
    options.get("ZMAX",         ZMAX,         1.0);
    
    zperiod = ROUND(1.0 / (ZMAX - ZMIN));
  }else {
    ZMIN = 0.0;
    ZMAX = 1.0 / (double) zperiod;
  }

  if(TwistShift) {
    output.write("Applying Twist-Shift condition. Interpolation: ");
    if(mesh->TwistOrder == 0) {
      output.write("FFT\n");
    }else
      output.write("%d-point\n", mesh->TwistOrder);
  }
  
  /// Number of grid cells is ng* = M*SUB + guard/boundary cells
  ngx = MXSUB + 2*MXG;
  ngy = MYSUB + 2*MYG;
  ngz = MZ;
  
  // Set local index ranges
  
  xstart = MXG;
  xend = MXG + MXSUB - 1;

  ystart = MYG;
  yend = MYG + MYSUB - 1;
  
  ///////////////////// TOPOLOGY //////////////////////////
  
  // separatrix location
  if(get(ixseps1, "ixseps1")) {
    ixseps1 = ngx;
    output.write("\tWARNING: Separatrix location 'ixseps1' not found. Setting to %d\n", ixseps1);
  }
  if(get(ixseps2, "ixseps2")) {
    ixseps2 = ngx;
    output.write("\tWARNING: Separatrix location 'ixseps2' not found. Setting to %d\n", ixseps2);
  }
  if(get(jyseps1_1,"jyseps1_1")) {
    jyseps1_1 = -1;
    output.write("\tWARNING: Branch-cut 'jyseps1_1' not found. Setting to %d\n", jyseps1_1);
  }
  if(get(jyseps1_2,"jyseps1_2")) {
    jyseps1_2 = ny/2;
    output.write("\tWARNING: Branch-cut 'jyseps1_2' not found. Setting to %d\n", jyseps1_2);
  }
  if(get(jyseps2_1,"jyseps2_1")) {
    jyseps2_1 = jyseps1_2;
    output.write("\tWARNING: Branch-cut 'jyseps2_1' not found. Setting to %d\n", jyseps2_1);
  }
  if(get(jyseps2_2,"jyseps2_2")) {
    jyseps2_2 = ny-1;
    output.write("\tWARNING: Branch-cut 'jyseps2_2' not found. Setting to %d\n", jyseps2_2);
  }

  if(get(ny_inner,"ny_inner")) {
    ny_inner = jyseps2_1;
    output.write("\tWARNING: Number of inner y points 'ny_inner' not found. Setting to %d\n", ny_inner);
  }
  
  /// Call topology to set layout of grid
  topology();
  
  ///////////////// DIFFERENCING QUANTITIES ///////////////
  
  if(get(dx, "dx")) {
    output.write("\tWARNING: differencing quantity 'dx' not found. Set to 1.0\n");
    dx = 1.0;
  }
  if(get(dy, "dy")) {
    output.write("\tWARNING: differencing quantity 'dy' not found. Set to 1.0\n");
    dy = 1.0;
  }
  
  if(non_uniform) {
    // Read correction for non-uniform meshes
    if(get(d2x, "d2x")) {
      output.write("\tWARNING: differencing quantity 'd2x' not found. Set to 0.0\n");
      d2x = 0.0;
    }
    if(get(d2y, "d2y")) {
      output.write("\tWARNING: differencing quantity 'd2y' not found. Set to 0.0\n");
      d2y = 0.0;
    }
  }
  
  zlength = (ZMAX-ZMIN)*TWOPI;
  dz = zlength/(ngz-1);

  ///////////////// DIFFERENTIAL GEOMETRY /////////////////
  
  // Diagonal components of metric tensor g^{ij} (default to 1)
  get(g11, "g11", 1.0);
  get(g22, "g22", 1.0);
  get(g33, "g33", 1.0);
  
  // Off-diagonal elements. Default to 0
  get(g12, "g12", 0.0);
  get(g13, "g13", 0.0);
  get(g23, "g23", 0.0);
  
  // Check input metrics
  if((!finite(g11)) || (!finite(g22)) || (!finite(g33))) {
    output.write("\tERROR: Diagonal metrics are not finite!\n");
    exit(1);
  }
  if((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    output.write("\tERROR: Diagonal metrics are negative!\n");
  }
  if((!finite(g12)) || (!finite(g13)) || (!finite(g23))) {
    output.write("\tERROR: Off-diagonal metrics are not finite!\n");
    exit(1);
  }

  /// Set shift for radial derivatives
  if(get(zShift, "zShift")) {
    output.write("\tWARNING: Z shift for radial derivatives not found\n");
    ShiftTorsion = zShift = 0.0;
  }else if(get(ShiftTorsion, "ShiftTorsion")) {
    output.write("\tWARNING: No Torsion specified for zShift. Derivatives may not be correct\n");
    ShiftTorsion = 0.0;
  }
  
  if(mesh->IncIntShear) {
    if(get(IntShiftTorsion, "IntShiftTorsion")) {
      output.write("\tWARNING: No Integrated torsion specified\n");
      IntShiftTorsion = 0.0;
    }
  }
  // Allocate some memory for twist-shift

  ShiftAngle  = rvector(ngx);

  // Try to read the shift angle from the grid file
  // NOTE: All processors should know the twist-shift angle (for invert_parderiv)
  GridDataSource* s = findSource("ShiftAngle");
  if(s) {
    s->open("ShiftAngle");
    s->setOrigin(XGLOBAL(0));
    if(!s->fetch(ShiftAngle,  "ShiftAngle", ngx)) {
      output.write("\tWARNING: Twist-shift angle 'ShiftAngle' not found. Setting to zero\n");
      for(int i=0;i<ngx;i++)
        ShiftAngle[i] = 0.0;
    }
    s->close();
  }else {
    output.write("\tWARNING: Twist-shift angle 'ShiftAngle' not found. Setting from zShift\n");

    if(YPROC(jyseps2_2) == PE_YIND) {
      for(int i=0;i<ngx;i++)
	ShiftAngle[i] = zShift[i][MYG+MYSUB-1] - zShift[i][MYG+MYSUB]; // Jump across boundary
   
    }else if(YPROC(jyseps1_1+1) == PE_YIND) {
      for(int i=0;i<mesh->ngx;i++)
	ShiftAngle[i] = mesh->zShift[i][MYG-1] - mesh->zShift[i][MYG]; // Jump across boundary
    }
    
    // In the core, need to set ShiftAngle everywhere for ballooning initial condition
    MPI_Group groupw;
    MPI_Comm_group(MPI_COMM_WORLD, &groupw); // Group of all processors
    
    int *ranks = new int[NYPE];
    int npcore = 0;
    for(int p = YPROC(jyseps1_1+1); p <= YPROC(jyseps2_2);p++) {
      ranks[npcore] = PROC_NUM(PE_XIND, p);
      npcore++;
    }
    
    MPI_Group grp;
    MPI_Group_incl(groupw, npcore, ranks, &grp); // Create group
    
    MPI_Comm core_comm;
    MPI_Comm_create(MPI_COMM_WORLD, grp, &core_comm); // Create communicator
    
    delete[] ranks;
    
    if(MYPE_IN_CORE)
      MPI_Bcast(ShiftAngle, ngx, PVEC_REAL_MPI_TYPE, npcore-1, core_comm);
    
  }

  /// Can have twist-shift in the private flux regions too
  bool twistshift_pf;
  options.setSection("");
  OPTION(twistshift_pf, false);
  if(twistshift_pf) {
    output << "Adding twist-shift in lower PF region" << endl;
    // Lower PF. Note by default no Twist-Shift used here, so need to switch on
    if(YPROC(jyseps1_1) == PE_YIND) {
      for(int i=0;i<ngx;i++) {
	ShiftAngle[i] = zShift[i][MYG+MYSUB-1] - zShift[i][MYG+MYSUB]; // Jump across boundary
      }
      TS_up_in = true; // Switch on twist-shift
      
    }else if(YPROC(jyseps2_2+1) == PE_YIND) {
      for(int i=0;i<mesh->ngx;i++) {
	ShiftAngle[i] = mesh->zShift[i][MYG-1] - mesh->zShift[i][MYG]; // Jump across boundary
      }
      TS_down_in = true;
    }
  }

  /// Calculate contravariant metric components
  if(calcCovariant())
    return 1;

  /// Calculate Jacobian and Bxy
  if(jacobian())
    return 1;
  
  // Attempt to read J from the grid file
  Field2D Jcalc = mesh->J;
  if(mesh->get(mesh->J, "J")) {
    output.write("\tWARNING: Jacobian 'J' not found. Calculating from metric tensor\n");
    mesh->J = Jcalc;
  }else {
    // Compare calculated and loaded values  
    output.write("\tMaximum difference in J is %e\n", max(abs(mesh->J - Jcalc)));
    
    // Re-evaluate Bxy using new J
    Bxy = sqrt(mesh->g_22)/mesh->J;
  }

  // Attempt to read Bxy from the grid file
  Field2D Bcalc = Bxy;
  if(mesh->get(Bxy, "Bxy")) {
    output.write("\tWARNING: Magnitude of B field 'Bxy' not found. Calculating from metric tensor\n");
    Bxy = Bcalc;
  }else {
    output.write("\tMaximum difference in Bxy is %e\n", max(abs(Bxy - Bcalc)));
    // Check Bxy
    if(!finite(Bxy)) {
      output.write("\tERROR: Bxy not finite everywhere!\n");
      exit(1);
    }
  }
  
  /// Calculate Christoffel symbols
  if(geometry()) {
    output << "  Differential geometry failed\n";
    return 1;
  }
  
  //////////////////////////////////////////////////////
  /// Communicators for Y gather/scatter
  
  //MPI_Comm comm_inner, comm_middle, comm_outer;
  
  MPI_Group group_world;
  MPI_Comm_group(MPI_COMM_WORLD, &group_world); // Get the entire group
  
  MPI_Group group;
  MPI_Group group_tmp1, group_tmp2;
  
  int proc[3]; // Processor range
  proc[2] = NXPE; // Stride in processor rank
  
  MPI_Comm comm_tmp;
  
  // Outer SOL regions
  if(jyseps1_2 == jyseps2_1) {
    // Single-null. All processors with same PE_XIND
    
    for(int i=0;i<NXPE;i++) {
      proc[0] = PROC_NUM(i, 0);
      proc[1] = PROC_NUM(i, NYPE-1);
      MPI_Group_range_incl(group_world, 1, &proc, &group);
      MPI_Comm_create(MPI_COMM_WORLD, group, &comm_tmp);
      if(i == PE_XIND) {
	// Should be part of this communicator
	if(comm_tmp == MPI_COMM_NULL) {
	  // error
	  bout_error("Single null outer SOL not correct\n");
	}
	comm_outer = comm_tmp;
      }else if(comm_tmp != MPI_COMM_NULL) {
	// Not part of this communicator so should be NULL
	bout_error("Single null outer SOL not correct\n");
      }
      MPI_Group_free(&group);
    }
  }else {
    // Double null
    
    for(int i=0;i<NXPE;i++) {
      // Inner SOL
      proc[0] = PROC_NUM(i, 0);
      proc[1] = PROC_NUM(i, YPROC(ny_inner-1));
      MPI_Comm_create(MPI_COMM_WORLD, group, &comm_tmp);
      if(comm_tmp != MPI_COMM_NULL)
	comm_outer = comm_tmp;
      MPI_Group_free(&group);

      // Outer SOL
      proc[0] = PROC_NUM(i, YPROC(ny_inner));
      proc[1] = PROC_NUM(i, NYPE-1);
      MPI_Group_range_incl(group_world, 1, &proc, &group);
      MPI_Comm_create(MPI_COMM_WORLD, group, &comm_tmp);
      if(comm_tmp != MPI_COMM_NULL)
	comm_outer = comm_tmp;
      MPI_Group_free(&group);
    }
  }

  for(int i=0;i<NXPE;i++) {
    // Lower PF region

    if((jyseps1_1 >= 0) || (jyseps2_2 < ny)) {
      // A lower PF region exists

      if(jyseps1_1 >= 0) {
	proc[0] = PROC_NUM(i, 0);
	proc[1] = PROC_NUM(i, YPROC(jyseps1_1));
	//output << "PF1 "<< proc[0] << ", " << proc[1] << endl;
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
      }else
	group_tmp1 = MPI_GROUP_EMPTY;
      
      if(jyseps2_2+1 < ny) {
	proc[0] = PROC_NUM(i, YPROC(jyseps2_2+1));
	proc[1] = PROC_NUM(i, NYPE-1);
	//output << "PF2 "<< proc[0] << ", " << proc[1] << endl;
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
      }else
	group_tmp2 = MPI_GROUP_EMPTY;
      
      MPI_Group_union(group_tmp1, group_tmp2, &group);
      MPI_Comm_create(MPI_COMM_WORLD, group, &comm_tmp);
      if(comm_tmp != MPI_COMM_NULL) {
	comm_inner = comm_tmp;
	if(ixseps_lower == ixseps_outer) {
	  // Between the separatrices is still in the PF region
	  comm_middle = comm_inner;
	}else
	  comm_middle = comm_outer;
      }
      MPI_Group_free(&group);
    }

    if(jyseps2_1 != jyseps1_2) {
      // Upper PF region
      // Note need to order processors so that a continuous surface is formed
      
      proc[0] = PROC_NUM(i, YPROC(ny_inner));
      proc[1] = PROC_NUM(i, YPROC(jyseps1_2));
      //output << "PF3 "<< proc[0] << ", " << proc[1] << endl;
      MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
      proc[0] = PROC_NUM(i, YPROC(jyseps2_1+1));
      proc[1] = PROC_NUM(i, YPROC(ny_inner-1));
      //output << "PF4 "<< proc[0] << ", " << proc[1] << endl;
      MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
      MPI_Group_union(group_tmp1, group_tmp2, &group);
      MPI_Comm_create(MPI_COMM_WORLD, group, &comm_tmp);
      if(comm_tmp != MPI_COMM_NULL) {
	comm_inner = comm_tmp;
	if(ixseps_upper == ixseps_outer) {
	  comm_middle = comm_inner;
	}else
	  comm_middle = comm_outer;
      }
      MPI_Group_free(&group);
    }
    
    // Core region
    
    proc[0] = PROC_NUM(i, YPROC(jyseps1_1+1));
    proc[1] = PROC_NUM(i, YPROC(jyseps2_1));
    //output << "CORE1 "<< proc[0] << ", " << proc[1] << endl;
    MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
    proc[0] = PROC_NUM(i, YPROC(jyseps1_2+1));
    proc[1] = PROC_NUM(i, YPROC(jyseps2_2));
    //output << "CORE2 "<< proc[0] << ", " << proc[1] << endl;
    MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
    MPI_Group_union(group_tmp1, group_tmp2, &group);
    MPI_Comm_create(MPI_COMM_WORLD, group, &comm_tmp);
    if(comm_tmp != MPI_COMM_NULL) {
      comm_inner = comm_tmp;
      
      if(ixseps_inner == ixseps_outer)
	comm_middle = comm_inner;
    }
  }
  
  if(ixseps_inner != ixseps_outer) {
    // Need to handle unbalanced double-null case
    
    if(ixseps_upper > ixseps_lower) {
      // middle is connected to the bottom
	
      for(int i=0;i<NXPE;i++) {
	proc[0] = PROC_NUM(i, 0);
	proc[1] = PROC_NUM(i, YPROC(jyseps2_1));
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
	proc[0] = PROC_NUM(i, YPROC(jyseps1_2+1));
	proc[1] = PROC_NUM(i, NYPE-1);
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
	MPI_Group_union(group_tmp1, group_tmp2, &group);
	MPI_Comm_create(MPI_COMM_WORLD, group, &comm_tmp);
	if(comm_tmp != MPI_COMM_NULL)
	  comm_middle = comm_tmp;
      }
    }else {
      // middle is connected to the top
	
      for(int i=0;i<NXPE;i++) {
	proc[0] = PROC_NUM(i, YPROC(ny_inner));
	proc[1] = PROC_NUM(i, YPROC(jyseps2_2));
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
	proc[0] = PROC_NUM(i, YPROC(jyseps1_1+1));
	proc[1] = PROC_NUM(i, YPROC(ny_inner-1));
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
	MPI_Group_union(group_tmp1, group_tmp2, &group);
	MPI_Comm_create(MPI_COMM_WORLD, group, &comm_tmp);
	if(comm_tmp != MPI_COMM_NULL)
	  comm_middle = comm_tmp;
      }
    }
  }
  // Now have communicators for all regions.
  
  //////////////////////////////////////////////////////
  // Boundary regions
  if(PE_XIND == 0) {
    // Inner either core or PF
    
    int yg = YGLOBAL(MYG); // Get a global index in this processor
    
    if( ((yg > jyseps1_1) && (yg <= jyseps2_1)) ||
	((yg > jyseps1_2) && (yg <= jyseps2_2)) ) {
      // Core
      boundary.push_back(new BoundaryRegionXIn("core", ystart, yend));
    }else {
      // PF region
      boundary.push_back(new BoundaryRegionXIn("pf", ystart, yend));
    }
  }
  if(PE_XIND == (NXPE-1)){
    // Outer SOL
    boundary.push_back(new BoundaryRegionXOut("sol", ystart, yend));
  }
  
  if((UDATA_INDEST < 0) && (UDATA_XSPLIT > xstart))
    boundary.push_back(new BoundaryRegionYUp("target", xstart, UDATA_XSPLIT-1));
  if((UDATA_OUTDEST < 0) && (UDATA_XSPLIT <= xend))
    boundary.push_back(new BoundaryRegionYUp("target", UDATA_XSPLIT, xend));
  
  if((DDATA_INDEST < 0) && (DDATA_XSPLIT > xstart))
    boundary.push_back(new BoundaryRegionYDown("target", xstart, UDATA_XSPLIT-1));
  if((DDATA_OUTDEST < 0) && (DDATA_XSPLIT <= xend))
    boundary.push_back(new BoundaryRegionYDown("target", UDATA_XSPLIT, xend));

  if(!boundary.empty()) {
    output << "Boundary regions in this processor: ";
    for(vector<BoundaryRegion*>::iterator it=boundary.begin(); it != boundary.end(); it++) {
      output << (*it)->label << ", ";
    }
    output << endl;
  }else {
    output << "No boundary regions in this processor" << endl;
  }
  
  output.write("\tdone\n");
  
#ifdef CHECK
  msg_stack.pop(msg);
#endif
  
  return 0;
}

/*****************************************************************************
 * get routines
 *****************************************************************************/

/// Get an integer
int BoutMesh::get(int &ival, const char *name)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Loading integer: BoutMesh::get(int, %s)", name);
#endif

  GridDataSource* s = findSource(name);
  if(s == NULL) {
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 1;
  }
  
  s->open(name);
  bool success = s->fetch(&ival, name);
  s->close();
  
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  if(!success) {
    return 2;
  }
  return 0;
}

/// A BoutReal number
int BoutMesh::get(BoutReal &rval, const char *name)
{
  GridDataSource* s = findSource(name);
  if(s == NULL)
    return 1;
  
  s->open(name);
  bool success = s->fetch(&rval, name);
  s->close();
  
  if(!success)
    return 2;
  return 0;
}

int BoutMesh::get(Field2D &var, const char *name, BoutReal def)
{
  if(name == NULL)
    return 1;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("Loading 2D field: BoutMesh::get(Field2D, %s)", name);
#endif
  
  GridDataSource *s = findSource(name);
  if(s == NULL) {
    output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name, def);
    var = def;
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 2;
  }
  
  BoutReal **data;
  int jy;
  
  var = 0;// Makes sure the memory is allocated
  
  data = var.getData(); // Get a pointer to the data
  
  // Send an open signal to the source
  s->open(name);
  
  // Read in data for bulk of points
  
  if(readgrid_2dvar(s, name,
		    YGLOBAL(MYG), // Start reading at global index for y=MYG
		    MYG,          // Insert data starting from y=MYG
		    MYSUB,        // Length of data is MYSUB
		    0, ngx,       // All x indices (local indices)
		    data)) {
    output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name, def);
    var = def;
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 1;
  }

  // Read in data for upper boundary
  // lowest (smallest y) part of UDATA_*DEST processor domain

  if((UDATA_INDEST != -1) && (UDATA_XSPLIT > 0)) {
    // Inner data exists and has a destination 
    
    if(readgrid_2dvar(s, name,
		      (UDATA_INDEST/NXPE)*MYSUB, // the "bottom" (y=1) of the destination processor
		      MYSUB+MYG,            // the same as the upper guard cell
		      MYG,                  // Only one y point
		      0, UDATA_XSPLIT,    // Just the inner cells
		      data)) {
      output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name, def);
      var = def;
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 2;
    }

  }else if(UDATA_XSPLIT > 0) {
    // Inner part exists but has no destination => boundary
    for(jy=0;jy<MYG;jy++)
      cpy_2d_data(MYSUB+MYG-1-jy, MYSUB+MYG+jy, 0, UDATA_XSPLIT, data); // copy values from last index into guard cell
  }
  
  if((UDATA_OUTDEST != -1) && (UDATA_XSPLIT < ngx)) { 
    
    if(readgrid_2dvar(s, name,
		      (UDATA_OUTDEST / NXPE)*MYSUB,
		      MYSUB+MYG,
		      MYG,
		      UDATA_XSPLIT, ngx, // the outer cells
		      data)) {
      output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name, def);
      var = def;
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 2;
    }
  }else if(UDATA_XSPLIT < MX) {
    // Inner data exists, but has no destination
    for(jy=0;jy<MYG;jy++)
      cpy_2d_data(MYSUB+MYG-1-jy, MYSUB+MYG+jy, UDATA_XSPLIT, ngx, data);
  }
  
  // Read in data for lower boundary
  if((DDATA_INDEST != -1) && (DDATA_XSPLIT > 0)) {
    //output.write("Reading DDEST: %d\n", (DDATA_INDEST+1)*MYSUB -1);

    if(readgrid_2dvar(s, name,
		      ((DDATA_INDEST/NXPE)+1)*MYSUB - MYG, // The "top" of the destination processor
		      0,  // belongs in the lower guard cell
		      MYG,  // just one y point
		      0, DDATA_XSPLIT, // just the inner data
		      data)) {
      output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name, def);
      var = def;
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 2;
    }
    
    }else if(DDATA_XSPLIT > 0) {
      for(jy=0;jy<MYG;jy++)
        cpy_2d_data(MYG+jy, MYG-1-jy, 0, DDATA_XSPLIT, data);
    }
    if((DDATA_OUTDEST != -1) && (DDATA_XSPLIT < ngx)) {

      if(readgrid_2dvar(s, name,
		      ((DDATA_OUTDEST/NXPE)+1)*MYSUB - MYG,
		      0,
		      MYG,
		      DDATA_XSPLIT, ngx,
		      data)) {
	output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name, def);
	var = def;
#ifdef CHECK
	msg_stack.pop(msg_pos);
#endif
        return 2;
      }
  }else if(DDATA_XSPLIT < ngx) {
    for(jy=0;jy<MYG;jy++)
      cpy_2d_data(MYG+jy, MYG-1-jy, DDATA_XSPLIT, ngx, data);
  }

#ifdef TRACK
  var.name = copy_string(name);
#endif
   
  // Close source
  s->close();
  
#ifdef CHECK
  // Check that the data is ok
  var.checkData(true);
  
  msg_stack.pop(msg_pos);
#endif
  
  return 0;
}

int BoutMesh::get(Field2D &var, const string &name, BoutReal def)
{
  return get(var, name.c_str());
}

/// Load a 3D variable from the grid file
/*!
  Data stored as toroidal FFTs in BoutReal space at each X-Y point.
  In toroidal direction, array must have an odd number of points.
  Format is:

  DC, r1,i1, r2,i2, ... , rn,in

  with the BoutReal and imaginary parts of each (positive) frequency
  up to the nyquist frequency.
 */
int BoutMesh::get(Field3D &var, const char *name)
{
  if(name == NULL)
    return 1;
  
#ifdef CHECK
  msg_stack.push("Loading 3D field: BoutMesh::get(Field3D, %s)", name);
#endif
  
  GridDataSource *s = findSource(name);
  if(s == NULL) {
    output.write("\tWARNING: Could not read '%s' from grid. Setting to zero\n", name);
    var = 0.0;
#ifdef CHECK
    msg_stack.pop();
#endif
    return 2;
  }

  BoutReal ***data;
  int jy;

  var = 0.0; // Makes sure the memory is allocated

  data = var.getData(); // Get a pointer to the data

  // Send open signal to data source
  s->open(name);

  // Read in data for bulk of points
  if(readgrid_3dvar(s, name,
		    YGLOBAL(MYG), // Start reading at global index for y=MYG
		    MYG,          // Insert data starting from y=MYG
		    MYSUB,        // Length of data is MYSUB
		    0, ngx,       // All x indices (local indices)
		    data)) {
    output.write("\tWARNING: Could not read '%s' from grid. Setting to zero\n", name);
    var = 0.0;
#ifdef CHECK
    msg_stack.pop();
#endif
    return 1;
  }

  // Read in data for upper boundary
  // lowest (smallest y) part of UDATA_*DEST processor domain

  if((UDATA_INDEST != -1) && (UDATA_XSPLIT > 0)) {
    // Inner data exists and has a destination 
    
    if(readgrid_3dvar(s, name,
		      (UDATA_INDEST/NXPE)*MYSUB, // the "bottom" (y=1) of the destination processor
		      MYSUB+MYG,            // the same as the upper guard cell
		      MYG,                  // Only one y point
		      0, UDATA_XSPLIT,    // Just the inner cells
		      data)) {
      output.write("\tWARNING: Could not read '%s' from grid. Setting to zero\n", name);
      var = 0.0;
#ifdef CHECK
      msg_stack.pop();
#endif
      return 2;
    }

  }else if(UDATA_XSPLIT > 0) {
    // Inner part exists but has no destination => boundary
    for(jy=0;jy<MYG;jy++)
      cpy_3d_data(MYSUB+MYG-1-jy, MYSUB+MYG+jy, 0, UDATA_XSPLIT, data); // copy values from last index into guard cell
  }
  
  if((UDATA_OUTDEST != -1) && (UDATA_XSPLIT < ngx)) { 
    
    if(readgrid_3dvar(s, name,
		      (UDATA_OUTDEST / NXPE)*MYSUB,
		      MYSUB+MYG,
		      MYG,
		      UDATA_XSPLIT, ngx, // the outer cells
		      data)) {
      output.write("\tWARNING: Could not read '%s' from grid. Setting to zero\n", name);
      var = 0.0;
#ifdef CHECK
      msg_stack.pop();
#endif
      return 2;
    }
  }else if(UDATA_XSPLIT < MX) {
    // Inner data exists, but has no destination
    for(jy=0;jy<MYG;jy++)
      cpy_3d_data(MYSUB+MYG-1-jy, MYSUB+MYG+jy, UDATA_XSPLIT, ngx, data);
  }
  
  // Read in data for lower boundary
  if((DDATA_INDEST != -1) && (DDATA_XSPLIT > 0)) {
    //output.write("Reading DDEST: %d\n", (DDATA_INDEST+1)*MYSUB -1);

    if(readgrid_3dvar(s, name,
		      ((DDATA_INDEST/NXPE)+1)*MYSUB - MYG, // The "top" of the destination processor
		      0,  // belongs in the lower guard cell
		      MYG,  // just one y point
		      0, DDATA_XSPLIT, // just the inner data
		      data)) {
      output.write("\tWARNING: Could not read '%s' from grid. Setting to zero\n", name);
      var = 0.0;
#ifdef CHECK
      msg_stack.pop();
#endif
      return 2;
    }
		     
   }else if(DDATA_XSPLIT > 0) {
     for(jy=0;jy<MYG;jy++)
       cpy_3d_data(MYG+jy, MYG-1-jy, 0, DDATA_XSPLIT, data);
   }
  if((DDATA_OUTDEST != -1) && (DDATA_XSPLIT < ngx)) {

    if(readgrid_3dvar(s, name,
		      ((DDATA_OUTDEST/NXPE)+1)*MYSUB - MYG,
		      0,
		      MYG,
		      DDATA_XSPLIT, ngx,
		      data)) {
      output.write("\tWARNING: Could not read '%s' from grid. Setting to zero\n", name);
      var = 0.0;
#ifdef CHECK
      msg_stack.pop();
#endif
      return 2;
    }
  }else if(DDATA_XSPLIT < ngx) {
    for(jy=0;jy<MYG;jy++)
      cpy_3d_data(MYG+jy, MYG-1-jy, DDATA_XSPLIT, ngx, data);
  }
  
#ifdef CHECK
  // Check that the data is ok
  var.checkData(true);
  
  msg_stack.pop();
#endif
  
  return 0;
}

int BoutMesh::get(Field3D &var, const string &name)
{
  return get(var, name.c_str());
}

/****************************************************************
 *                 COMMUNICATIONS
 ****************************************************************/

const int IN_SENT_UP    = 0;
const int OUT_SENT_UP   = 1;
const int IN_SENT_DOWN  = 2;
const int OUT_SENT_DOWN = 3;
// X communication signals
const int IN_SENT_OUT = 4;
const int OUT_SENT_IN  = 5;

int BoutMesh::communicate(FieldGroup &g)
{
  comm_handle c = send(g);
  return wait(c);
}

void BoutMesh::post_receive(CommHandle &ch)
{
  BoutReal *inbuff;
  int len;
  
  /// Post receive data from above (y+1)

  len = 0;
  if(UDATA_INDEST != -1) {
    len = msg_len(ch.var_list, 0, UDATA_XSPLIT, 0, MYG);
    MPI_Irecv(ch.umsg_recvbuff,
	      len,
	      PVEC_REAL_MPI_TYPE,
	      UDATA_INDEST,
	      IN_SENT_DOWN,
	      MPI_COMM_WORLD,
	      &ch.request[0]);
  }
  if(UDATA_OUTDEST != -1) {
    inbuff = &ch.umsg_recvbuff[len]; // pointer to second half of the buffer
    MPI_Irecv(inbuff,
	      msg_len(ch.var_list, UDATA_XSPLIT, ngx, 0, MYG),
	      PVEC_REAL_MPI_TYPE,
	      UDATA_OUTDEST,
	      OUT_SENT_DOWN,
	      MPI_COMM_WORLD,
	      &ch.request[1]);
  }
  
  /// Post receive data from below (y-1)

  len = 0;

  if(DDATA_INDEST != -1) { // If sending & recieving data from a processor
    len = msg_len(ch.var_list, 0, DDATA_XSPLIT, 0, MYG);
    MPI_Irecv(ch.dmsg_recvbuff, 
	      len,
	      PVEC_REAL_MPI_TYPE,
	      DDATA_INDEST,
	      IN_SENT_UP,
	      MPI_COMM_WORLD,
	      &ch.request[2]);
  }
  if(DDATA_OUTDEST != -1) {
    inbuff = &ch.dmsg_recvbuff[len];
    MPI_Irecv(inbuff,
	      msg_len(ch.var_list, DDATA_XSPLIT, ngx, 0, MYG),
	      PVEC_REAL_MPI_TYPE,
	      DDATA_OUTDEST,
	      OUT_SENT_UP,
	      MPI_COMM_WORLD,
	      &ch.request[3]);
  }

  /// Post receive data from left (x-1)
  
  if(IDATA_DEST != -1) {
    MPI_Irecv(ch.imsg_recvbuff,
	      msg_len(ch.var_list, 0, MXG, 0, MYSUB),
	      PVEC_REAL_MPI_TYPE,
	      IDATA_DEST,
	      OUT_SENT_IN,
	      MPI_COMM_WORLD,
	      &ch.request[4]);
  }

  // Post receive data from right (x+1)

  if(ODATA_DEST != -1) {
    MPI_Irecv(ch.omsg_recvbuff,
	      msg_len(ch.var_list, 0, MXG, 0, MYSUB),
	      PVEC_REAL_MPI_TYPE,
	      ODATA_DEST,
	      IN_SENT_OUT,
	      MPI_COMM_WORLD,
	      &ch.request[5]);
  }
}

comm_handle BoutMesh::send(FieldGroup &g)
{ 
  /// Record starting wall-time
  BoutReal t = MPI_Wtime();
  
  /// Get the list of variables to send
  vector<FieldData*> var_list = g.get();
  
  /// Work out length of buffer needed
  int xlen = msg_len(var_list, 0, MXG, 0, MYSUB);
  int ylen = msg_len(var_list, 0, ngx, 0, MYG);

  /// Get a communications handle of (at least) the needed size
  CommHandle *ch = get_handle(xlen, ylen);
  ch->var_list = var_list;

  /// Post receives
  post_receive(*ch);
  
  //////////////////////////////////////////////////
  
  /// Send data going up (y+1)
  
  int len = 0;
  BoutReal *outbuff;
  
  if(UDATA_INDEST != -1) { // If there is a destination for inner x data
    len = pack_data(var_list, 0, UDATA_XSPLIT, MYSUB, MYSUB+MYG, ch->umsg_sendbuff);
    // Send the data to processor UDATA_INDEST

    if(async_send) {
      MPI_Isend(ch->umsg_sendbuff,   // Buffer to send
		len,             // Length of buffer in BoutReals
		PVEC_REAL_MPI_TYPE,  // Real variable type
		UDATA_INDEST,        // Destination processor
		IN_SENT_UP,          // Label (tag) for the message
		MPI_COMM_WORLD,
		&(ch->sendreq[0]));
    }else
      MPI_Send(ch->umsg_sendbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       UDATA_INDEST,
	       IN_SENT_UP,
	       MPI_COMM_WORLD);
  }
  if(UDATA_OUTDEST != -1) { // if destination for outer x data
    outbuff = &(ch->umsg_sendbuff[len]); // A pointer to the start of the second part
                                   // of the buffer 
    len = pack_data(var_list, UDATA_XSPLIT, ngx, MYSUB, MYSUB+MYG, outbuff);
    // Send the data to processor UDATA_OUTDEST
    if(async_send) {
      MPI_Isend(outbuff, 
		len, 
		PVEC_REAL_MPI_TYPE,
		UDATA_OUTDEST,
		OUT_SENT_UP,
		MPI_COMM_WORLD,
		&(ch->sendreq[1]));
    }else
      MPI_Send(outbuff, 
	       len, 
	       PVEC_REAL_MPI_TYPE,
	       UDATA_OUTDEST,
	       OUT_SENT_UP,
	       MPI_COMM_WORLD);
  }
    
  /// Send data going down (y-1)

  len = 0;
  if(DDATA_INDEST != -1) { // If there is a destination for inner x data
    len = pack_data(var_list, 0, DDATA_XSPLIT, MYG, 2*MYG, ch->dmsg_sendbuff);    
    // Send the data to processor DDATA_INDEST
    if(async_send) {
      MPI_Isend(ch->dmsg_sendbuff, 
		len,
		PVEC_REAL_MPI_TYPE,
		DDATA_INDEST,
		IN_SENT_DOWN,
		MPI_COMM_WORLD,
		&(ch->sendreq[2]));
    }else
      MPI_Send(ch->dmsg_sendbuff, 
	       len,
	       PVEC_REAL_MPI_TYPE,
	       DDATA_INDEST,
	       IN_SENT_DOWN,
	       MPI_COMM_WORLD);
  }
  if(DDATA_OUTDEST != -1) { // if destination for outer x data
    outbuff = &(ch->dmsg_sendbuff[len]); // A pointer to the start of the second part
			           // of the buffer
    len = pack_data(var_list, DDATA_XSPLIT, ngx, MYG, 2*MYG, outbuff);
    // Send the data to processor DDATA_OUTDEST

    if(async_send) {
      MPI_Isend(outbuff,
		len,
		PVEC_REAL_MPI_TYPE,
		DDATA_OUTDEST,
		OUT_SENT_DOWN,
		MPI_COMM_WORLD,
		&(ch->sendreq[3]));
    }else
      MPI_Send(outbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       DDATA_OUTDEST,
	       OUT_SENT_DOWN,
	       MPI_COMM_WORLD);
  }

  /// Send to the left (x-1)
  
  if(IDATA_DEST != -1) {
    len = pack_data(var_list, MXG, 2*MXG, MYG, MYG+MYSUB, ch->imsg_sendbuff);
    if(async_send) {
      MPI_Isend(ch->imsg_sendbuff,
		len,
		PVEC_REAL_MPI_TYPE,
		IDATA_DEST,
		IN_SENT_OUT,
		MPI_COMM_WORLD,
		&(ch->sendreq[4]));
    }else
      MPI_Send(ch->imsg_sendbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       IDATA_DEST,
	       IN_SENT_OUT,
	       MPI_COMM_WORLD);
  }

  /// Send to the right (x+1)

  if(ODATA_DEST != -1) {
    len = pack_data(var_list, MXSUB, MXSUB+MXG, MYG, MYG+MYSUB, ch->omsg_sendbuff);
    if(async_send) {
      MPI_Isend(ch->omsg_sendbuff,
		len,
		PVEC_REAL_MPI_TYPE,
		ODATA_DEST,
		OUT_SENT_IN,
		MPI_COMM_WORLD,
		&(ch->sendreq[5]));
    }else
      MPI_Send(ch->omsg_sendbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       ODATA_DEST,
	       OUT_SENT_IN,
	       MPI_COMM_WORLD);
  }
  
  /// Mark communication handle as in progress
  ch->in_progress = true;

  /// Add the elapsed wall-time to wtime
  wtime_comms += MPI_Wtime() - t;
  return (void*) ch;
}

int BoutMesh::wait(comm_handle handle)
{
  if(handle == NULL)
    return 1;
  
  CommHandle *ch = (CommHandle*) handle;

  if(!ch->in_progress)
    return 2;

  /// Record starting time
  BoutReal t = MPI_Wtime();
  
  ///////////// WAIT FOR DATA //////////////
  
  int ind, len;
  MPI_Status status;

  if(ch->var_list.size() == 0) {
    // Just waiting for a single MPI request
    MPI_Wait(ch->request, &status);
    free_handle(ch);

    /// Add the time elapsed to the communications wall time
    wtime_comms += MPI_Wtime() - t;
  }

  do {
    MPI_Waitany(6, ch->request, &ind, &status);
    switch(ind) {
    case 0: { // Up, inner
      unpack_data(ch->var_list, 0, UDATA_XSPLIT, MYSUB+MYG, MYSUB+2*MYG, ch->umsg_recvbuff);
      break;
    }
    case 1: { // Up, outer
      len = msg_len(ch->var_list, 0, UDATA_XSPLIT, 0, MYG);
      unpack_data(ch->var_list, UDATA_XSPLIT, ngx, MYSUB+MYG, MYSUB+2*MYG, &(ch->umsg_recvbuff[len]));
      break;
    }
    case 2: { // Down, inner
      unpack_data(ch->var_list, 0, DDATA_XSPLIT, 0, MYG, ch->dmsg_recvbuff);
      break;
    }
    case 3: { // Down, outer
      len = msg_len(ch->var_list, 0, DDATA_XSPLIT, 0, MYG);
      unpack_data(ch->var_list, DDATA_XSPLIT, ngx, 0, MYG, &(ch->dmsg_recvbuff[len]));
      break;
    }
    case 4: { // inner
      unpack_data(ch->var_list, 0, MXG, MYG, MYG+MYSUB, ch->imsg_recvbuff);
      break;
    }
    case 5: { // outer
      unpack_data(ch->var_list, MXSUB+MXG, MXSUB+2*MXG, MYG, MYG+MYSUB, ch->omsg_recvbuff);
      break;
    }
    }
    if(ind != MPI_UNDEFINED)
      ch->request[ind] = MPI_REQUEST_NULL;
    
  }while(ind != MPI_UNDEFINED);
  
  if(async_send) {
    /// Asyncronous sending: Need to check if sends have completed (frees MPI memory)
    MPI_Status status;

    if(UDATA_INDEST != -1)
      MPI_Wait(ch->sendreq, &status);
    if(UDATA_OUTDEST != -1)
      MPI_Wait(ch->sendreq+1, &status);
    if(DDATA_INDEST != -1)
      MPI_Wait(ch->sendreq+2, &status);
    if(DDATA_OUTDEST != -1)
      MPI_Wait(ch->sendreq+3, &status);
    if(IDATA_DEST != -1)
      MPI_Wait(ch->sendreq+4, &status);
    if(ODATA_DEST != -1)
      MPI_Wait(ch->sendreq+5, &status);
  }

  // TWIST-SHIFT CONDITION
  if(TwistShift && (mesh->TwistOrder == 0)) {
    int jx, jy;
    
    // Perform Twist-shift using shifting method (rather than in setStencil)
    for(std::vector<FieldData*>::iterator it = ch->var_list.begin(); it != ch->var_list.end(); it++)
      if((*it)->is3D()) {
	
	// Lower boundary

	if(TS_down_in && (DDATA_INDEST  != -1)) {
	  for(jx=0;jx<DDATA_XSPLIT;jx++)
	    for(jy=0;jy != MYG; jy++)
	      (*it)->shiftZ(jx, jy, ShiftAngle[jx]);
      
	}
	if(TS_down_out && (DDATA_OUTDEST  != -1)) {
	  for(jx=DDATA_XSPLIT;jx<ngx; jx++)
	    for(jy=0;jy != MYG; jy++)
	      (*it)->shiftZ(jx, jy, ShiftAngle[jx]);
	  
	}
	
	// Upper boundary
	
	if(TS_up_in && (UDATA_INDEST  != -1)) {
	  for(jx=0;jx<UDATA_XSPLIT; jx++)
	    for(jy=ngy-MYG;jy != ngy; jy++)
	      (*it)->shiftZ(jx, jy, -ShiftAngle[jx]);
	  
	}
	if(TS_up_out && (UDATA_OUTDEST  != -1)) {
	  for(jx=UDATA_XSPLIT;jx<ngx; jx++)
	    for(jy=ngy-MYG;jy != ngy; jy++)
	      (*it)->shiftZ(jx, jy, -ShiftAngle[jx]);
	  
	}
      }
  }

#ifdef CHECK
  // Keeping track of whether communications have been done
  for(std::vector<FieldData*>::iterator it = ch->var_list.begin(); it != ch->var_list.end(); it++)
    (*it)->doneComms();
#endif

  free_handle(ch);

  /// Add the time elapsed to the communications wall time
  wtime_comms += MPI_Wtime() - t;
  
  return 0;
}

/****************************************************************
 *                 X COMMUNICATIONS
 * 
 * Intended mainly to handle the perpendicular inversion operators
 ****************************************************************/

bool BoutMesh::firstX()
{
  return PE_XIND == 0;
}

bool BoutMesh::lastX()
{
  return PE_XIND == NXPE-1;
}

int BoutMesh::sendXOut(BoutReal *buffer, int size, int tag)
{
  if(PE_XIND == NXPE-1)
    return 1;
  
  BoutReal t = MPI_Wtime();

  MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   PROC_NUM(PE_XIND+1, PE_YIND),
	   tag,
	   MPI_COMM_WORLD);
  
  wtime_comms += MPI_Wtime() - t;

  return 0;
}

int BoutMesh::sendXIn(BoutReal *buffer, int size, int tag)
{
  if(PE_XIND == 0)
    return 1;
  
  BoutReal t = MPI_Wtime();

  MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   PROC_NUM(PE_XIND-1, PE_YIND),
	   tag,
	   MPI_COMM_WORLD);

  wtime_comms += MPI_Wtime() - t;

  return 0;
}

comm_handle BoutMesh::irecvXOut(BoutReal *buffer, int size, int tag)
{
  if(PE_XIND == NXPE-1)
    return NULL;

  BoutReal t = MPI_Wtime();
  
  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0,0);

  MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    PROC_NUM(PE_XIND+1, PE_YIND),
	    tag,
	    MPI_COMM_WORLD,
	    ch->request);
  
  ch->in_progress = true;

  wtime_comms += MPI_Wtime() - t;

  return (void*) ch;
}

comm_handle BoutMesh::irecvXIn(BoutReal *buffer, int size, int tag)
{
  if(PE_XIND == 0)
    return NULL;
  
  BoutReal t = MPI_Wtime();

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0,0);

  MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    PROC_NUM(PE_XIND-1, PE_YIND),
	    tag,
	    MPI_COMM_WORLD,
	    ch->request);
  
  ch->in_progress = true;

  wtime_comms += MPI_Wtime() - t;

  return (void*) ch;
}

/****************************************************************
 *                 GRID INDEX ROUTINES
 * 
 * These routines translate between local and global coordinates
 * This so that more flexible indexing (e.g. not equal on processors)
 * can be implemented later (maybe)
 ****************************************************************/

/// Returns the processor number, given X and Y processor indices.
/*!
 * If out of range returns -1 (no processor)
 */

int BoutMesh::PROC_NUM(int xind, int yind)
{
  if((xind >= NXPE) || (xind < 0))
    return -1;
  if((yind >= NYPE) || (yind < 0))
    return -1;
  
  return yind * NXPE + xind;
}

/// Returns true if the given grid-point coordinates are in this processor
bool BoutMesh::IS_MYPROC(int xind, int yind)
{
  return ((xind / MXSUB) == PE_XIND) && ((yind / MYSUB) == PE_YIND);
}

/// Returns the global X index given a local index
int BoutMesh::XGLOBAL(int xloc)
{
  return xloc + PE_XIND * MXSUB;
}

/// Returns a local X index given a global index
int BoutMesh::XLOCAL(int xglo)
{
  return xglo - PE_XIND * MXSUB;
}

/// Returns the global Y index given a local index
int BoutMesh::YGLOBAL(int yloc)
{
  return yloc + PE_YIND*MYSUB - MYG;
}

/// Global Y index given local index and processor
int BoutMesh::YGLOBAL(int yloc, int yproc)
{
  return yloc + yproc*MYSUB - MYG;
}

/// Returns a local Y index given a global index
int BoutMesh::YLOCAL(int yglo)
{
  return yglo - PE_YIND*MYSUB + MYG;
}

int BoutMesh::YLOCAL(int yglo, int yproc)
{
  return yglo - yproc*MYSUB + MYG;
}

/// Return the Y processor number given a global Y index
int BoutMesh::YPROC(int yind)
{
  return yind / MYSUB;
}

/// Return the X processor number given a global X index
int BoutMesh::XPROC(int xind)
{
  return (xind >= MXG) ? (xind - MXG) / MXSUB : 0;
}

/****************************************************************
 *                       CONNECTIONS
 ****************************************************************/

/// Connection initialisation: Set processors in a simple 2D grid
void BoutMesh::default_connections()
{
  DDATA_XSPLIT = UDATA_XSPLIT = 0;  // everything by default outside (arb. choice)
  DDATA_INDEST = UDATA_INDEST = -1; // since nothing inside

  DDATA_OUTDEST = PROC_NUM(PE_XIND, PE_YIND - 1);
  UDATA_OUTDEST = PROC_NUM(PE_XIND, PE_YIND + 1);
  
  IDATA_DEST = PROC_NUM(PE_XIND - 1, PE_YIND);
  ODATA_DEST = PROC_NUM(PE_XIND + 1, PE_YIND);
  
  TS_up_in = TS_up_out = TS_down_in = TS_down_out = false; // No twist-shifts

  /// Check if X is periodic
  bool xperiodic;
  options.setSection("");
  options.get("xperiodic", xperiodic, false);
  if(xperiodic) {
    if(PE_XIND == (NXPE-1))
      ODATA_DEST = PROC_NUM(0, PE_YIND);
    
    if(PE_XIND == 0)
      IDATA_DEST = PROC_NUM(NXPE-1, PE_YIND);
  }
}

/// Add a topology connection
/*!
 * Set ypos1 and ypos2 to be neighbours in the range xge <= x < xlt.
 * Optional argument ts sets whether to use twist-shift condition
 */
void BoutMesh::set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts)
{
  int ype1, ype2; // the two Y processor indices
  int ypeup, ypedown;
  int yind1, yind2;

  if(xlt <= xge)
    return;

  if((ypos1 < 0) || (ypos1 >= MY)) {
    output.write("WARNING adding connection: poloidal index %d out of range\n", ypos1);
    return;
  }
  if((ypos2 < 0) || (ypos2 >= MY)) {
    output.write("WARNING adding connection: poloidal index %d out of range\n", ypos2);
    return;
  }
  
  ype1 = YPROC(ypos1);
  ype2 = YPROC(ypos2);

  /* y index within processors */
  yind1 = YLOCAL(ypos1, ype1);
  yind2 = YLOCAL(ypos2, ype2);

  /* Check which boundary the connection is on */
  if((yind1 == MYG) && (yind2 == MYSUB+MYG-1)) {
    ypeup = ype2; /* processor sending data up (+ve y) */
    ypedown = ype1; /* processor sending data down (-ve y) */
  }else if((yind2 == MYG) && (yind1 == MYSUB+MYG-1)) {
    ypeup = ype1;
    ypedown = ype2;
  }else {
    output.write("ERROR adding connection: y index %d or %d not on processor boundary\n", ypos1, ypos2);
    exit(1);
  }

  /* check the x ranges are possible */
  if((xge != 0) && (xlt != MX)) {
    output.write("ERROR adding connection(%d,%d,%d,%d): can only divide X domain in 2\n",
		 ypos1, ypos2, xge, xlt);
    exit(1);
  }

  output.write("Connection between top of Y processor %d and bottom of %d in range %d <= x < %d\n",
	       ypeup, ypedown, xge, xlt);

  // Convert X coordinates into local indices

  xge = XLOCAL(xge);
  xlt = XLOCAL(xlt);
  
  if(( xge >= ngx ) || ( xlt <= 0 )) {
    return; // Not in this x domain
  }
  
  if(xge < 0)   xge = 0;
  if(xlt > ngx) xlt = ngx;

  if(MYPE == PROC_NUM(PE_XIND, ypeup)) { /* PROCESSOR SENDING +VE Y */
    /* Set the branch cut x position */
    if(xge <= MXG) {
      /* Connect on the inside */
      UDATA_XSPLIT = xlt;
      UDATA_INDEST = PROC_NUM(PE_XIND, ypedown);
      if(UDATA_XSPLIT == ngx)
	UDATA_OUTDEST = -1;

      TS_up_in = ts; // Twist-shift
      
      output.write("=> This processor sending in up\n");
    }else {
      /* Connect on the outside */
      UDATA_XSPLIT = xge;
      UDATA_OUTDEST = PROC_NUM(PE_XIND, ypedown);
      if(UDATA_XSPLIT == 0)
	UDATA_INDEST = -1;

      TS_up_out = ts;
      output.write("=> This processor sending out up\n");
    }
  }
  
  if(MYPE == PROC_NUM(PE_XIND, ypedown)) { /* PROCESSOR SENDING -VE Y */
    /* Set the branch cut x position */
    if(xge <= MXG) {
      /* Connect on the inside */
      DDATA_XSPLIT = xlt;
      DDATA_INDEST = PROC_NUM(PE_XIND, ypeup);
      if(DDATA_XSPLIT == ngx)
	DDATA_OUTDEST = -1;

      TS_down_in = ts;

      output.write("=> This processor sending in down\n");
    }else {
      /* Connect on the outside */
      DDATA_XSPLIT = xge;
      DDATA_OUTDEST = PROC_NUM(PE_XIND, ypeup);
      if(DDATA_XSPLIT == 0)
	DDATA_INDEST = -1;

      TS_down_out = ts;

      output.write("=> This processor sending out down\n");
    }
  }
}

/****************************************************************
 *                MAIN TOPOLOGY FUNCTION
 ****************************************************************/

void BoutMesh::topology()
{ 
  // Perform checks common to all topologies

  if (NPES != NXPE*NYPE) {
    output.write("\tTopology error: npes=%d is not equal to NXPE*NYPE=%d\n",
		 NPES,NXPE*NYPE);
    exit(1);
  }
  if(MYSUB * NYPE != MY) {
    output.write("\tTopology error: MYSUB[%d] * NYPE[%d] != MY[%d]\n",MYSUB,NYPE,MY);
    exit(1);
  }
  if(MXSUB * NXPE != MX) {
    output.write("\tTopology error: MXSUB[%d] * NXPE[%d] != MX[%d]\n",MXSUB,NXPE,MX);
    exit(1);
  }

  if((NXPE > 1) && (MXSUB < MXG)) {
    output.write("\tERROR: Grid X size must be >= guard cell size\n");
    exit(1);
  }
  if(MYSUB < MYG) {
    output.write("\tERROR: Grid Y size must be >= guard cell size\n");
    exit(1);
  }
  
  if(jyseps2_1 == jyseps1_2) {
    /********* SINGLE NULL OPERATION *************/
    output.write("\tEQUILIBRIUM IS SINGLE NULL (SND) \n");

    /* Set separatrices - x location all the same */
    ixseps_inner = ixseps_outer = ixseps_upper = ixseps_lower = ixseps1;

    default_connections();
    set_connection(jyseps1_1+1, jyseps2_2, 0, ixseps1, true); // Twist-shift this connection
    set_connection(jyseps1_1, jyseps2_2+1, 0, ixseps1); // No twist-shift in PF region
    
  }else {
    /*************** DOUBLE NULL OPERATION *******************/
    /* UPPER LEGS: Do not have to be the same length as each
       other or lower legs, but do have to have an integer number
       of processors */
    if((ny_inner-jyseps2_1-1) % MYSUB != 0) {
      output.write("\tTopology error: Upper inner leg does not have integer number of processors\n");
      exit(1);
    }
    if((jyseps1_2-ny_inner+1) % MYSUB != 0) {
      output.write("\tTopology error: Upper outer leg does not have integer number of processors\n");
    }

    if(ixseps1 == ixseps2) {
      /*************** CONNECTED (balanced) DOUBLE NULL ******************/
      output.write("\tEQUILIBRIUM IS CONNECTED DOUBLE NULL (CDND)\n");
      /* all separatrix indices the same */
      ixseps_inner = ixseps_outer = ixseps_upper = ixseps_lower = ixseps1;
        
    }else if(ixseps1 < ixseps2){
      /*************** LOWER DOUBLE NULL **********************/
      output.write("\tEQUILIBRIUM IS LOWER DOUBLE NULL (LDND)\n");
      ixseps_inner = ixseps_lower = ixseps1;
      ixseps_outer = ixseps_upper = ixseps2;
    }else {
      /*************** UPPER DOUBLE NULL **********************/
      output.write("\tEQUILIBRIUM IS UPPER DOUBLE NULL (UDND)\n");
      ixseps_inner = ixseps_upper = ixseps2;
      ixseps_outer = ixseps_lower = ixseps1;
    }
    
    /* Following code works for any Double Null */

    /********* DND CONNECTIONS **********/
    default_connections();
    /* Lower x-point */
    set_connection(jyseps1_1+1, jyseps2_2  , 0, ixseps_lower, ixseps1 <= ixseps2); /* Core */
    set_connection(jyseps1_1  , jyseps2_2+1, 0, ixseps_lower); /* PF   */
    /* Upper x-point */
    set_connection(jyseps2_1  , jyseps1_2+1, 0, ixseps_upper, ixseps1 > ixseps2); /* Core */
    set_connection(jyseps2_1+1, jyseps1_2  , 0, ixseps_upper); /* PF   */
  }

  MYPE_IN_CORE = 0; // processor not in core
  if( (ixseps_inner > 0) && ( ((PE_YIND*MYSUB > jyseps1_1) && (PE_YIND*MYSUB <= jyseps2_1)) || ((PE_YIND*MYSUB > jyseps1_2) && (PE_YIND*MYSUB <= jyseps2_2)) ) ) {
    MYPE_IN_CORE = 1; /* processor is in the core */
  }

  // Print out settings
  output.write("\tMYPE_IN_CORE = %d\n", MYPE_IN_CORE);
  output.write("\tDXS = %d, DIN = %d. DOUT = %d\n", 
	       DDATA_XSPLIT, DDATA_INDEST, DDATA_OUTDEST);
  output.write("\tUXS = %d, UIN = %d. UOUT = %d\n", 
	       UDATA_XSPLIT, UDATA_INDEST, UDATA_OUTDEST);
  output.write("\tXIN = %d, XOUT = %d\n",
	       IDATA_DEST, ODATA_DEST);

  output.write("\tTwist-shift: ");
  if(TS_down_in)
    output.write("DI ");
  if(TS_down_out)
    output.write("DO ");
  if(TS_up_in)
    output.write("UI ");
  if(TS_up_out)
    output.write("UO ");
  output.write("\n");
}

/****************************************************************
 *                     Communication handles
 ****************************************************************/

BoutMesh::CommHandle* BoutMesh::get_handle(int xlen, int ylen)
{
  if(comm_list.empty()) {
    //Allocate a new CommHandle
    
    CommHandle* ch = new CommHandle;
    for(int i=0;i<6;i++)
      ch->request[i] = MPI_REQUEST_NULL;
    
    if(ylen > 0) {
      ch->umsg_sendbuff = new BoutReal[ylen];
      ch->dmsg_sendbuff = new BoutReal[ylen];
      ch->umsg_recvbuff = new BoutReal[ylen];
      ch->dmsg_recvbuff = new BoutReal[ylen];
    }
    
    if(xlen > 0) {
      ch->imsg_sendbuff = new BoutReal[xlen];
      ch->omsg_sendbuff = new BoutReal[xlen];
      ch->imsg_recvbuff = new BoutReal[xlen];
      ch->omsg_recvbuff = new BoutReal[xlen];
    }
    
    ch->xbufflen = xlen;
    ch->ybufflen = ylen;
    
    ch->in_progress = false;
    
    return ch;
  }
  
  // Pop first pointer off the list
  CommHandle* ch = comm_list.front();
  comm_list.pop_front();

  // Check that the buffers are big enough (NOTE: Could search list for bigger buffers)
  if(ch->ybufflen < ylen) {
    if(ch->ybufflen > 0) {
      delete[] ch->umsg_sendbuff;
      delete[] ch->umsg_recvbuff;
      delete[] ch->dmsg_sendbuff;
      delete[] ch->dmsg_recvbuff;
    }
    
    ch->umsg_sendbuff = new BoutReal[ylen];
    ch->dmsg_sendbuff = new BoutReal[ylen];
    ch->umsg_recvbuff = new BoutReal[ylen];
    ch->dmsg_recvbuff = new BoutReal[ylen];
    
    ch->ybufflen = ylen;
  }
  if(ch->xbufflen < xlen) {
    if(ch->ybufflen > 0) {
      delete[] ch->imsg_sendbuff;
      delete[] ch->imsg_recvbuff;
      delete[] ch->omsg_sendbuff;
      delete[] ch->omsg_recvbuff;
    }
    
    ch->imsg_sendbuff = new BoutReal[xlen];
    ch->omsg_sendbuff = new BoutReal[xlen];
    ch->imsg_recvbuff = new BoutReal[xlen];
    ch->omsg_recvbuff = new BoutReal[xlen];
    
    ch->xbufflen = xlen;
  }
  
  ch->in_progress = false;
  
  return ch;
}

void BoutMesh::free_handle(CommHandle *h)
{
  h->var_list.clear();
  comm_list.push_front(h);
}

void BoutMesh::clear_handles()
{
  while(!comm_list.empty()) {
    CommHandle *ch = comm_list.front();
    delete[] ch->umsg_sendbuff;
    delete[] ch->dmsg_sendbuff;
    delete[] ch->imsg_sendbuff;
    delete[] ch->omsg_sendbuff;
    
    delete[] ch->umsg_recvbuff;
    delete[] ch->dmsg_recvbuff;
    delete[] ch->imsg_recvbuff;
    delete[] ch->omsg_recvbuff;
    
    delete ch;
    
    comm_list.pop_front();
  }
}

/****************************************************************
 *                   Communication utilities
 ****************************************************************/

int BoutMesh::pack_data(vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt, BoutReal *buffer)
{
  int jx, jy, jz;
  int len = 0;
  std::vector<FieldData*>::iterator it;
  
  for(jx=xge; jx != xlt; jx++) {
    
    /// Loop over variables
    for(it = var_list.begin(); it != var_list.end(); it++) {
      if((*it)->is3D()) {
	// 3D variable
	
	for(jy=yge;jy < ylt;jy++)
	  for(jz=0;jz < mesh->ngz-1;jz++)
	    len += (*it)->getData(jx,jy,jz,buffer+len);
	
      }else {
	// 2D variable
	for(jy=yge;jy < ylt;jy++)
	  len += (*it)->getData(jx,jy,0,buffer+len);
      }
    }
    
  }
  
  return(len);
}

int BoutMesh::unpack_data(vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt, BoutReal *buffer)
{
  int jx, jy, jz;
  int len = 0;
  std::vector<FieldData*>::iterator it;

  for(jx=xge; jx != xlt; jx++) {

    /// Loop over variables
    for(it = var_list.begin(); it != var_list.end(); it++) {
      if((*it)->is3D()) {
	// 3D variable
   
	for(jy=yge;jy < ylt;jy++)
	  for(jz=0;jz < mesh->ngz-1;jz++) {
	    len += (*it)->setData(jx,jy,jz,buffer+len);
	  }
	
      }else {
	// 2D variable
	for(jy=yge;jy < ylt;jy++)
	  len += (*it)->setData(jx,jy,0,buffer+len);
      }
    }
    
  }
  
  return(len);
}

int BoutMesh::msg_len(vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt)
{
  int len = 0;

  /// Loop over variables
  for(std::vector<FieldData*>::iterator it = var_list.begin(); it != var_list.end(); it++) {
    if((*it)->is3D()) {
      len += (xlt - xge) * (ylt - yge) * (mesh->ngz-1) * (*it)->BoutRealSize();
    }else
      len += (xlt - xge) * (ylt - yge) * (*it)->BoutRealSize();
  }
  
  return len;
}

/****************************************************************
 *                   Private data reading routines
 ****************************************************************/

/// Reads in a portion of the X-Y domain
int BoutMesh::readgrid_3dvar(GridDataSource *s, const char *name, 
	                     int yread, int ydest, int ysize, 
                             int xge, int xlt, BoutReal ***var)
{
  /// Check the arguments make sense
  if((yread < 0) || (ydest < 0) || (ysize < 0) || (xge < 0) || (xlt < 0))
    return 1;
  
  /// Check the size of the data
  vector<int> size = s->getSize(name);
  
  if(size.size() != 3) {
    output.write("\tWARNING: Number of dimensions of %s incorrect\n", name);
    return 1;
  }
  
  if((size[0] != nx) || (size[1] != ny)) {
    output.write("\tWARNING: X or Y size of %s incorrect\n", name);
    return 1;
  }

  if((size[2] & 1) != 1) {
    output.write("\tWARNING: Z size of %s should be odd\n", name);
    return 1;
  }

  int maxmode = (size[2] - 1)/2; ///< Maximum mode-number n

  int ncz = mesh->ngz-1;

  // Print out which modes are going to be read in
  if(zperiod > maxmode) {
    // Domain is too small: Only DC
    output.write(" => Only reading n = 0 component\n");
  }else {
    // Get maximum mode in the input which is a multiple of zperiod
    int mm = ((int) (maxmode/zperiod))*zperiod;
    if( (ncz/2)*zperiod < mm )
      mm = (ncz/2)*zperiod; // Limited by Z resolution
    
    if(mm == zperiod) {
      output.write(" => Reading n = 0, %d\n", zperiod);
    }else
      output.write(" => Reading n = 0, %d ... %d\n", zperiod, mm);
  }

  /// Data for FFT. Only positive frequencies
  dcomplex* fdata = new dcomplex[ncz/2 + 1];
  BoutReal* zdata = new BoutReal[size[2]];

  for(int jx=xge;jx<xlt;jx++) {
    // Set the global X index
    
    for(int jy=0; jy < ysize; jy++) {
      /// Read data
      
      int yind = yread + jy; // Global location to read from
      
      s->setOrigin(XGLOBAL(jx), yind);
      if(!s->fetch(zdata, name, 1, 1, size[2]))
	return 1;
      
      /// Load into dcomplex array
      
      fdata[0] = zdata[0]; // DC component

      for(int i=1;i<=ncz/2;i++) {
	int modenr = i*zperiod; // Z mode number
	
	if(modenr <= maxmode) {
	  // Have data for this mode
	  fdata[i] = dcomplex(zdata[modenr*2 - 1], zdata[modenr*2]);
	}else {
	  fdata[i] = 0.0;
	}
      }
      
      // Inverse FFT, shifting in the z direction
      for(int jz=0;jz<=ncz/2;jz++) {
	BoutReal kwave;
	
	kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
      
	// Multiply by EXP(ik*zoffset)
	fdata[jz] *= dcomplex(cos(kwave*zShift[jx][jy]) , sin(kwave*zShift[jx][jy]));
      }
      
      irfft(fdata, ncz, var[jx][ydest+jy]);
    }
  }

  s->setOrigin();

  // free data
  delete[] zdata;
  delete[] fdata;
  
  return 0;
}

/// Copies a section of a 3D variable
void BoutMesh::cpy_3d_data(int yfrom, int yto, int xge, int xlt, BoutReal ***var)
{
  int i, k;
  for(i=xge;i!=xlt;i++)
    for(k=0;k<ngz;k++)
      var[i][yto][k] = var[i][yfrom][k];
}

/// Helper routine for reading in grid data
/*!
  Reads in a single 2D variable varname[xge:xlt-1][yread:yread+ysize-1]
  and puts it into var[xge:xlt-1][ydest:ydesy+ysize-1]
  
  July 2008: Adapted to take into account X offsets
*/
int BoutMesh::readgrid_2dvar(GridDataSource *s, const char *varname, 
                             int yread, int ydest, int ysize, 
                             int xge, int xlt, BoutReal **var)
{
  for(int i=xge;i!=xlt;i++) { // go through all the x indices 
    // Set the indices to read in this x position 
    s->setOrigin(XGLOBAL(i), yread);
    // Read in the block of data for this x value (C ordering)
    if(!s->fetch(&(var[i][ydest]), varname, 1, ysize))
      return 1;
  }
  
  s->setOrigin();
  
  return 0;
}

void BoutMesh::cpy_2d_data(int yfrom, int yto, int xge, int xlt, BoutReal **var)
{
  int i;
  for(i=xge;i!=xlt;i++)
    var[i][yto] = var[i][yfrom];
}


/****************************************************************
 *                 SURFACE ITERATION
 ****************************************************************/

SurfaceIter* BoutMesh::iterateSurfaces()
{
  //return new BoutSurfaceIter(this);
  return (SurfaceIter*) NULL;
}

bool BoutMesh::surfaceClosed(int jx, BoutReal &ts)
{
  ts = 0.;
  if(jx < ixseps_inner) {
    if(TwistShift)
      ts = ShiftAngle[jx];
    return true;
  }
  return false;
}

// Define MPI operation to sum 2D fields over y.
// NB: Don't sum in y boundary regions
void ysum_op(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
    BoutReal *rin = (BoutReal*) invec;
    BoutReal *rinout = (BoutReal*) inoutvec;
    for(int x=0;x<mesh->ngx;x++) {
	BoutReal val = 0.;
	// Sum values
	for(int y=mesh->ystart;y<=mesh->yend;y++) {
	    val += rin[x*mesh->ngy + y] + rinout[x*mesh->ngy + y];
	}
	// Put into output (spread over y)
	val /= mesh->yend - mesh->ystart + 1;
	for(int y=0;y<mesh->ngy;y++)
	    rinout[x*mesh->ngy + y] = val;
    }
}

const Field2D BoutMesh::averageY(const Field2D &f)
{
  static MPI_Op op;
  static bool opdefined = false;

#ifdef CHECK
  msg_stack.push("averageY(Field2D)");
#endif

  if(!opdefined) {
    MPI_Op_create(ysum_op, 1, &op);
    opdefined = true;
  }

  Field2D result;
  result.allocate();
  
  BoutReal **fd, **rd;
  fd = f.getData();
  rd = result.getData();
  
  MPI_Allreduce(*fd, *rd, mesh->ngx*mesh->ngy, MPI_DOUBLE, op, comm_inner);
  
  result /= (BoutReal) NYPE;

#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}

/*
BoutSurfaceIter::BoutSurfaceIter(BoutMesh* mi)
{
  m = mi;
}

int BoutSurfaceIter::ysize()
{
  int xglobal = m->XGLOBAL(xpos);
  int yglobal = m->YGLOBAL(MYG);
  
  if((xglobal < m->ixseps_lower) && ((yglobal <= m->jyseps1_1) || (yglobal > m->jyseps2_2))) {
    // Lower PF region
    return (m->jyseps1_1 + 1) + (ny - m->jyseps2_2);
    
  }else if((xglobal < m->ixseps_upper) && (yglobal > m->jyseps2_1) && (yglobal >= m->jyseps1_2)) {
    // Upper PF region
    return m->jyseps1_2 - m->jyseps2_1;
    
  }else if(xglobal < m->ixseps_inner) {
    // Core
    return (m->jyseps2_1 - m->jyseps1_1) + (m->jyseps2_2 - m->jyseps1_2);
    
  }else if(m->jyseps2_1 == m->jyseps1_2) {
    // Single null, so in the SOL
    return m->ny;
    
  }else if((xglobal >= m->ixseps_inner) && (xglobal < m->ixseps_outer)) {
    // Intermediate SOL in DND
    
    if(m->ixseps_lower < m->ixseps_upper) {
      // Connects to lower divertor
      return (m->jyseps2_1 + 1) + (ny - m->jyseps1_2);
    }else {
      // Connects to upper divertor
      return m->jyseps2_2 - m->jyseps1_1;
    }
  }else if(yglobal < m->ny_inner) {
    // Inner SOL
    return m->ny_inner;
  }
  // Outer SOL
  return m->ny - m->ny_inner;
}

bool BoutSurfaceIter::closed(BoutReal &ts)
{
  int xglobal = m->XGLOBAL(xpos);
  ts = 0.;
  if(TwistShift) {
    ts = ShiftAngle[xpos];
  }
  return m->MYPE_IN_CORE && (xglobal < m->ixseps_inner);
}

void BoutSurfaceIter::first()
{
  xpos = m->PE_XIND;
  if(xpos > m->ngx-1)
    xpos = -1; // Nothing to do
}

void BoutSurfaceIter::next()
{
  if(xpos < 0)
    return;
  
  xpos += NXPE;
  if(xpos > m->ngx-1)
    xpos = -1; // Nothing to do 
}

bool BoutSurfaceIter::isDone()
{
  return xpos < 0;
}

int BoutSurfaceIter::gather(const Field2D &f, BoutReal *data)
{
  
}

int BoutSurfaceIter::gather(const Field3D &f, BoutReal **data)
{
  
}

int BoutSurfaceIter::scatter(BoutReal *data, Field2D &f)
{
  
}

int BoutSurfaceIter::scatter(BoutReal **data, Field3D &f)
{
  
}

*/

/****************************************************************
 *                 Range iteration
 ****************************************************************/

RangeIter* BoutMesh::iterateBndryLowerY()
{
  int xs = 0;
  int xe = ngx-1;
  if((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart))
    xs = DDATA_XSPLIT;
  if((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend+1))
    xe = DDATA_XSPLIT-1;

  if(xs < xstart)
    xs = xstart;
  if(xe > xend)
    xe = xend;

  return new BoutRangeIter(xs, xe);
}

RangeIter* BoutMesh::iterateBndryUpperY()
{
  int xs = 0;
  int xe = ngx-1;
  if((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart))
    xs = UDATA_XSPLIT;
  if((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend+1))
    xe = UDATA_XSPLIT-1;

  if(xs < xstart)
    xs = xstart;
  if(xe > xend)
    xe = xend;

  return new BoutRangeIter(xs, xe);
}


vector<BoundaryRegion*> BoutMesh::getBoundaries()
{
  return boundary;
}

BoutRangeIter::BoutRangeIter(int start, int end)
{
  s = start;
  e = end;
}

void BoutRangeIter::first()
{
  ind = s;
}

void BoutRangeIter::next()
{
  ind++;
}

bool BoutRangeIter::isDone()
{
  return ind > e;
}

BoutReal BoutMesh::GlobalX(int jx)
{
  return ((BoutReal) XGLOBAL(jx)) / ((BoutReal) MX);
}

BoutReal BoutMesh::GlobalY(int jy)
{
  int ly = YGLOBAL(jy); // global poloidal index across subdomains
  int nycore = (jyseps1_2 - jyseps1_1) + (jyseps2_2 - jyseps2_1);

  if(MYPE_IN_CORE) {
    // Turn ly into an index over the core cells only
    if(ly < jyseps1_2) {
      ly -= jyseps1_1+1;
    }else
      ly -= jyseps1_1+1 + (jyseps2_1 - jyseps1_2);
  }else {
    // Not in core. Need to get the last "core" value
    if(ly <= jyseps1_1) {
      // Inner lower leg
      ly = 0;
    }else if(ly > jyseps2_2) {
      // Outer lower leg
      ly = nycore-1;
    }
  }
  return ((BoutReal) ly) / ((BoutReal) nycore);
}

void BoutMesh::outputVars(Datafile &file)
{
  file.add(MXSUB, "MXSUB", 0);
  file.add(MYSUB, "MYSUB", 0);
  file.add(MXG,   "MXG",   0);
  file.add(MYG,   "MYG",   0);
  file.add(ngz,   "MZ",    0);
  file.add(NXPE,  "NXPE",  0);
  file.add(NYPE,  "NYPE",  0);
  file.add(ZMAX,  "ZMAX",  0);
  file.add(ZMIN,  "ZMIN",  0);
}
