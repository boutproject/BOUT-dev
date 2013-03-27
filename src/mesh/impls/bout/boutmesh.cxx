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

#include "boutmesh.hxx"

#include <utils.hxx>
#include <fft.hxx>
#include <derivs.hxx>
#include <boutcomm.hxx>
#include <dcomplex.hxx>
#include <options.hxx>
#include <boutexception.hxx>
#include <output.hxx>
#include <bout/sys/timer.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>

#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

//#define COMMDEBUG 1   // Uncomment to print communications debugging information

BoutMesh::BoutMesh(GridDataSource *s, Options *options) : Mesh(s) {
  if(options == NULL)
    options = Options::getRoot()->getSection("mesh");
  
  OPTION(options, symmetricGlobalX,  false);

  comm_x = MPI_COMM_NULL;
  comm_inner = MPI_COMM_NULL;
  comm_middle = MPI_COMM_NULL;
  comm_outer = MPI_COMM_NULL;
}

BoutMesh::~BoutMesh() {
  // Delete the communication handles
  clear_handles();
  
  // Delete the boundary regions
  for(vector<BoundaryRegion*>::iterator it = boundary.begin(); it != boundary.end(); it++)
    delete (*it);

  delete[] ShiftAngle;
  
  
  if(comm_x != MPI_COMM_NULL)
    MPI_Comm_free(&comm_x);
  if(comm_inner != MPI_COMM_NULL)
    MPI_Comm_free(&comm_inner);
  //if(comm_middle != MPI_COMM_NULL)
  //  MPI_Comm_free(&comm_middle); // Already freed
  
  if(comm_outer != MPI_COMM_NULL)
    MPI_Comm_free(&comm_outer);
  
}

int BoutMesh::load() {
#ifdef CHECK
  int msg = msg_stack.push("BoutMesh::load()");
#endif
  
  output << "Loading mesh" << endl;
  
  // Use root level options
  Options *options = Options::getRoot();

  //////////////
  // Number of processors
  
  MPI_Comm_size(BoutComm::get(), &NPES);
  MPI_Comm_rank(BoutComm::get(), &MYPE);

  //////////////
  // Grid sizes
  
  if(Mesh::get(nx, "nx"))
    throw BoutException("Mesh must contain nx");
  
  if(Mesh::get(ny, "ny"))
    throw BoutException("Mesh must contain ny");
  
  int MZ;
  OPTION(options, MZ,           65);
  if(!is_pow2(MZ-1)) {
    if(is_pow2(MZ)) {
      MZ++;
      output.write("WARNING: Number of toroidal points increased to %d\n", MZ);
    }else {
      throw BoutException("Error: Number of toroidal points must be 2^n + 1");
    }
  }

  output << "\tGrid size: " << nx << " x " << ny << " x " << MZ << endl;

  // Set global grid sizes
  GlobalNx = nx;
  GlobalNy = ny + 4;
  GlobalNz = MZ;

  options->get("MXG", MXG, 2);
  options->get("MYG", MYG, 2);
  
  // separatrix location
  if(Mesh::get(ixseps1, "ixseps1")) {
    ixseps1 = GlobalNx;
    output.write("\tWARNING: Separatrix location 'ixseps1' not found. Setting to %d\n", ixseps1);
  }
  if(Mesh::get(ixseps2, "ixseps2")) {
    ixseps2 = GlobalNx;
    output.write("\tWARNING: Separatrix location 'ixseps2' not found. Setting to %d\n", ixseps2);
  }
  if(Mesh::get(jyseps1_1,"jyseps1_1")) {
    jyseps1_1 = -1;
    output.write("\tWARNING: Branch-cut 'jyseps1_1' not found. Setting to %d\n", jyseps1_1);
  }
  if(Mesh::get(jyseps1_2,"jyseps1_2")) {
    jyseps1_2 = ny/2;
    output.write("\tWARNING: Branch-cut 'jyseps1_2' not found. Setting to %d\n", jyseps1_2);
  }
  if(Mesh::get(jyseps2_1,"jyseps2_1")) {
    jyseps2_1 = jyseps1_2;
    output.write("\tWARNING: Branch-cut 'jyseps2_1' not found. Setting to %d\n", jyseps2_1);
  }
  if(Mesh::get(jyseps2_2,"jyseps2_2")) {
    jyseps2_2 = ny-1;
    output.write("\tWARNING: Branch-cut 'jyseps2_2' not found. Setting to %d\n", jyseps2_2);
  }

  if(Mesh::get(ny_inner,"ny_inner")) {
    ny_inner = jyseps2_1;
    output.write("\tWARNING: Number of inner y points 'ny_inner' not found. Setting to %d\n", ny_inner);
  }

  /// Check inputs
  if(jyseps1_1 < -1) {
    output.write("\tWARNING: jyseps1_1 (%d) must be >= -1. Setting to -1\n", jyseps1_1);
    jyseps1_1 = -1;
  }
  
  if(jyseps2_1 <= jyseps1_1) {
    output.write("\tWARNING: jyseps2_1 (%d) must be > jyseps1_1 (%d). Setting to %d\n",
                 jyseps2_1, jyseps1_1, jyseps1_1+1);
    jyseps2_1 = jyseps1_1 + 1;
  }
  if(jyseps1_2 < jyseps2_1) {
    output.write("\tWARNING: jyseps1_2 (%d) must be >= jyseps2_1 (%d). Setting to %d\n",
                 jyseps1_2, jyseps2_1, jyseps2_1);
    jyseps1_2 = jyseps2_1;
  }
  if(jyseps2_2 >= ny) {
    output.write("\tWARNING: jyseps2_2 (%d) must be < ny (%d). Setting to %d\n",
                 jyseps2_2, ny, ny-1);
    jyseps2_2 = ny - 1;
  }

  if(options->isSet("NXPE")) { // Specified NXPE
    options->get("NXPE", NXPE, 1); // Decomposition in the radial direction
    if((NPES % NXPE) != 0) {
      throw BoutException("Number of processors (%d) not divisible by NPs in x direction (%d)\n",
                          NPES, NXPE);
    }
    
    NYPE = NPES / NXPE;
  }else {
    // Choose NXPE
    
    MX = nx - 2*MXG;
    
    NXPE = -1; // Best option 
    
    BoutReal ideal = sqrt(MX * NPES / ((double) ny)); // Results in square domains

    output.write("Finding value for NXPE\n");

    for(int i=1; i<= NPES; i++) { // Loop over all possibilities
      //output.write("Testing %d: %d, %d, %d, %d, %d\n",
      //             i, NPES % i, MX % i, MX / i, ny % (NPES/i), ny / (NPES/i));
      if( (NPES % i == 0) &&      // Processors divide equally
          (MX % i == 0) &&        // Mesh in X divides equally
    //      (MX / i >= MXG) &&      // Resulting mesh is large enough
          (ny % (NPES/i) == 0) ) { // Mesh in Y divides equally
        
        output.write("\tCandidate value: %d\n", i);
        
        int nyp = NPES/i;
        int ysub = ny / nyp;
        
        // Check size of Y mesh
        if(ysub < MYG) {
          output.write("\t -> ny/NYPE (%d/%d = %d) must be >= MYG (%d)\n", ny, nyp, ysub, MYG);
          continue;
        }
        // Check branch cuts
        if( (jyseps1_1+1) % ysub != 0 ) {
          output.write("\t -> Leg region jyseps1_1+1 (%d) must be a multiple of MYSUB (%d)\n", jyseps1_1+1, ysub);
          continue;
        }
        
        if(jyseps2_1 != jyseps1_2) {
          // Double Null
          
          if( (jyseps2_1-jyseps1_1) % ysub != 0 ) {
            output.write("\t -> Core region jyseps2_1-jyseps1_1 (%d-%d = %d) must be a multiple of MYSUB (%d)\n", 
                         jyseps2_1, jyseps1_1, jyseps2_1-jyseps1_1, ysub);
            continue;
          }
          
          if( (jyseps2_2 - jyseps1_2) % ysub != 0 ) {
            output.write("\t -> Core region jyseps2_2-jyseps1_2 (%d-%d = %d) must be a multiple of MYSUB (%d)\n", 
                         jyseps2_2, jyseps1_2, jyseps2_2-jyseps1_2, ysub);
            continue;
          }
          
          // Check upper legs
          if( (ny_inner - jyseps2_1-1) % ysub != 0 ) {
            output.write("\t -> leg region ny_inner-jyseps2_1-1 (%d-%d-1 = %d) must be a multiple of MYSUB (%d)\n", 
                         ny_inner, jyseps2_1, ny_inner-jyseps2_1-1, ysub);
            continue;
          }
          if( (jyseps1_2-ny_inner+1) % ysub != 0 ) {
            output.write("\t -> leg region jyseps1_2-ny_inner+1 (%d-%d+1 = %d) must be a multiple of MYSUB (%d)\n", 
                         jyseps1_2, ny_inner, jyseps1_2-ny_inner+1, ysub);
            continue;
          }
        }else {
          // Single Null
          if( (jyseps2_2-jyseps1_1) % ysub != 0 ) {
            output.write("\t -> Core region jyseps2_2-jyseps1_1 (%d-%d = %d) must be a multiple of MYSUB (%d)\n", 
                         jyseps2_2, jyseps1_1, jyseps2_2-jyseps1_1, ysub);
            continue;
          }
        }
        
        if( (ny - jyseps2_2 - 1) % ysub != 0) {
          output.write("\t -> leg region ny-jyseps2_2-1 (%d-%d-1 = %d) must be a multiple of MYSUB (%d)\n", 
                       ny, jyseps2_2, ny - jyseps2_2 - 1, ysub);
          continue;
        }
        output.write("\t -> Good value\n");
        // Found an acceptable value
        if((NXPE < 1) || 
           (fabs(ideal - i) < fabs(ideal - NXPE)))
          NXPE = i; // Keep value nearest to the ideal
      }
    }
    
    if(NXPE < 1)
      throw BoutException("Could not find a valid value for NXPE");
    
    NYPE = NPES / NXPE;
    
    output.write("\tDomain split (%d, %d) into domains (%d, %d)\n",
                 NXPE, NYPE, MX / NXPE, ny / NYPE);
  }
  
  /// Get X and Y processor indices
  PE_YIND = MYPE / NXPE;
  PE_XIND = MYPE % NXPE;
  
  // Work out other grid size quantities

  /// MXG at each end needed for edge boundary regions
  MX = nx - 2*MXG;
  
  /// Split MX points between NXPE processors
  MXSUB = MX / NXPE;
  if((MX % NXPE) != 0) {
    throw BoutException("Cannot split %d X points equally between %d processors\n",
                            MX, NXPE);
  }

  /// NOTE: No grid data reserved for Y boundary cells - copy from neighbours
  MY = ny;
  MYSUB = MY / NYPE;
  if((MY % NYPE) != 0) {
    throw BoutException("\tERROR: Cannot split %d Y points equally between %d processors\n",
                        MY, NYPE);
  }
  
  /// Get mesh options
  OPTION(options, non_uniform,  false);
  OPTION(options, TwistShift,   false);
  OPTION(options, TwistOrder,   0);
  OPTION(options, ShiftOrder,   0);
  OPTION(options, ShiftXderivs, false);
  OPTION(options, IncIntShear,  false);
  OPTION(options, BoundaryOnCell, false); // Determine location of boundary
  OPTION(options, StaggerGrids,   false); // Stagger grids
  OPTION(options, periodicX, false); // Periodic in X
  
  OPTION(options, async_send, false); // Whether to use asyncronous sends
  
  // Set global offsets
  
  OffsetX = PE_XIND*MXSUB;
  OffsetY = PE_YIND*MYSUB;
  OffsetZ = 0;
  
  if(ShiftXderivs) {
    output.write("Using shifted X derivatives. Interpolation: ");
    if(ShiftOrder == 0) {
      output.write("FFT\n");
    }else
      output.write("%d-point\n", ShiftOrder);
  }
  
  if(options->isSet("zperiod")) {
    OPTION(options, zperiod, 1);
    ZMIN = 0.0;
    ZMAX = 1.0 / (double) zperiod;
  }else {
    OPTION(options, ZMIN, 0.0);
    OPTION(options, ZMAX, 1.0);
    
    zperiod = ROUND(1.0 / (ZMAX - ZMIN));
  }

  if(TwistShift) {
    output.write("Applying Twist-Shift condition. Interpolation: ");
    if(TwistOrder == 0) {
      output.write("FFT\n");
    }else
      output.write("%d-point\n", TwistOrder);
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
    throw BoutException("\tERROR: Diagonal metrics are not finite!\n");
  }
  if((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal metrics are negative!\n");
  }
  if((!finite(g12)) || (!finite(g13)) || (!finite(g23))) {
    throw BoutException("\tERROR: Off-diagonal metrics are not finite!\n");
  }

  /// Set shift for radial derivatives
  if(get(zShift, "zShift")) {
    output.write("\tWARNING: Z shift for radial derivatives not found\n");
    ShiftTorsion = zShift = 0.0;
  }else if(get(ShiftTorsion, "ShiftTorsion")) {
    output.write("\tWARNING: No Torsion specified for zShift. Derivatives may not be correct\n");
    ShiftTorsion = 0.0;
  }
  
  if(IncIntShear) {
    if(get(IntShiftTorsion, "IntShiftTorsion")) {
      output.write("\tWARNING: No Integrated torsion specified\n");
      IntShiftTorsion = 0.0;
    }
  }
  // Allocate some memory for twist-shift

  ShiftAngle  = new BoutReal[ngx];

  // Try to read the shift angle from the grid file
  // NOTE: All processors should know the twist-shift angle (for invert_parderiv)
  if(source->hasVar("ShiftAngle")) {
    source->open("ShiftAngle");
    source->setGlobalOrigin(XGLOBAL(0));
    if(!source->fetch(ShiftAngle,  "ShiftAngle", ngx)) {
      output.write("\tWARNING: Twist-shift angle 'ShiftAngle' not found. Setting to zero\n");
      for(int i=0;i<ngx;i++)
        ShiftAngle[i] = 0.0;
    }
    source->close();
  }else if(MYG > 0) {
    output.write("\tWARNING: Twist-shift angle 'ShiftAngle' not found. Setting from zShift\n");

    if(YPROC(jyseps2_2) == PE_YIND) {
      for(int i=0;i<ngx;i++)
	ShiftAngle[i] = zShift[i][MYG+MYSUB-1] - zShift[i][MYG+MYSUB]; // Jump across boundary
   
    }else if(YPROC(jyseps1_1+1) == PE_YIND) {
      for(int i=0;i<ngx;i++)
	ShiftAngle[i] = zShift[i][MYG-1] - zShift[i][MYG]; // Jump across boundary
    }
    
    msg_stack.push("Creating core_comm for ShiftAngle");
    
    // In the core, need to set ShiftAngle everywhere for ballooning initial condition
    MPI_Group groupw;
    MPI_Comm_group(BoutComm::get(), &groupw); // Group of all processors
    
    int *ranks = new int[NYPE];
    int npcore = 0;
    for(int p = YPROC(jyseps1_1+1); p <= YPROC(jyseps2_2);p++) {
      ranks[npcore] = PROC_NUM(PE_XIND, p);
      npcore++;
    }
    
    MPI_Group grp;
    MPI_Group_incl(groupw, npcore, ranks, &grp); // Create group
    
    MPI_Comm core_comm;
    MPI_Comm_create(BoutComm::get(), grp, &core_comm); // Create communicator
    
    delete[] ranks;
    
    if(MYPE_IN_CORE)
      MPI_Bcast(ShiftAngle, ngx, PVEC_REAL_MPI_TYPE, npcore-1, core_comm);
    
    // Free MPI handles
    if(core_comm != MPI_COMM_NULL)
      MPI_Comm_free(&core_comm);
    MPI_Group_free(&grp);
    MPI_Group_free(&groupw);
    
    msg_stack.pop();
  }

  /// Can have twist-shift in the private flux regions too
  bool twistshift_pf;
  OPTION(options, twistshift_pf, false);
  if(twistshift_pf) {
    output << "Adding twist-shift in lower PF region" << endl;
    // Lower PF. Note by default no Twist-Shift used here, so need to switch on
    if(YPROC(jyseps1_1) == PE_YIND) {
      for(int i=0;i<ngx;i++) {
	ShiftAngle[i] = zShift[i][MYG+MYSUB-1] - zShift[i][MYG+MYSUB]; // Jump across boundary
      }
      TS_up_in = true; // Switch on twist-shift
      
    }else if(YPROC(jyseps2_2+1) == PE_YIND) {
      for(int i=0;i<ngx;i++) {
	ShiftAngle[i] = zShift[i][MYG-1] - zShift[i][MYG]; // Jump across boundary
      }
      TS_down_in = true;
    }
  }

  /// Calculate contravariant metric components
  if(calcCovariant())
    throw BoutException("Error in calcCovariant call");

  /// Calculate Jacobian and Bxy
  if(jacobian())
    throw BoutException("Error in jacobian call");
  
  // Attempt to read J from the grid file
  Field2D Jcalc = J;
  if(get(J, "J")) {
    output.write("\tWARNING: Jacobian 'J' not found. Calculating from metric tensor\n");
    J = Jcalc;
  }else {
    // Compare calculated and loaded values  
    output.write("\tMaximum difference in J is %e\n", max(abs(J - Jcalc)));
    
    // Re-evaluate Bxy using new J
    Bxy = sqrt(g_22)/J;
  }

  // Attempt to read Bxy from the grid file
  Field2D Bcalc = Bxy;
  if(get(Bxy, "Bxy")) {
    output.write("\tWARNING: Magnitude of B field 'Bxy' not found. Calculating from metric tensor\n");
    Bxy = Bcalc;
  }else {
    output.write("\tMaximum difference in Bxy is %e\n", max(abs(Bxy - Bcalc)));
    // Check Bxy
    if(!finite(Bxy))
      throw BoutException("\tERROR: Bxy not finite everywhere!\n");
  }

  //////////////////////////////////////////////////////
  /// Communicator
  
  MPI_Group group_world;
  MPI_Comm_group(BoutComm::get(), &group_world); // Get the entire group
  
  //////////////////////////////////////////////////////
  /// Communicator in X
  
  MPI_Group group;
  MPI_Comm comm_tmp;
  
  int proc[3]; // Processor range
  
  for(int yp = 0; yp < NYPE; yp++) {
    proc[0] = PROC_NUM(0, yp);      // First 
    proc[1] = PROC_NUM(NXPE-1, yp); // Last
    proc[2] = 1;                    // stride
    
#ifdef COMMDEBUG
    output << "XCOMM " << proc[0] << ", " << proc[1] << endl;
#endif
    if(MPI_Group_range_incl(group_world, 1, &proc, &group) != MPI_SUCCESS)
      throw BoutException("Could not create X communication group for yp=%d (xind=%d,yind=%d)\n",
			  yp, PE_XIND, PE_YIND);
    if(MPI_Comm_create(BoutComm::get(), group, &comm_tmp) != MPI_SUCCESS)
      throw BoutException("Could not create X communicator for yp=%d (xind=%d,yind=%d)\n", 
			  yp, PE_XIND, PE_YIND);
    MPI_Group_free(&group);
    
    if(yp == PE_YIND) {
      // Should be in this group
      if(comm_tmp == MPI_COMM_NULL)
        throw BoutException("X communicator null");
      
      comm_x = comm_tmp;
    }else {
      if(comm_tmp != MPI_COMM_NULL)
        throw BoutException("X communicator should be null");
    }
  }
  
  //////////////////////////////////////////////////////
  /// Communicators for Y gather/scatter
  
  MPI_Group group_tmp1, group_tmp2;
  
  proc[2] = NXPE; // Stride in processor rank

  // Outer SOL regions
  if(jyseps1_2 == jyseps2_1) {
    // Single-null. All processors with same PE_XIND

    msg_stack.push("Creating Outer SOL communicators for Single Null operation");
    
    for(int i=0;i<NXPE;i++) {
      proc[0] = PROC_NUM(i, 0);
      proc[1] = PROC_NUM(i, NYPE-1);
#ifdef COMMDEBUG
    output << "Outer SOL " << proc[0] << ", " << proc[1] << endl;
#endif
      MPI_Group_range_incl(group_world, 1, &proc, &group);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if(i == PE_XIND) {
	// Should be part of this communicator
	if(comm_tmp == MPI_COMM_NULL) {
	  // error
	  throw BoutException("Single null outer SOL not correct\n");
	}
	comm_outer = comm_tmp;
      }else if(comm_tmp != MPI_COMM_NULL) {
	// Not part of this communicator so should be NULL
	throw BoutException("Single null outer SOL not correct\n");
      }
      MPI_Group_free(&group);
    }
    msg_stack.pop();
  }else {
    // Double null
    
    msg_stack.push("Creating Outer SOL communicators for Double Null operation");
    
    for(int i=0;i<NXPE;i++) {
      // Inner SOL
      proc[0] = PROC_NUM(i, 0);
      proc[1] = PROC_NUM(i, YPROC(ny_inner-1));
#ifdef COMMDEBUG
    output << "Double Null inner SOL " << proc[0] << ", " << proc[1] << endl;
#endif
      if(MPI_Group_range_incl(group_world, 1, &proc, &group) != MPI_SUCCESS)
	throw BoutException("MPI_Group_range_incl failed for xp = %d", NXPE);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if(comm_tmp != MPI_COMM_NULL)
	comm_outer = comm_tmp;
      MPI_Group_free(&group);

      // Outer SOL
      proc[0] = PROC_NUM(i, YPROC(ny_inner));
      proc[1] = PROC_NUM(i, NYPE-1);
#ifdef COMMDEBUG
    output << "Double Null outer SOL " << proc[0] << ", " << proc[1] << endl;
#endif
      MPI_Group_range_incl(group_world, 1, &proc, &group);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if(comm_tmp != MPI_COMM_NULL) {
	comm_outer = comm_tmp;
      }
      MPI_Group_free(&group);
    }
    
    msg_stack.pop();
  }

  for(int i=0;i<NXPE;i++) {
    // Lower PF region

    if((jyseps1_1 >= 0) || (jyseps2_2+1 < ny)) {
      // A lower PF region exists

#ifdef COMMDEBUG
      output << "Creating lower PF communicators for xp = " << i << endl;
#endif

      msg_stack.push("Creating lower PF communicators for xp=%d", i);

      if(jyseps1_1 >= 0) {
	proc[0] = PROC_NUM(i, 0);
	proc[1] = PROC_NUM(i, YPROC(jyseps1_1));
#ifdef COMMDEBUG
	output << "PF1 "<< proc[0] << ", " << proc[1] << endl;
#endif
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
      }else
	group_tmp1 = MPI_GROUP_EMPTY;
      
      if(jyseps2_2+1 < ny) {
	proc[0] = PROC_NUM(i, YPROC(jyseps2_2+1));
	proc[1] = PROC_NUM(i, NYPE-1);
#ifdef COMMDEBUG
	output << "PF2 "<< proc[0] << ", " << proc[1] << endl;
#endif
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
      }else
	group_tmp2 = MPI_GROUP_EMPTY;
      
      MPI_Group_union(group_tmp1, group_tmp2, &group);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if(comm_tmp != MPI_COMM_NULL) {
	comm_inner = comm_tmp;
	if(ixseps_lower == ixseps_outer) {
	  // Between the separatrices is still in the PF region
#ifdef COMMDEBUG
          output << "-> Inner and middle\n";
#endif
          comm_middle = comm_inner;
          //MPI_Comm_dup(comm_inner, &comm_middle);
	}else {
#ifdef COMMDEBUG
          output << "-> Outer and middle\n";
#endif
          comm_middle = comm_outer;
          //MPI_Comm_dup(comm_outer, &comm_middle); // Error! Needs to be collective on comm_outer
        }
      }
#ifdef COMMDEBUG
      output << "Freeing\n";
#endif
      MPI_Group_free(&group);
      if(group_tmp1 != MPI_GROUP_EMPTY) {
        MPI_Group_free(&group_tmp1);
      }
      if(group_tmp2 != MPI_GROUP_EMPTY) {
        MPI_Group_free(&group_tmp2);
      }
#ifdef COMMDEBUG
      output << "done lower PF\n";
#endif
      msg_stack.pop();
    }

    if(jyseps2_1 != jyseps1_2) {
      // Upper PF region
      // Note need to order processors so that a continuous surface is formed
      
#ifdef COMMDEBUG
      output << "Creating upper PF communicators for xp = " << i << endl;
#endif
      msg_stack.push("Creating upper PF communicators for xp=%d", i);

      proc[0] = PROC_NUM(i, YPROC(ny_inner));
      proc[1] = PROC_NUM(i, YPROC(jyseps1_2));
#ifdef COMMDEBUG
      output << "PF3 "<< proc[0] << ", " << proc[1] << endl;
#endif
      MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
      proc[0] = PROC_NUM(i, YPROC(jyseps2_1+1));
      proc[1] = PROC_NUM(i, YPROC(ny_inner-1));
#ifdef COMMDEBUG
      output << "PF4 "<< proc[0] << ", " << proc[1] << endl;
#endif
      MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
      MPI_Group_union(group_tmp1, group_tmp2, &group);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if(comm_tmp != MPI_COMM_NULL) {
	comm_inner = comm_tmp;
	if(ixseps_upper == ixseps_outer) {
#ifdef COMMDEBUG
          output << "-> Inner and middle\n";
#endif
          comm_middle = comm_inner;
          //MPI_Comm_dup(comm_inner, &comm_middle);
	}else {
#ifdef COMMDEBUG
          output << "-> Outer and middle\n";
#endif
          comm_middle = comm_outer;
          //MPI_Comm_dup(comm_outer, &comm_middle);
        }
      }
#ifdef COMMDEBUG
      output << "Freeing\n";
#endif
      MPI_Group_free(&group);
      if(group_tmp1 != MPI_GROUP_EMPTY)
        MPI_Group_free(&group_tmp1);
      if(group_tmp2 != MPI_GROUP_EMPTY)
        MPI_Group_free(&group_tmp2);
      msg_stack.pop();
#ifdef COMMDEBUG
      output << "done upper PF\n";
#endif
    }
    
    // Core region
    msg_stack.push("Creating core communicators");
    proc[0] = PROC_NUM(i, YPROC(jyseps1_1+1));
    proc[1] = PROC_NUM(i, YPROC(jyseps2_1));
#ifdef COMMDEBUG
    output << "CORE1 "<< proc[0] << ", " << proc[1] << endl;
#endif 
    if( (proc[0] < 0) || (proc[1] < 0) )
      throw BoutException("Invalid processor range for core processors");
    MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
    
    proc[0] = PROC_NUM(i, YPROC(jyseps1_2+1));
    proc[1] = PROC_NUM(i, YPROC(jyseps2_2));
#ifdef COMMDEBUG
    output << "CORE2 "<< proc[0] << ", " << proc[1] << endl;
#endif
    if( (proc[0] < 0) || (proc[1] < 0) ) {
      group_tmp2 = MPI_GROUP_EMPTY;
    }else {
      MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
    }
    
    MPI_Group_union(group_tmp1, group_tmp2, &group);
    MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
    if(comm_tmp != MPI_COMM_NULL) {
      comm_inner = comm_tmp;
      
      if(ixseps_inner == ixseps_outer)
        MPI_Comm_dup(comm_inner, &comm_middle);
    }
    
    if(group_tmp1 != MPI_GROUP_EMPTY)
      MPI_Group_free(&group_tmp1);
    if(group_tmp2 != MPI_GROUP_EMPTY)
      MPI_Group_free(&group_tmp2);
    MPI_Group_free(&group);
    
    msg_stack.pop();
  }
  
  if(ixseps_inner == ixseps_outer) {
    // Balanced null, so no middle
    MPI_Comm_dup(comm_inner, &comm_middle);
  }else {
    // Need to handle unbalanced double-null case
    
#ifdef COMMDEBUG
    output << "Unbalanced " << endl;
#endif

    if(ixseps_upper > ixseps_lower) {
      // middle is connected to the bottom
      
      msg_stack.push("Creating unbalanced lower communicators");
      
      for(int i=0;i<NXPE;i++) {
	proc[0] = PROC_NUM(i, 0);
	proc[1] = PROC_NUM(i, YPROC(jyseps2_1));
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
	proc[0] = PROC_NUM(i, YPROC(jyseps1_2+1));
	proc[1] = PROC_NUM(i, NYPE-1);
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
	MPI_Group_union(group_tmp1, group_tmp2, &group);
	MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
	if(comm_tmp != MPI_COMM_NULL)
	  comm_middle = comm_tmp;
        
        if(group_tmp1 != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_tmp1);
        if(group_tmp2 != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_tmp2);
        MPI_Group_free(&group);
      }
      msg_stack.pop();
    }else {
      // middle is connected to the top
      
      msg_stack.push("Creating unbalanced upper communicators");
      for(int i=0;i<NXPE;i++) {
	proc[0] = PROC_NUM(i, YPROC(ny_inner));
	proc[1] = PROC_NUM(i, YPROC(jyseps2_2));
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
	proc[0] = PROC_NUM(i, YPROC(jyseps1_1+1));
	proc[1] = PROC_NUM(i, YPROC(ny_inner-1));
	MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
	MPI_Group_union(group_tmp1, group_tmp2, &group);
	MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
	if(comm_tmp != MPI_COMM_NULL)
	  comm_middle = comm_tmp;
        
        if(group_tmp1 != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_tmp1);
        if(group_tmp2 != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_tmp2);
        MPI_Group_free(&group);
      }
      msg_stack.pop();
    }
  }
  MPI_Group_free(&group_world);
  // Now have communicators for all regions.

#ifdef COMMDEBUG
  output << "Got communicators" << endl;
#endif

  //////////////////////////////////////////////////////
  /// Calculate Christoffel symbols. Needs communication
  if(geometry()) {
    throw BoutException("Differential geometry failed\n");
  }

  if(periodicX) {
    FieldGroup g;
    g.add(zShift, dx);
    communicate(g);
  }

  //////////////////////////////////////////////////////
  /// Non-uniform meshes. Need to use DDX, DDY
  
  if(non_uniform) {
    Field2D d2x, d2y; // d^2 x / d i^2
    // Read correction for non-uniform meshes
    if(get(d2x, "d2x")) {
      output.write("\tWARNING: differencing quantity 'd2x' not found. Calculating from dx\n");
      d1_dx = DDX(1./dx)*dx; // d/di(1/dx)
    }else
      d1_dx = -d2x / (dx*dx);
    
    if(get(d2y, "d2y")) {
      output.write("\tWARNING: differencing quantity 'd2y' not found. Calculating from dy\n");
      d1_dy = DDY(1./dy)*dy; // d/di(1/dy)
    }else
      d1_dy = -d2y / (dy*dy);
  }

  //////////////////////////////////////////////////////
  // Boundary regions
  if(!periodicX) {
    // Need boundaries in X if not periodic
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
  }
  
  if((UDATA_INDEST < 0) && (UDATA_XSPLIT > xstart))
    boundary.push_back(new BoundaryRegionYUp("target", xstart, UDATA_XSPLIT-1));
  if((UDATA_OUTDEST < 0) && (UDATA_XSPLIT <= xend))
    boundary.push_back(new BoundaryRegionYUp("target", UDATA_XSPLIT, xend));
  
  if((DDATA_INDEST < 0) && (DDATA_XSPLIT > xstart))
    boundary.push_back(new BoundaryRegionYDown("target", xstart, DDATA_XSPLIT-1));
  if((DDATA_OUTDEST < 0) && (DDATA_XSPLIT <= xend))
    boundary.push_back(new BoundaryRegionYDown("target", DDATA_XSPLIT, xend));
    
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

int BoutMesh::get(Field2D &var, const char *name, BoutReal def) {
  if(name == NULL)
    return 1;
  
  int msg_pos = msg_stack.push("Loading 2D field: BoutMesh::get(Field2D, %s)", name);
  
  if(!source->hasVar(name)) {
    output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name, def);
    var = def;
    msg_stack.pop(msg_pos);
    return 2;
  }
  
  BoutReal **data;
  int jy;
  
  var = 0;// Makes sure the memory is allocated
  
  data = var.getData(); // Get a pointer to the data
  
  // Send an open signal to the source
  source->open(name);
  
  // Get the size of the variable
  vector<int> size = source->getSize(name);
  switch(size.size()) {
  case 1: {
    // 0 or 1 dimension
    if(size[0] != 1) {
      output.write("Expecting a 2D variable, but '%s' is 1D with %d elements\n", name, size[0]);
      source->close();
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 1;
    }
    BoutReal val;
    if(!source->fetch(&val, name)) {
      output.write("Couldn't read 0D variable '%s'\n", name);
      source->close();
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 1;
    }
    
    var = val;
    
    // Close source
    source->close();
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 0;
  }
  case 2: {
    if((size[0] != nx) || (size[1] != ny)) {
      output.write("Error: Variable '%s' has dimensions [%d,%d]. Expecting [%d,%d]\n",
                   name, size[0], size[1], nx, ny);
      source->close();
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 1;
    }
    break;
  }
  default: {
    output.write("Error: Variable '%s' should be 2D, but has %d dimensions\n", 
                 name, size.size());
    source->close();
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 1;
  }
  }
  
  // Read in data for bulk of points
  
  if(readgrid_2dvar(source, name,
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
    
    if(readgrid_2dvar(source, name,
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
    
    if(readgrid_2dvar(source, name,
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

    if(readgrid_2dvar(source, name,
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

    if(readgrid_2dvar(source, name,
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
  
  /*
  if(IDATA_DEST >= 0) {
    int xfrom = (IDATA_DEST % NXPE)*MXSUB + MXSUB;
    int yfrom = (IDATA_DEST/NXPE)*MYSUB;
    for(int i=0;i<MXG;i++) {
      output.write("in: (%d,%d) -> (%d,%d)\n", xfrom+i, yfrom, i, MYG);
      s->setGlobalOrigin(xfrom+i, yfrom);
      s->fetch(&(data[i][MYG]), name, 1, MYSUB);
    }
    s->setGlobalOrigin();
  }
  if(ODATA_DEST >= 0) {
    int xfrom = (ODATA_DEST % NXPE)*MXSUB + MXG;
    int yfrom = (ODATA_DEST/NXPE)*MYSUB;
    for(int i=0;i<MXG;i++) {
      output.write("out: (%d,%d) -> (%d,%d)\n", xfrom+i, yfrom, MXG+MXSUB+i, MYG);
      s->setGlobalOrigin(xfrom+i, yfrom);
      s->fetch(&(data[MXG+MXSUB+i][MYG]), name, 1, MYSUB);
    }
    s->setGlobalOrigin();
  }
  */
#ifdef TRACK
  var.name = copy_string(name);
#endif
   
  // Close source
  source->close();
  
#ifdef CHECK
  // Check that the data is ok
  var.checkData(true);
  
  msg_stack.pop(msg_pos);
#endif
  
  return 0;
}

int BoutMesh::get(Field2D &var, const string &name, BoutReal def) {
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
int BoutMesh::get(Field3D &var, const char *name) {
  if(name == NULL)
    return 1;
  
#ifdef CHECK
  msg_stack.push("Loading 3D field: BoutMesh::get(Field3D, %s)", name);
#endif
  
  if(!source->hasVar(name)) {
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
  source->open(name);

  // Read in data for bulk of points
  if(readgrid_3dvar(source, name,
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
    
    if(readgrid_3dvar(source, name,
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
    
    if(readgrid_3dvar(source, name,
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

    if(readgrid_3dvar(source, name,
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

    if(readgrid_3dvar(source, name,
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

int BoutMesh::get(Field3D &var, const string &name) {
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

int BoutMesh::communicate(FieldGroup &g) {
  comm_handle c = send(g);
  return wait(c);
}

void BoutMesh::post_receive(CommHandle &ch) {
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
	      BoutComm::get(),
	      &ch.request[0]);
  }
  if(UDATA_OUTDEST != -1) {
    inbuff = &ch.umsg_recvbuff[len]; // pointer to second half of the buffer
    MPI_Irecv(inbuff,
	      msg_len(ch.var_list, UDATA_XSPLIT, ngx, 0, MYG),
	      PVEC_REAL_MPI_TYPE,
	      UDATA_OUTDEST,
	      OUT_SENT_DOWN,
	      BoutComm::get(),
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
	      BoutComm::get(),
	      &ch.request[2]);
  }
  if(DDATA_OUTDEST != -1) {
    inbuff = &ch.dmsg_recvbuff[len];
    MPI_Irecv(inbuff,
	      msg_len(ch.var_list, DDATA_XSPLIT, ngx, 0, MYG),
	      PVEC_REAL_MPI_TYPE,
	      DDATA_OUTDEST,
	      OUT_SENT_UP,
	      BoutComm::get(),
	      &ch.request[3]);
  }

  /// Post receive data from left (x-1)
  
  if(IDATA_DEST != -1) {
    MPI_Irecv(ch.imsg_recvbuff,
	      msg_len(ch.var_list, 0, MXG, 0, MYSUB),
	      PVEC_REAL_MPI_TYPE,
	      IDATA_DEST,
	      OUT_SENT_IN,
	      BoutComm::get(),
	      &ch.request[4]);
  }

  // Post receive data from right (x+1)

  if(ODATA_DEST != -1) {
    MPI_Irecv(ch.omsg_recvbuff,
	      msg_len(ch.var_list, 0, MXG, 0, MYSUB),
	      PVEC_REAL_MPI_TYPE,
	      ODATA_DEST,
	      IN_SENT_OUT,
	      BoutComm::get(),
	      &ch.request[5]);
  }
}

comm_handle BoutMesh::send(FieldGroup &g) { 
  /// Start timer
  Timer timer("comms");
  
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
		BoutComm::get(),
		&(ch->sendreq[0]));
    }else
      MPI_Send(ch->umsg_sendbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       UDATA_INDEST,
	       IN_SENT_UP,
	       BoutComm::get());
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
		BoutComm::get(),
		&(ch->sendreq[1]));
    }else
      MPI_Send(outbuff, 
	       len, 
	       PVEC_REAL_MPI_TYPE,
	       UDATA_OUTDEST,
	       OUT_SENT_UP,
	       BoutComm::get());
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
		BoutComm::get(),
		&(ch->sendreq[2]));
    }else
      MPI_Send(ch->dmsg_sendbuff, 
	       len,
	       PVEC_REAL_MPI_TYPE,
	       DDATA_INDEST,
	       IN_SENT_DOWN,
	       BoutComm::get());
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
		BoutComm::get(),
		&(ch->sendreq[3]));
    }else
      MPI_Send(outbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       DDATA_OUTDEST,
	       OUT_SENT_DOWN,
	       BoutComm::get());
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
		BoutComm::get(),
		&(ch->sendreq[4]));
    }else
      MPI_Send(ch->imsg_sendbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       IDATA_DEST,
	       IN_SENT_OUT,
	       BoutComm::get());
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
		BoutComm::get(),
		&(ch->sendreq[5]));
    }else
      MPI_Send(ch->omsg_sendbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       ODATA_DEST,
	       OUT_SENT_IN,
	       BoutComm::get());
  }
  
  /// Mark communication handle as in progress
  ch->in_progress = true;
  
  return (void*) ch;
}

int BoutMesh::wait(comm_handle handle) {
  if(handle == NULL)
    return 1;
  
  CommHandle *ch = (CommHandle*) handle;

  if(!ch->in_progress)
    return 2;

  /// Start timer
  Timer timer("comms");
  
  ///////////// WAIT FOR DATA //////////////
  
  int ind, len;
  MPI_Status status;

  if(ch->var_list.size() == 0) {
    
    // Just waiting for a single MPI request
    MPI_Wait(ch->request, &status);
    free_handle(ch);
    
    return 0;
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
  if(TwistShift && (TwistOrder == 0)) {
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
  
  return 0;
}

/***************************************************************
 *             Non-Local Communications
 ***************************************************************/

MPI_Request BoutMesh::sendToProc(int xproc, int yproc, BoutReal *buffer, int size, int tag) {
  Timer timer("comms");
  
  MPI_Request request;
  
  MPI_Isend(buffer, size, PVEC_REAL_MPI_TYPE,
	   PROC_NUM(xproc,yproc),
	   tag,
	   BoutComm::get(),
	   &request);
	   
  return request;
}

comm_handle BoutMesh::receiveFromProc(int xproc, int yproc, BoutReal *buffer, int size, int tag) {
  Timer timer("comms");
  
  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0,0);
  
  MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    PROC_NUM(xproc,yproc),
	    tag,
	    BoutComm::get(),
	    ch->request);
	    
  ch->in_progress = true;
  
  return (comm_handle) ch;
}

int BoutMesh::getNXPE() {
  return NXPE;
}

int BoutMesh::getNYPE() {
  return NYPE;
}

int BoutMesh::getXProcIndex() {
  return PE_XIND;
}

int BoutMesh::getYProcIndex() {
  return PE_YIND;
}

/****************************************************************
 *                 X COMMUNICATIONS
 * 
 * Intended mainly to handle the perpendicular inversion operators
 ****************************************************************/

bool BoutMesh::firstX() {
  return PE_XIND == 0;
}

bool BoutMesh::lastX() {
  return PE_XIND == NXPE-1;
}

int BoutMesh::sendXOut(BoutReal *buffer, int size, int tag) {
  if(PE_XIND == NXPE-1)
    return 1;
  
  Timer timer("comms");

  MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   PROC_NUM(PE_XIND+1, PE_YIND),
	   tag,
	   BoutComm::get());

  return 0;
}

int BoutMesh::sendXIn(BoutReal *buffer, int size, int tag) {
  if(PE_XIND == 0)
    return 1;
  
  Timer timer("comms");

  MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   PROC_NUM(PE_XIND-1, PE_YIND),
	   tag,
	   BoutComm::get());

  return 0;
}

comm_handle BoutMesh::irecvXOut(BoutReal *buffer, int size, int tag) {
  if(PE_XIND == NXPE-1)
    return NULL;

  Timer timer("comms");
  
  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0,0);

  MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    PROC_NUM(PE_XIND+1, PE_YIND),
	    tag,
	    BoutComm::get(),
	    ch->request);
  
  ch->in_progress = true;

  return (comm_handle) ch;
}

comm_handle BoutMesh::irecvXIn(BoutReal *buffer, int size, int tag) {
  if(PE_XIND == 0)
    return NULL;
  
  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0,0);

  MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    PROC_NUM(PE_XIND-1, PE_YIND),
	    tag,
	    BoutComm::get(),
	    ch->request);
  
  ch->in_progress = true;

  return (comm_handle) ch;
}

/****************************************************************
 *                 Y COMMUNICATIONS
 * 
 * Intended mainly to handle the non-local heat flux integrations
 ****************************************************************/

bool BoutMesh::firstY() {
  return PE_YIND == 0;
}

bool BoutMesh::lastY() {
  return PE_YIND == NYPE-1;
}

int BoutMesh::UpXSplitIndex() {
  return UDATA_XSPLIT;
}

int BoutMesh::DownXSplitIndex() {
  return DDATA_XSPLIT;
}

int BoutMesh::sendYOutIndest(BoutReal *buffer, int size, int tag) {
  if(PE_YIND == NYPE-1)
    return 1;
  
  Timer timer("comms");

  if (UDATA_INDEST != -1)
    MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   UDATA_INDEST,
	   tag,
	   BoutComm::get());
  else throw BoutException("Expected UDATA_INDEST to exist, but it does not.");
  return 0;
}

int BoutMesh::sendYOutOutdest(BoutReal *buffer, int size, int tag) {
  if(PE_YIND == NYPE-1)
    return 1;
  
  Timer timer("comms");
	   
  if (UDATA_OUTDEST != -1)
    MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   UDATA_OUTDEST,
	   tag,
	   BoutComm::get());
  else throw BoutException("Expected UDATA_OUTDEST to exist, but it does not.");

  return 0;
}

int BoutMesh::sendYInIndest(BoutReal *buffer, int size, int tag) {
  if(PE_YIND == 0)
    return 1;
  
  Timer timer("comms");

  if (DDATA_INDEST != -1)
    MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   DDATA_INDEST,
	   tag,
	   BoutComm::get());
  else throw BoutException("Expected DDATA_INDEST to exist, but it does not.");

  return 0;
}

int BoutMesh::sendYInOutdest(BoutReal *buffer, int size, int tag) {
  if(PE_YIND == 0)
    return 1;
  
  Timer timer("comms");

  if (DDATA_OUTDEST != -1)
    MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   DDATA_OUTDEST,
	   tag,
	   BoutComm::get());
  else throw BoutException("Expected DDATA_OUTDEST to exist, but it does not.");

  return 0;
}

comm_handle BoutMesh::irecvYOutIndest(BoutReal *buffer, int size, int tag) {
  if(PE_YIND == NYPE-1)
    return NULL;

  Timer timer("comms");
  
  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0,0);

  if (UDATA_INDEST != -1)
    MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    UDATA_INDEST,
	    tag,
	    BoutComm::get(),
	    ch->request);
  else throw BoutException("Expected UDATA_INDEST to exist, but it does not.");
  
  ch->in_progress = true;

  return (comm_handle) ch;
}

comm_handle BoutMesh::irecvYOutOutdest(BoutReal *buffer, int size, int tag) {
  if(PE_YIND == NYPE-1)
    return NULL;

  Timer timer("comms");
  
  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0,0);

  if (UDATA_OUTDEST != -1)
    MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    UDATA_OUTDEST,
	    tag,
	    BoutComm::get(),
	    ch->request);
  else throw BoutException("Expected UDATA_OUTDEST to exist, but it does not.");
  
  ch->in_progress = true;

  return (comm_handle) ch;
}

comm_handle BoutMesh::irecvYInIndest(BoutReal *buffer, int size, int tag) {
  if(PE_YIND == 0)
    return NULL;
  
  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0,0);

  if (DDATA_INDEST != -1)
    MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    DDATA_INDEST,
	    tag,
	    BoutComm::get(),
	    ch->request);
  else throw BoutException("Expected DDATA_INDEST to exist, but it does not.");
  
  ch->in_progress = true;

  return (comm_handle) ch;
}

comm_handle BoutMesh::irecvYInOutdest(BoutReal *buffer, int size, int tag) {
  if(PE_YIND == 0)
    return NULL;
  
  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0,0);

  if (DDATA_OUTDEST != -1)
    MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    DDATA_OUTDEST,
	    tag,
	    BoutComm::get(),
	    ch->request);
  else throw BoutException("Expected DDATA_OUTDEST to exist, but it does not.");
  
  ch->in_progress = true;

  return (comm_handle) ch;
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
int BoutMesh::XGLOBAL(int xloc) const {
  return xloc + PE_XIND * MXSUB;
}

/// Returns a local X index given a global index
int BoutMesh::XLOCAL(int xglo) const {
  return xglo - PE_XIND * MXSUB;
}

/// Returns the global Y index given a local index
int BoutMesh::YGLOBAL(int yloc) const {
  return yloc + PE_YIND*MYSUB - MYG;
}

/// Global Y index given local index and processor
int BoutMesh::YGLOBAL(int yloc, int yproc) const {
  return yloc + yproc*MYSUB - MYG;
}

/// Returns a local Y index given a global index
int BoutMesh::YLOCAL(int yglo) const {
  return yglo - PE_YIND*MYSUB + MYG;
}

int BoutMesh::YLOCAL(int yglo, int yproc) const {
  return yglo - yproc*MYSUB + MYG;
}

/// Return the Y processor number given a global Y index
int BoutMesh::YPROC(int yind) {
  if((yind < 0) || (yind > ny))
    return -1;
  return yind / MYSUB;
}

/// Return the X processor number given a global X index
int BoutMesh::XPROC(int xind) {
  return (xind >= MXG) ? (xind - MXG) / MXSUB : 0;
}

/****************************************************************
 *                       CONNECTIONS
 ****************************************************************/

/// Connection initialisation: Set processors in a simple 2D grid
void BoutMesh::default_connections() {
  DDATA_XSPLIT = UDATA_XSPLIT = 0;  // everything by default outside (arb. choice)
  DDATA_INDEST = UDATA_INDEST = -1; // since nothing inside

  DDATA_OUTDEST = PROC_NUM(PE_XIND, PE_YIND - 1);
  UDATA_OUTDEST = PROC_NUM(PE_XIND, PE_YIND + 1);
  
  IDATA_DEST = PROC_NUM(PE_XIND - 1, PE_YIND);
  ODATA_DEST = PROC_NUM(PE_XIND + 1, PE_YIND);
  
  TS_up_in = TS_up_out = TS_down_in = TS_down_out = false; // No twist-shifts

  /// Check if X is periodic
  if(periodicX) {
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
void BoutMesh::set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts) {
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
    throw BoutException("ERROR adding connection: y index %d or %d not on processor boundary\n", ypos1, ypos2);
  }

  /* check the x ranges are possible */
  if((xge != 0) && (xlt != MX)) {
    throw BoutException("ERROR adding connection(%d,%d,%d,%d): can only divide X domain in 2\n",
                            ypos1, ypos2, xge, xlt);
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
      if(UDATA_XSPLIT <= 0)
        UDATA_INDEST = UDATA_OUTDEST;
      UDATA_XSPLIT = xge;
      UDATA_OUTDEST = PROC_NUM(PE_XIND, ypedown);
      if(UDATA_XSPLIT <= 0)
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
      if(DDATA_XSPLIT <= 0)
        DDATA_INDEST = DDATA_OUTDEST;
      DDATA_XSPLIT = xge;
      DDATA_OUTDEST = PROC_NUM(PE_XIND, ypeup);
      if(DDATA_XSPLIT == 0)
	DDATA_INDEST = -1;

      TS_down_out = ts;

      output.write("=> This processor sending out down\n");
    }
  }
}

/// Add a divertor target or limiter
/*!
 * ypos is the y index which will become an upper target
 * ypos+1 will become a lower target.
 * Target created in the range xge <= x < xlt.
 */
void BoutMesh::add_target(int ypos, int xge, int xlt) {
  if(xlt <= xge)
    return;

  if((ypos < 0) || (ypos >= MY)) {
    output.write("WARNING adding target: poloidal index %d out of range\n", ypos);
    return;
  }
  
  int ypeup = YPROC(ypos);
  int ypedown = YPROC(ypos+1);
  if(ypeup == ypedown) {
    throw BoutException("Adding target at y=%d in middle of processor %d\n", ypos, ypeup);
  }

  output.write("Target at top of Y processor %d and bottom of %d in range %d <= x < %d\n", ypeup, ypedown, xge, xlt);

  // Convert X coordinates into local indices
  xge = XLOCAL(xge);
  xlt = XLOCAL(xlt);
  if(( xge >= ngx ) || ( xlt <= 0 )) {
    return; // Not in this x domain
  }

  if(MYPE == PROC_NUM(PE_XIND, ypeup)) {
    // Target on upper processor boundary
    if(xge <= MXG) {
      // Target on inside
      UDATA_XSPLIT = xlt;
      UDATA_INDEST = -1;
      if(xlt >= ngx)
        UDATA_OUTDEST = -1;
      output.write("=> This processor has target upper inner\n");
    }else {
      // Target on outside
      if(UDATA_XSPLIT <= 0)
        UDATA_INDEST = UDATA_OUTDEST;
      UDATA_XSPLIT = xge;
      UDATA_OUTDEST = -1;
      if(xge <= 0)
        UDATA_INDEST = -1;
      output.write("=> This processor has target upper outer\n");
    }
  }
  if(MYPE == PROC_NUM(PE_XIND, ypedown)) {
    // Target on upper processor boundary
    if(xge <= MXG) {
      // Target on inside
      DDATA_XSPLIT = xlt;
      DDATA_INDEST = -1;
      if(xlt >= ngx)
        DDATA_OUTDEST = -1;
      output.write("=> This processor has target lower inner\n");
    }else {
      // Target on outside
      if(DDATA_XSPLIT <= 0)
        DDATA_INDEST = DDATA_OUTDEST;
      DDATA_XSPLIT = xge;
      DDATA_OUTDEST = -1;
      if(xge <= 0)
        DDATA_INDEST = -1;
      output.write("=> This processor has target lower outer\n");
    }
  }
}

/****************************************************************
 *                MAIN TOPOLOGY FUNCTION
 ****************************************************************/

void BoutMesh::topology() {
  // Perform checks common to all topologies

  if (NPES != NXPE*NYPE) {
    throw BoutException("\tTopology error: npes=%d is not equal to NXPE*NYPE=%d\n",
                            NPES,NXPE*NYPE);
  }
  if(MYSUB * NYPE != MY) {
    throw BoutException("\tTopology error: MYSUB[%d] * NYPE[%d] != MY[%d]\n",MYSUB,NYPE,MY);
  }
  if(MXSUB * NXPE != MX) {
    throw BoutException("\tTopology error: MXSUB[%d] * NXPE[%d] != MX[%d]\n",MXSUB,NXPE,MX);
  }

  if((NXPE > 1) && (MXSUB < MXG)) {
    throw BoutException("\tERROR: Grid X size must be >= guard cell size\n");
  }
  if(MYSUB < MYG) {
    throw BoutException("\tERROR: Grid Y size must be >= guard cell size\n");
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
      throw BoutException("\tTopology error: Upper inner leg does not have integer number of processors\n");
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
    
    // Add target plates at the top
    add_target(ny_inner-1, 0, nx);
  }

  MYPE_IN_CORE = 0; // processor not in core
  if( (ixseps_inner > 0) && ( ((PE_YIND*MYSUB > jyseps1_1) && (PE_YIND*MYSUB <= jyseps2_1)) || ((PE_YIND*MYSUB > jyseps1_2) && (PE_YIND*MYSUB <= jyseps2_2)) ) ) {
    MYPE_IN_CORE = 1; /* processor is in the core */
  }
  
  if(DDATA_XSPLIT > ngx)
    DDATA_XSPLIT = ngx;
  if(UDATA_XSPLIT > ngx)
    UDATA_XSPLIT = ngx;

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

BoutMesh::CommHandle* BoutMesh::get_handle(int xlen, int ylen) {
  if(comm_list.empty()) {
    // Allocate a new CommHandle
    
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
    if(ch->xbufflen > 0) {
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
  
  ch->var_list.clear();

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
	  for(jz=0;jz < ngz-1;jz++)
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
	  for(jz=0;jz < ngz-1;jz++) {
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

/****************************************************************
 *                   Private data reading routines
 ****************************************************************/

/// Reads in a portion of the X-Y domain
int BoutMesh::readgrid_3dvar(GridDataSource *s, const char *name, 
	                     int yread, int ydest, int ysize, 
                             int xge, int xlt, BoutReal ***var) {
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

  int ncz = ngz-1;

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
      
      s->setGlobalOrigin(XGLOBAL(jx), yind);
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

  s->setGlobalOrigin();

  // free data
  delete[] zdata;
  delete[] fdata;
  
  return 0;
}

/// Copies a section of a 3D variable
void BoutMesh::cpy_3d_data(int yfrom, int yto, int xge, int xlt, BoutReal ***var) {
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
                             int xge, int xlt, BoutReal **var) {
  for(int i=xge;i!=xlt;i++) { // go through all the x indices 
    // Set the indices to read in this x position 
    s->setGlobalOrigin(XGLOBAL(i), yread);
    // Read in the block of data for this x value (C ordering)
    if(!s->fetch(&(var[i][ydest]), varname, 1, ysize))
      return 1;
  }
  
  s->setGlobalOrigin();
  
  return 0;
}

void BoutMesh::cpy_2d_data(int yfrom, int yto, int xge, int xlt, BoutReal **var) {
  msg_stack.push("cpy_2d_data(%d,%d,%d,%d)", yfrom, yto, xge, xlt);
  for(int i=xge;i<xlt;i++)
    var[i][yto] = var[i][yfrom];
  msg_stack.pop();
}

/****************************************************************
 *                 SURFACE ITERATION
 ****************************************************************/

bool BoutMesh::periodicY(int jx) const {
  return (XGLOBAL(jx) < ixseps_inner) && MYPE_IN_CORE;
}

bool BoutMesh::periodicY(int jx, BoutReal &ts) const {
  ts = 0.;
  if( (XGLOBAL(jx) < ixseps_inner) && MYPE_IN_CORE) {
    if(TwistShift)
      ts = ShiftAngle[jx];
    return true;
  }
  return false;
}

const Field2D BoutMesh::averageY(const Field2D &f) {
  static BoutReal *input = NULL, *result;
 
#ifdef CHECK
  msg_stack.push("averageY(Field2D)");
#endif
 
  if(input == NULL) {
    input = new BoutReal[ngx];
    result = new BoutReal[ngx];
  }
 
  BoutReal **fd = f.getData();
  
  // Average on this processor
  for(int x=0;x<ngx;x++) {
    input[x] = 0.;
    // Sum values, not including boundaries
    for(int y=ystart;y<=yend;y++) {
      input[x] += fd[x][y];
    }
    input[x] /= yend - ystart + 1;
  }

  Field2D r;
  r.allocate();
  BoutReal **rd = r.getData();

  int np;
  MPI_Comm_size(comm_inner, &np);
  
  if(np == 1) {
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        rd[x][y] = input[x];
  }else {
    MPI_Allreduce(input, result, ngx, MPI_DOUBLE, MPI_SUM, comm_inner);
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        rd[x][y] = result[x] / (BoutReal) np;
  }

#ifdef CHECK
  msg_stack.pop();
#endif
  
  return r;
}

const Field3D BoutMesh::averageY(const Field3D &f) {
  static BoutReal **input = NULL, **result;

#ifdef CHECK
  msg_stack.push("averageY(Field3D)");
#endif

  if(input == NULL) {
    input = rmatrix(ngx, ngz);
    result = rmatrix(ngx, ngz);
  }
  
  BoutReal ***fd = f.getData();
  
  // Average on this processor
  for(int x=0;x<ngx;x++)
    for(int z=0;z<ngz;z++) {
      input[x][z] = 0.;
      // Sum values, not including boundaries
      for(int y=ystart;y<=yend;y++) {
        input[x][z] += fd[x][y][z];
      }
      input[x][z] /= yend - ystart + 1;
    }
  
  Field3D r;
  r.allocate();
  BoutReal ***rd = r.getData();

  int np;
  MPI_Comm_size(comm_inner, &np);
  if(np > 1) {
    MPI_Allreduce(*input, *result, ngx*ngz, MPI_DOUBLE, MPI_SUM, comm_inner);
    
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
          rd[x][y][z] = result[x][z] / (BoutReal) np;
        }
  }else {
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
          rd[x][y][z] = input[x][z];
        }
  }

#ifdef CHECK
  msg_stack.pop();
#endif
  
  return r;
}

int BoutMesh::ySize(int xpos) const {
  int xglobal = XGLOBAL(xpos);
  int yglobal = YGLOBAL(MYG);
  
  if((xglobal < ixseps_lower) && ((yglobal <= jyseps1_1) || (yglobal > jyseps2_2))) {
    // Lower PF region
    return (jyseps1_1 + 1) + (ny - jyseps2_2);
    
  }else if((xglobal < ixseps_upper) && (yglobal > jyseps2_1) && (yglobal >= jyseps1_2)) {
    // Upper PF region
    return jyseps1_2 - jyseps2_1;
    
  }else if(xglobal < ixseps_inner) {
    // Core
    return (jyseps2_1 - jyseps1_1) + (jyseps2_2 - jyseps1_2);
    
  }else if(jyseps2_1 == jyseps1_2) {
    // Single null, so in the SOL
    return ny;
    
  }else if((xglobal >= ixseps_inner) && (xglobal < ixseps_outer)) {
    // Intermediate SOL in DND
    
    if(ixseps_lower < ixseps_upper) {
      // Connects to lower divertor
      return (jyseps2_1 + 1) + (ny - jyseps1_2);
    }else {
      // Connects to upper divertor
      return jyseps2_2 - jyseps1_1;
    }
  }else if(yglobal < ny_inner) {
    // Inner SOL
    return ny_inner;
  }
  // Outer SOL
  return ny - ny_inner;
}

MPI_Comm BoutMesh::getYcomm(int xpos) const {
  int xglobal = XGLOBAL(xpos);
  
  if(xglobal < ixseps_inner) {
    return comm_inner;
  }else if(xglobal < ixseps_outer) {
    return comm_middle;
  }
  return comm_outer;
}

/****************************************************************
 *                 Range iteration
 ****************************************************************/

const RangeIterator BoutMesh::iterateBndryLowerY() const {
  
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

  return RangeIterator(xs, xe);
}

const RangeIterator BoutMesh::iterateBndryUpperY() const {
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

  return RangeIterator(xs, xe);
}


vector<BoundaryRegion*> BoutMesh::getBoundaries() {
  return boundary;
}

const Field3D BoutMesh::smoothSeparatrix(const Field3D &f) {
  Field3D result(f);
  if((ixseps_inner > 0) && (ixseps_inner < nx-1)) {
    result.allocate();
    if(XPROC(ixseps_inner) == PE_XIND) {
      int x = XLOCAL(ixseps_inner);
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
	  result[x][y][z] = 0.5*(f[x][y][z] + f[x-1][y][z]);
          //result[x][y][z] = f[x-1][y][z] - 2.*f[x][y][z] + f[x+1][y][z];
        }
    }
    if(XPROC(ixseps_inner-1) == PE_XIND) {
      int x = XLOCAL(ixseps_inner-1);
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
	  result[x][y][z] = 0.5*(f[x][y][z] + f[x+1][y][z]);
          //result[x][y][z] = f[x-1][y][z] - 2.*f[x][y][z] + f[x+1][y][z];
        }
    }
  }
  if((ixseps_outer > 0) && (ixseps_outer < nx-1) && (ixseps_outer != ixseps_inner)) {
    result.allocate();
    if(XPROC(ixseps_outer) == PE_XIND) {
      int x = XLOCAL(ixseps_outer);
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
          result[x][y][z] = 0.5*(f[x][y][z] + f[x-1][y][z]); //f[x-1][y][z] - 2.*f[x][y][z] + f[x+1][y][z];
        }
    }
    if(XPROC(ixseps_outer-1) == PE_XIND) {
      int x = XLOCAL(ixseps_outer-1);
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
          result[x][y][z] = 0.5*(f[x][y][z] + f[x+1][y][z]); //f[x-1][y][z] - 2.*f[x][y][z] + f[x+1][y][z];
        }
    }
  }
  return result;
}

BoutReal BoutMesh::GlobalX(int jx) const {
  if(symmetricGlobalX) {
    // Symmetric X index, mainly for reconnection studies
    return ((BoutReal) XGLOBAL(jx)) / ((BoutReal) nx-1);
  }
  return ((BoutReal) XGLOBAL(jx)) / ((BoutReal) MX);
}

BoutReal BoutMesh::GlobalY(int jy) const {
  int ly = YGLOBAL(jy); // global poloidal index across subdomains
  int nycore = (jyseps2_1 - jyseps1_1) + (jyseps2_2 - jyseps1_2);

  if(MYPE_IN_CORE) {
    // Turn ly into an index over the core cells only
    if(ly <= jyseps2_1) {
      ly -= jyseps1_1+1;
    }else
      ly -= jyseps1_1+1 + (jyseps1_2 - jyseps2_1);
  }else {
    // Not in core. Need to get the last "core" value
    if(ly <= jyseps1_1) {
      // Inner lower leg
      ly = 0;
    }else if((ly > jyseps2_1) && (ly <= jyseps1_2)) {
      // Upper legs
      ly = jyseps2_1 - jyseps1_1;
    }else if(ly > jyseps2_2) {
      // Outer lower leg
      ly = nycore;
    }
  }
  
  //output.write("GlobalY: %d, %d, %d, %d -> %e\n", jy, YGLOBAL(jy), ly, nycore, ((BoutReal) ly) / ((BoutReal) nycore));

  return ((BoutReal) ly) / ((BoutReal) nycore);
}

void BoutMesh::outputVars(Datafile &file) {
  file.add(MXSUB, "MXSUB", 0);
  file.add(MYSUB, "MYSUB", 0);
  file.add(MXG,   "MXG",   0);
  file.add(MYG,   "MYG",   0);
  file.add(ngz,   "MZ",    0);
  file.add(NXPE,  "NXPE",  0);
  file.add(NYPE,  "NYPE",  0);
  file.add(ZMAX,  "ZMAX",  0);
  file.add(ZMIN,  "ZMIN",  0);
  
  file.add(dx,    "dx",    0);
  file.add(dy,    "dy",    0);
  file.add(dz,    "dz",    0);
  
  file.add(g11,   "g11",   0);
  file.add(g22,   "g22",   0);
  file.add(g33,   "g33",   0);
  file.add(g12,   "g12",   0);
  file.add(g13,   "g13",   0);
  file.add(g23,   "g23",   0);
  
  file.add(g_11,  "g_11",  0);
  file.add(g_22,  "g_22",  0);
  file.add(g_33,  "g_33",  0);
  file.add(g_12,  "g_12",  0);
  file.add(g_13,  "g_13",  0);
  file.add(g_23,  "g_23",  0);
  
  file.add(J,     "J",     0);
}


//================================================================
// Poloidal lowpass filter for n=0 mode, keeping 0<=m<=mmax
// Developed by T. Rhee and S. S. Kim
//================================================================

void BoutMesh::slice_r_y(BoutReal *fori, BoutReal * fxy, int ystart, int ncy)
{
  int i,j;
  for(i=0;i<ncy;i++)
      fxy[i]=fori[i+ystart];
}

void BoutMesh::get_ri( dcomplex * ayn, int ncy, BoutReal * ayn_Real, BoutReal * ayn_Imag)
{
  for(int i=0;i<ncy;i++){
    ayn_Real[i]=ayn[i].Real();
    ayn_Imag[i]=ayn[i].Imag();
  }
}

void BoutMesh::set_ri( dcomplex * ayn, int ncy, BoutReal * ayn_Real, BoutReal * ayn_Imag)
{
  for(int i=0;i<ncy;i++){
    ayn[i]=dcomplex(ayn_Real[i],ayn_Imag[i]);
  }
}

// Lowpass filter for n=0 mode, keeping poloidal mode number 0<=m<=mmax
const Field2D BoutMesh::lowPass_poloidal(const Field2D &var,int mmax)
{
  Field2D result;
  static BoutReal *f1d = (BoutReal *) NULL;
  static dcomplex *aynall = (dcomplex*)NULL;
  static BoutReal *aynall_Real = (BoutReal *) NULL;
  static BoutReal *aynall_Imag = (BoutReal *) NULL;
  static dcomplex *ayn = (dcomplex*) NULL;
  static BoutReal *aynReal = (BoutReal *) NULL;
  static BoutReal *aynImag = (BoutReal *) NULL;

  int ncx, ncy;
  int jx, jy;
  int ncyall,ncyall2;   // nype is number of processors in the Y dir.
  int t1; //return value for setData of Field2D
  int i,j,k;
  int mmax1;

  mmax1 = mmax+1;
  ncx = ngx;
  ncy = yend - ystart + 1;
  ncyall = ncy*NYPE;
  ncyall2 = ncyall/2;

  result.allocate();  //initialize
  if(f1d == (BoutReal*) NULL)
     f1d = new BoutReal[ncy];  
  if(ayn == (dcomplex*) NULL)
     ayn = new dcomplex[ncy];  
  if(aynall == (dcomplex*) NULL)
     aynall = new dcomplex[ncyall];  
  if(aynall_Real == (BoutReal*) NULL)
     aynall_Real = new BoutReal[ncyall];  
  if(aynall_Imag == (BoutReal*) NULL)
     aynall_Imag = new BoutReal[ncyall];  
  if(aynReal == (BoutReal*) NULL)
     aynReal = new BoutReal[ncy];  
  if(aynImag == (BoutReal*) NULL)
     aynImag = new BoutReal[ncy];  

  for(jx=0;jx<ncx;jx++){ //start x
     //saving the real 2D data
     slice_r_y(*(var.getData()+jx),f1d,ystart,ncy); // 2d -> 1d

     for(jy=0;jy<ncy;jy++)
       ayn[jy]=dcomplex(f1d[jy],0.);

     // allgather from ayn to aynall
     get_ri(ayn,ncy,aynReal,aynImag); 
     MPI_Allgather(aynReal, ncy, MPI_DOUBLE, aynall_Real, 
                      ncy, MPI_DOUBLE, comm_inner);
     MPI_Allgather(aynImag, ncy, MPI_DOUBLE, aynall_Imag, 
                      ncy, MPI_DOUBLE, comm_inner);
     set_ri(aynall,ncyall,aynall_Real,aynall_Imag); 

     // FFT in y over extended domain
     cfft(aynall, ncyall, -1); 

     // lowpass filter
     for(jy=mmax1;jy<=(ncyall-mmax1);jy++)
       aynall[jy]= dcomplex(0.,0.); 

     // inverse FFT in y over extended domain
     cfft(aynall, ncyall, 1);

     // scatter data
     get_ri(aynall,ncyall,aynall_Real,aynall_Imag); 
     MPI_Scatter(aynall_Real, ncy, MPI_DOUBLE, aynReal, 
                      ncy, MPI_DOUBLE, 0,comm_inner);
     MPI_Scatter(aynall_Imag, ncy, MPI_DOUBLE, aynImag, 
                      ncy, MPI_DOUBLE, 0,comm_inner);
     set_ri(ayn,ncy,aynReal,aynImag);

     for(jy=0;jy<ncy;jy++)
      f1d[jy]=ayn[jy].Real();

    for(jy=0;jy<ncy;jy++)
      t1=result.setData(jx,jy+ystart,1,f1d+jy); 
  }//end of x

  return result; 
} 

/*================================================================
// Volume integral of Field2D variable
// Developed by T. Rhee and S. S. Kim
//================================================================*/

//BoutReal Vol_Average(const Field2D &var);
BoutReal BoutMesh::Average_XY(const Field2D &var)
{
  Field2D result;
  BoutReal Vol_Loc, Vol_Glb;
  int i;
  result.allocate();  //initialize
  result=averageY(var);

  Vol_Loc = 0.;
  Vol_Glb = 0.;

  for (i=xstart;i<=xend;i++)
      Vol_Loc +=  result[i][0];

  MPI_Allreduce(&Vol_Loc,&Vol_Glb,1,MPI_DOUBLE,MPI_SUM,comm_x);
  Vol_Glb /= (BoutReal)(nx - 2*MXG);

  return Vol_Glb;
}
BoutReal BoutMesh::Vol_Integral(const Field2D &var)
{
  Field2D result;
  BoutReal Int_Glb;
  result.allocate();  //initialize
  result = J * var * dx * dy;

  Int_Glb = 0.;
  Int_Glb = Average_XY(result);
  Int_Glb *= (BoutReal) ( (nx - 2*MXG)*ny )*PI * 2.;

  return Int_Glb;
}


const Field3D BoutMesh::Switch_YZ(const Field3D &var)
{
  static BoutReal **ayz = (BoutReal **) NULL;
  static BoutReal **ayz_all = (BoutReal **) NULL;
  Field3D  result;
  int ncy, ncy_all,ncz;
  int i,j,ix;
  ncy = yend - ystart + 1;
  ncy_all = MY;
  ncz = ngz-1 ;

  if(MY != ngz-1){
    throw new BoutException("Y and Z dimension is not same in Switch_YZ code"); }

  //memory allocation
  result.allocate();
  if(ayz == (BoutReal**) NULL)
    ayz=rmatrix(ncy,ncz);
  if(ayz_all == (BoutReal**) NULL)
    ayz_all=rmatrix(ncy_all,ncz);

  for(ix=xstart;ix<=xend;ix++){
    //Field 3D to rmatrix of local
    for (i=0;i<ncy;i++)
      for(j = 0;j<ncz;j++)
        ayz[i][j]=var[ix][i+ystart][j];

    //Collect to rmatrix of global from local
    MPI_Allgather(ayz[0],ncy*ncz,MPI_DOUBLE,ayz_all[0],ncy*ncz,
                MPI_DOUBLE,comm_inner);

    //Y 2 Z switch
    for(i=ystart;i<=yend;i++)
      for(j=0;j<ncz;j++)
        result[ix][i][j]=ayz_all[j][YGLOBAL(i)];

    //boundary at ngz
    for(i=0;i<ncy;i++)
        result[ix][i+ystart][ncz]=result[ix][i+ystart][0];
  }

  return result;
}

