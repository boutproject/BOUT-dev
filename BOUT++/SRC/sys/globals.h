/************************************************************************//**
 * \brief Global variables for BOUT++
 * 
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

#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include "mpi.h"

#include "bout_types.h"
#include "field2d.h"
#include "options.h"
#include "output.h"
#include "msg_stack.h" 

#include "datafile.h"
#include "grid.h"
#include "mesh.h"

#ifndef GLOBALORIGIN
#define GLOBAL extern
#define SETTING(name, val) extern name;
#else
#define GLOBAL
#define SETTING(name, val) name = val;
#endif

GLOBAL Mesh *mesh; ///< The mesh object

const real PI = 3.141592653589793;
const real TWOPI = 6.2831853071795;

// Grid sizes/indexes set in grid.cpp
GLOBAL int nx, ny;        ///< Size of the grid in the input file
GLOBAL int MX, MY;        ///< size of the grid excluding boundary regions
GLOBAL int MYSUB, MXSUB;  ///< Size of the grid on this processor
GLOBAL int ngx, ngy, ngz; ///< Total domain size on this processor including guard/boundary cells
GLOBAL int ncx, ncy, ncz;

GLOBAL int xstart, xend, jstart, jend; // local index range

GLOBAL int NPES; ///< Number of processors (bout++.cpp)
GLOBAL int MYPE; ///< Rank of this processor (bout++.cpp)
GLOBAL int PE_YIND; ///< Y index of this processor (bout++.cpp)
GLOBAL int PE_XIND; ///< X index (bout++.cpp)

GLOBAL int NXPE; ///< Number of processors in X direction (bout++.cpp)
GLOBAL int NYPE; ///< Number in Y direction (bout++.cpp)

// Communicators
GLOBAL MPI_Comm comm_x, comm_y;

// separatrix indices (read in grid.cpp)
GLOBAL int ixseps1, ixseps2, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2; 
GLOBAL int ixseps_inner, ixseps_outer, ixseps_upper, ixseps_lower;
GLOBAL int ny_inner;

// Twist-shift switches (topology.cpp)
GLOBAL bool TS_up_in, TS_up_out, TS_down_in, TS_down_out;

// Communication parameters calculated by topology
GLOBAL int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
GLOBAL int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
GLOBAL int IDATA_DEST, ODATA_DEST; // X inner and outer destinations

GLOBAL int MYPE_IN_CORE; // 1 if processor in core (topology.cpp)

///////////////// TWIST-SHIFT CONDITION ///////////////////////

#ifndef METRIC3D
GLOBAL real *ShiftAngle;  // angle a field-line moves in z for a y circuit (radians)
#else
GLOBAL FieldPerp ShiftAngle; // X-Z (shift now depends on z)
#endif

///////////////// SHEARED X DERIVATIVES ///////////////////////

#ifndef METRIC3D
GLOBAL Field2D zShift; // Z shift for each point (radians)
GLOBAL Field2D ShiftTorsion; // d <pitch angle> / dx. Needed for vector differentials (Curl)
GLOBAL Field2D IntShiftTorsion; // Integrated shear (I in BOUT notation)
#else
// 3D metrics
GLOBAL Field3D zShift; 
GLOBAL Field3D ShiftTorsion;
#endif

///////////////// DIFFERENCING QUANTITIES /////////////////////

// These used for differential operators 
#ifndef METRIC3D
GLOBAL Field2D dx, dy;      // Read in grid.cpp
GLOBAL Field2D d2x, d2y;    // 2nd-order correction for non-uniform meshes
GLOBAL real zlength, dz;    // Derived from options in grid.cpp (in radians)
#else
GLOBAL Field3D dx, dy, dz;
GLOBAL real zlength;
#endif

////////////////// DIFFERENTIAL GEOMETRY //////////////////////
#ifndef METRIC3D // Normal 2D metric case

// Contravariant metric tensor (g^{ij})
GLOBAL Field2D g11, g22, g33, g12, g13, g23; // These are read in grid.cpp

// rest of these quantities are derived from the metric tensor
// in geometry.cpp

GLOBAL Field2D J; // Jacobian
GLOBAL Field2D Bxy; // Magnitude of B = nabla z times nabla x

// Covariant metric tensor
GLOBAL Field2D g_11, g_22, g_33, g_12, g_13, g_23;

// Christoffel symbol of the second kind (connection coefficients)
GLOBAL Field2D G1_11, G1_22, G1_33, G1_12, G1_13;
GLOBAL Field2D G2_11, G2_22, G2_33, G2_12, G2_23;
GLOBAL Field2D G3_11, G3_22, G3_33, G3_13, G3_23;

GLOBAL Field2D G1, G2, G3;

#else // 3D metric case
// Contravariant metric tensor (g^{ij})
GLOBAL Field3D g11, g22, g33, g12, g13, g23; // These are read in grid.cpp

// rest of these quantities are derived from the metric tensor
// in geometry.cpp

GLOBAL Field3D J; // Jacobian
GLOBAL Field3D Bxy; // Magnitude of B = nabla z times nabla x

// Covariant metric tensor
GLOBAL Field3D g_11, g_22, g_33, g_12, g_13, g_23;

// Christoffel symbol of the second kind (connection coefficients)
GLOBAL Field3D G1_11, G1_22, G1_33, G1_12, G1_13;
GLOBAL Field3D G2_11, G2_22, G2_33, G2_12, G2_23;
GLOBAL Field3D G3_11, G3_22, G3_33, G3_13, G3_23;

GLOBAL Field3D G1, G2, G3;
#endif

///////////////////////////////////////////////////////////////

/// Options file object
GLOBAL OptionFile options;

/// Define for reading options which passes the variable name
#define OPTION(var, def) options.get(#var, var, def)

/// Output object
GLOBAL Output output;

/// Dump file object
GLOBAL Datafile dump;

/// Status message stack. Used for debugging messages
GLOBAL MsgStack msg_stack;

/// Define for reading a variable from the grid
#define GRID_LOAD(var) mesh->get(var, #var)

// Settings (read from file, bout++.cpp)

GLOBAL bool restarting;   // Specifies whether code is restarting

GLOBAL bool ShiftXderivs; // Use shifted X derivatives
GLOBAL bool IncIntShear;  // Include integrated shear (if shifting X)
GLOBAL bool TwistShift;   // Use a twist-shift condition in core?
GLOBAL int  ShiftOrder;   // Order of shifted X derivative interpolation
GLOBAL int  TwistOrder;   // Order of twist-shift interpolation
GLOBAL int  MZ;           // Number of points in the Z direction
GLOBAL int  zperiod;      // Number of z domains in 2 pi
GLOBAL real ZMIN;
GLOBAL real ZMAX;
GLOBAL int  MXG;
GLOBAL int  MYG;
GLOBAL bool BoundaryOnCell;  ///< Boundary is on the last "real" point. Otherwise between points (old method)

GLOBAL bool StaggerGrids;    ///< Enable staggered grids (Centre, Lower). Otherwise all vars are cell centred (default).

// Timing information
GLOBAL real wtime_invert; //< Time spent performing inversions

GLOBAL bool non_uniform;  // Use corrections for non-uniform meshes

// Error handling (bout++.cpp)
void bout_error();
void bout_error(const char *str);

#undef GLOBAL
#undef SETTING

#endif // __GLOBALS_H__
