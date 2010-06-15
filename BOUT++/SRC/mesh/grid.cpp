/*!
 * \file grid.cpp
 * \brief Functions to read variables from a grid file
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
 */

#include "mpi.h"

#include "globals.h"
#include "grid.h"
#include "topology.h"
#include "utils.h"
#include "geometry.h"

#include "dataformat.h" // Abstract data class

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "fft.h" // For reading 3D variables
#include "dcomplex.h"

#include "bout_types.h"

/*******************************************************************************
 * GridFile class
 * 
 * This is a fairly thin wrapper around the Datafile class. Only reading
 * routines are required though. Takes care of opening and closing file
 * then just passes read requests through.
 *******************************************************************************/

/// Default constructor
GridFile::GridFile()
{
  file = NULL; ///< Signal that the file is invalid
  isOpen = false;
}

/// Constructor passing a file format and name
GridFile::GridFile(DataFormat *format, const char *gridfilename)
{
  setFile(format, gridfilename);
}

void GridFile::setFile(DataFormat *format, const char *gridfilename)
{
  file = format;
  filename = string(gridfilename);
  isOpen = false;
}

bool GridFile::hasVar(const char *name)
{
  if(file == NULL)
    return false;

  /// Try to get the size of the variable
  vector<int> s = getSize(name);
  
  /// Test if the variable has zero size
  return s.size() != 0;
}

/// Return the size of variable
vector<int> GridFile::getSize(const char *name)
{
  vector<int> s;
  
  if(file == NULL)
    return s;
  
  if(!isOpen) {
    // File not open, so need to open
    if(!file->openr(filename))
      return s;
  }

  s = file->getSize(name);
  
  if(!isOpen)
    file->close(); // If was closed before then close again
  
  return s;
}

bool GridFile::setOrigin(int x, int y, int z)
{
  if(file == NULL)
    return false;
  
  return file->setOrigin(x,y,z);
}

bool GridFile::fetch(int *var, const char *name, int lx, int ly, int lz)
{
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

bool GridFile::fetch(int *var, const string &name, int lx, int ly, int lz)
{
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

bool GridFile::fetch(real *var, const char *name, int lx, int ly, int lz)
{
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

bool GridFile::fetch(real *var, const string &name, int lx, int ly, int lz)
{
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

void GridFile::open(const char *name)
{
  file->openr(filename);
  if(!file->is_valid()) {
    output << "\tERROR: Could not open file " << filename << endl;
  }
  isOpen = true;
}

void GridFile::close()
{
  file->close();
  isOpen = false;
}

/*******************************************************************************
 * GridData class
 *******************************************************************************/

GridData::GridData()
{
  
}

GridData::GridData(GridDataSource *source)
{
  addSource(source);
}

bool GridData::addSource(GridDataSource *source)
{
  source_list.push_front(source);
  return true;
}

bool GridData::loadTopology()
{
#ifdef CHECK
  msg_stack.push("grid::loadTopology");
#endif
  
  output << "Loading topology" << endl;
  
  //////////////
  // Grid sizes
  
  if(get(nx, "nx"))
    return false;
  
  if(get(ny, "ny"))
    return false;

  output << "\tGrid size: " << nx << " by " << ny << endl;

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
  
  /// Number of grid cells is ng* = M*SUB + guard/boundary cells
  ngx = MXSUB + 2*MXG;
  ngy = MYSUB + 2*MYG;
  ngz = MZ;
  
  ncx = ngx - 1;
  ncy = ngy - 1;
  ncz = ngz - 1;

  // Set local index ranges
  
  jstart = MYG;
  jend = MYG + MYSUB - 1;

  xstart = MXG;
  xend = MXG + MXSUB - 1;
  
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
  dz = zlength/ncz;

  ///////////////// DIFFERENTIAL GEOMETRY /////////////////
  
  // Diagonal components of metric tensor g^{ij} (default to 1)
  get(g11, "g11", 1.0);
  get(g22, "g22", 1.0);
  get(g33, "g33", 1.0);
  
  // Off-diagonal elements. Default to 0
  get(g12, "g12", 0.0);
  get(g13, "g13", 0.0);
  get(g23, "g23", 0.0);
  
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

  ShiftAngle  = rvector(ngx);

  // Try to read the shift angle from the grid file
  // NOTE: All processors should know the twist-shift angle (for invert_parderiv)
  GridDataSource* s = find_source("ShiftAngle");
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
    output.write("\tWARNING: Twist-shift angle 'ShiftAngle' not found. Setting to zero\n");
    for(int i=0;i<ngx;i++)
      ShiftAngle[i] = 0.0;
  }

  /// Call geometry routines to set the rest of the coefficients
  geometry(1);    // Calculate metrics, christoffel symbols
  
  output.write("\tdone\n");
#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

/// Get an integer
int GridData::get(int &ival, const char *name)
{
#ifdef CHECK
  msg_stack.push("Loading integer: grid::get(int, %s)", name);
#endif

  GridDataSource* s = find_source(name);
  if(s == NULL) {
#ifdef CHECK
    msg_stack.pop();
#endif
    return 1;
  }
  
  s->open(name);
  bool success = s->fetch(&ival, name);
  s->close();
  
#ifdef CHECK
  msg_stack.pop();
#endif

  if(!success) {
    return 2;
  }
  return 0;
}

/// A real number
int GridData::get(real &rval, const char *name)
{
  GridDataSource* s = find_source(name);
  if(s == NULL)
    return 1;
  
  s->open(name);
  bool success = s->fetch(&rval, name);
  s->close();
  
  if(!success)
    return 2;
  return 0;
}

int GridData::get(Field2D &var, const char *name, real def)
{
  if(name == NULL)
    return 1;
  
#ifdef CHECK
  msg_stack.push("Loading 2D field: grid::get(Field2D, %s)", name);
#endif
  
  GridDataSource *s = find_source(name);
  if(s == NULL) {
    output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name, def);
    var = def;
#ifdef CHECK
    msg_stack.pop();
#endif
    return 2;
  }
  
  real **data;
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
    msg_stack.pop();
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
      msg_stack.pop();
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
      msg_stack.pop();
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
      msg_stack.pop();
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
	msg_stack.pop();
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

   
  // Close signal
  s->close();
  
#ifdef CHECK
  msg_stack.pop();
#endif
  
  return 0;
}

int GridData::get(Field2D &var, const string &name, real def)
{
  return get(var, name.c_str());
}

/// Load a 3D variable from the grid file
/*!
  Data stored as toroidal FFTs in real space at each X-Y point.
  In toroidal direction, array must have an odd number of points.
  Format is:

  DC, r1,i1, r2,i2, ... , rn,in

  with the real and imaginary parts of each (positive) frequency
  up to the nyquist frequency.
 */
int GridData::get(Field3D &var, const char *name)
{
  if(name == NULL)
    return 1;
  
#ifdef CHECK
  msg_stack.push("Loading 3D field: grid::get(Field3D, %s)", name);
#endif
  
  GridDataSource *s = find_source(name);
  if(s == NULL) {
    output.write("\tWARNING: Could not read '%s' from grid. Setting to zero\n", name);
    var = 0.0;
#ifdef CHECK
    msg_stack.pop();
#endif
    return 2;
  }

  real ***data;
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
    msg_stack.pop();
#endif
  
  return 0;
}

int GridData::get(Field3D &var, const string &name)
{
  return get(var, name.c_str());
}

int GridData::get(Vector2D &var, const char *name)
{
  return get(var, string(name));
}

int GridData::get(Vector3D &var, const char *name)
{
  return get(var, string(name));
}

int GridData::get(Vector2D &var, const string &name)
{
#ifdef CHECK
  msg_stack.push("Loading 2D vector: grid::get(Vector2D, %s)", name.c_str());
#endif

  if(var.covariant) {
    output << "\tReading covariant vector " << name << endl;
    
    get(var.x, name+"_x");
    get(var.y, name+"_y");
    get(var.z, name+"_z");
    
  }else {
    output << "\tReading contravariant vector " << name << endl;
    
    get(var.x, name+"x");
    get(var.y, name+"y");
    get(var.z, name+"z");
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif  

  return 0;
}

int GridData::get(Vector3D &var, const string &name)
{
#ifdef CHECK
  msg_stack.push("Loading 3D vector: grid::get(Vector3D, %s)", name.c_str());
#endif

  if(var.covariant) {
    output << "\tReading covariant vector " << name << endl;
    
    get(var.x, name+"_x");
    get(var.y, name+"_y");
    get(var.z, name+"_z");
    
  }else {
    output << "\tReading contravariant vector " << name << endl;
    
    get(var.x, name+"x");
    get(var.y, name+"y");
    get(var.z, name+"z");
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif  

  return 0;
}

/*******************************************************************************
 * GridDat private functions
 *******************************************************************************/

GridDataSource* GridData::find_source(const char *name)
{
  for(std::list<GridDataSource*>::iterator it = source_list.begin(); 
      it != source_list.end(); it++) {
    
    // Query this source
    if((*it) != NULL)
      if((*it)->hasVar(name)) {
	return *it;
      }
  }
  return NULL;
}

/// Reads in a portion of the X-Y domain
int GridData::readgrid_3dvar(GridDataSource *s, const char *name, 
	                     int yread, int ydest, int ysize, 
                             int xge, int xlt, real ***var)
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
  real* zdata = new real[size[2]];

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
	real kwave;
	
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
void GridData::cpy_3d_data(int yfrom, int yto, int xge, int xlt, real ***var)
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
int GridData::readgrid_2dvar(GridDataSource *s, const char *varname, 
                             int yread, int ydest, int ysize, 
                             int xge, int xlt, real **var)
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

void GridData::cpy_2d_data(int yfrom, int yto, int xge, int xlt, real **var)
{
  int i;
  for(i=xge;i!=xlt;i++)
    var[i][yto] = var[i][yfrom];
}

/*******************************************************************************
 * OBSOLETE INTERFACE FUNCTIONS FOR BACKWARDS COMPATIBILITY
 *******************************************************************************/

/// Main grid function
/*!
  Reads the variables needed by the core BOUT++ code
  then calls topology
*/
bool grid_read(DataFormat *format, const char *gridfilename)
{
  /// Add a grid file source
  grid.addSource(new GridFile(format, gridfilename));
  /// Read the basic topology data
  return grid.loadTopology();
}

int grid_load(real &r, const char *name)
{
  return grid.get(r, name);
}

int grid_load(int &i, const char *name)
{
  return grid.get(i, name);
}

/// Loads a 2D field from grid file
int grid_load2d(Field2D &var, const char *name)
{
  return grid.get(var, name);
}

/// Loads a 2D vector from grid file
int grid_load2d(Vector2D &var, const char *name)
{
  return grid_load2d(var, string(name));
}

/// Loads a 2D vector from grid file
int grid_load2d(Vector2D &var, const string &name)
{
  return grid.get(var, name);
}

bool grid_load3d(Field3D &var, const char *name)
{
  
  if(grid.get(var, name))
    return false;
  return true;
}

bool grid_load3d(Vector3D &var, const char *name)
{
  return grid_load3d(var, string(name));
}

/// Loads a 3D vector from grid file
bool grid_load3d(Vector3D &var, const string &name)
{
  if(grid.get(var, name))
    return false;
  return true;
}

