/*******************************************************************************
 * Functions to read variables from a uedge.grd.pdb file
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
 *******************************************************************************/

#ifndef __GRID_H__
#define __GRID_H__

#include "field2d.h"
#include "vector2d.h"

#include "dataformat.h"
#include "bout_types.h"

#include <list>

/// Interface class to serve grid data
/*!
 * Provides a generic interface for sources of
 * equilibrium data. 
 * Could be used to simplify interfacing between BOUT++ and other codes
 */
class GridDataSource {
 public:
  virtual ~GridDataSource() { }
  
  virtual bool hasVar(const char *name) = 0; ///< Test if source can supply a variable
  
  virtual vector<int> getSize(const char *name) = 0; ///< Get size of the variable

  /// Set the (x,y,z) origin for all subsequent calls
  virtual bool setOrigin(int x = 0, int y = 0, int z = 0) = 0;
  
  /// Get data from the source
  virtual bool fetch(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool fetch(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool fetch(real *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool fetch(real *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  
  // Optional functions for efficiency

  /// Called before fetching a variable
  virtual void open(const char *name = NULL) {  }
  /// Called after fetching variables
  virtual void close() { }
};

/// Interface to grid data in a file
/*!
 * This is a thin wrapper around a DataFormat object. Only needs to implement 
 * reading routines.
 */
class GridFile : public GridDataSource {
 public:
  GridFile();
  GridFile(DataFormat *format, const char *gridfilename);
  
  void setFile(DataFormat *format, const char *gridfilename);
  
  virtual bool hasVar(const char *name);
  
  virtual vector<int> getSize(const char *name);

  virtual bool setOrigin(int x = 0, int y = 0, int z = 0);

  virtual bool fetch(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  virtual bool fetch(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0);
  virtual bool fetch(real *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  virtual bool fetch(real *var, const string &name, int lx = 1, int ly = 0, int lz = 0);
  
  virtual void open(const char *name = NULL);
  virtual void close();
 private:
  DataFormat *file;
  string filename;
  
  bool isOpen;
};

/// Class to handle equilibrium grid quantities
/*!
 * Physics code requests data from this class.
 * Depending on the variable, this then finds the
 * correct source.
 */
class GridData {
 public:
  GridData();
  GridData(GridDataSource *source); ///< Set a default input source
  
  /// Add a data source
  bool addSource(GridDataSource *source);
  
  /// Topology initialisation routine
  bool loadTopology();
  
  // Get routines to request data
  int get(int &ival, const char *name); ///< Get an integer
  int get(real &rval, const char *name); ///< Get a real number
  
  int get(Field2D &var, const char *name, real def=0.0);
  int get(Field2D &var, const string &name, real def=0.0);
  int get(Field3D &var, const char *name);
  int get(Field3D &var, const string &name);
  
  int get(Vector2D &var, const char *name);
  int get(Vector3D &var, const char *name);
  int get(Vector2D &var, const string &name);
  int get(Vector3D &var, const string &name);
  
 private:
  
  std::list<GridDataSource*> source_list; ///< List of sources
  
  GridDataSource *find_source(const char *name);
  
  /// Read in a portion of the X-Y domain
  int readgrid_3dvar(GridDataSource *s, const char *name, 
	             int yread, int ydest, int ysize, 
                     int xge, int xlt, real ***var);
  
  /// Copy a section of a 3D variable
  void cpy_3d_data(int yfrom, int yto, int xge, int xlt, real ***var);

  int readgrid_2dvar(GridDataSource *s, const char *varname, 
                     int yread, int ydest, int ysize, 
                     int xge, int xlt, real **var);
  void cpy_2d_data(int yfrom, int yto, int xge, int xlt, real **var);
  
};

bool grid_read(DataFormat *format, const char *gridfilename);

int grid_load(real &r, const char *name);
int grid_load(int &i, const char *name);

int grid_load2d(Field2D &var, const char *name);
int grid_load2d(Vector2D &var, const char *name);
int grid_load2d(Vector2D &var, const string &name);

bool grid_load3d(Field3D &var, const char *name);
bool grid_load3d(Vector3D &var, const char *name);
bool grid_load3d(Vector3D &var, const string &name);

#endif // __GRID_H__
