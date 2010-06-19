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

#endif // __GRID_H__
