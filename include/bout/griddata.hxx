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


class GridDataSource;

#ifndef __GRIDDATA_H__
#define __GRIDDATA_H__

#include "options.hxx"

#include "dataformat.hxx"
#include "bout_types.hxx"

#include "mesh.hxx"

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
  
  virtual bool isValid() {return true;}  ///< Check if the data is valid

  virtual bool hasVar(const string &name) = 0; ///< Test if source can supply a variable
  
  virtual vector<int> getSize(const string &name) = 0; ///< Get size of the variable

  /// Set the (x,y,z) origin for all subsequent calls
  virtual bool setGlobalOrigin(int x = 0, int y = 0, int z = 0) = 0;

  /// Specify a mesh. Required for GridFromOptions
  virtual bool setMesh(Mesh *m) {return true;}
  
  /// Get data from the source
  virtual bool fetch(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool fetch(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool fetch(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool fetch(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  
  // Optional functions for efficiency

  /// Called before fetching a variable
  virtual void open(const char *name = NULL) {  }
  
  /// Called after fetching variables
  virtual void close() { }
  
  void open(const string &name) { open(name.c_str()); }
};

/// Interface to grid data in a file
/*!
 * This is a thin wrapper around a DataFormat object. Only needs to implement 
 * reading routines.
 */
class GridFile : public GridDataSource {
 public:
  GridFile();
  ~GridFile();
  GridFile(DataFormat *format, const char *gridfilename);
  
  bool isValid();

  void setFile(DataFormat *format, const char *gridfilename);
  
  virtual bool hasVar(const string &name);
  
  virtual vector<int> getSize(const string &name);

  virtual bool setGlobalOrigin(int x = 0, int y = 0, int z = 0);

  // Fetch the data. Returns true on success
  virtual bool fetch(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  virtual bool fetch(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0);
  virtual bool fetch(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  virtual bool fetch(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0);

  void open(const char* file);
  void close();
 private:
  DataFormat *file;
  string filename;
  
  bool isOpen;
};

class GridFromOptions : GridDataSource {
public:
  GridFromOptions(Options *opt = NULL); 
  
  bool hasVar(const string &name);
  
  vector<int> getSize(const string &name); ///< Get size of the variable

  /// Set the (x,y,z) origin for all subsequent calls
  bool setGlobalOrigin(int x = 0, int y = 0, int z = 0) {return true;}

  bool setMesh(Mesh *m) {fieldmesh = m; return true;}
  
  /// Get data from the source
  bool fetch(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool fetch(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0);
  bool fetch(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool fetch(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0);
private:
  Options *options;
  Mesh *fieldmesh;
};

class GridDataGroup : GridDataSource {
public:
  GridDataGroup() {}
  GridDataGroup(GridDataSource *a, GridDataSource *b = NULL, GridDataSource *c = NULL, GridDataSource *d = NULL) {
    add(a); add(b); add(c); add(d);
  }
  
  /// Add a data source
  void add(GridDataSource &source) {add(&source);}
  void add(GridDataSource *source);
  
private:
  std::list<GridDataSource*> source_list; ///< List of sources
  
  GridDataSource *findSource(const char *name);
  GridDataSource *findSource(const string &name) {return findSource(name.c_str());}
};

#endif // __GRIDDATA_H__
