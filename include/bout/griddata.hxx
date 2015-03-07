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

#include <field2d.hxx>
#include <field3d.hxx>

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

  virtual bool hasVar(const string &name) = 0; ///< Test if source can supply a variable

  virtual bool get(Mesh *m, int &ival,      const string &name) = 0; ///< Get an integer
  virtual bool get(Mesh *m, BoutReal &rval, const string &name) = 0; ///< Get a BoutReal number
  virtual bool get(Mesh *m, Field2D &var,   const string &name, BoutReal def=0.0) = 0;
  virtual bool get(Mesh *m, Field3D &var,   const string &name, BoutReal def=0.0) = 0;

  enum Direction {X=1, Y=2, Z=3};
  virtual bool get(Mesh *m, vector<int> &var,      const string &name, int len, int offset=0, Direction dir = GridDataSource::X) = 0;
  virtual bool get(Mesh *m, vector<BoutReal> &var, const string &name, int len, int offset=0, Direction dir = GridDataSource::X) = 0;
};

/// Interface to grid data in a file
/*!
 * This is a thin wrapper around a DataFormat object. Only needs to implement 
 * reading routines.
 */
class GridFile : public GridDataSource {
 public:
  GridFile(DataFormat *format, const string gridfilename);
  ~GridFile();
  
  bool hasVar(const string &name);

  bool get(Mesh *m, int &ival,      const string &name); ///< Get an integer
  bool get(Mesh *m, BoutReal &rval, const string &name); ///< Get a BoutReal number
  bool get(Mesh *m, Field2D &var,   const string &name, BoutReal def=0.0);
  bool get(Mesh *m, Field3D &var,   const string &name, BoutReal def=0.0);

  bool get(Mesh *m, vector<int> &var,      const string &name, int len, int offset=0, GridDataSource::Direction dir = GridDataSource::X);
  bool get(Mesh *m, vector<BoutReal> &var, const string &name, int len, int offset=0, GridDataSource::Direction dir = GridDataSource::X);
  
 private:
  GridFile();
  
  DataFormat *file;
  string filename;

  bool readgrid_3dvar_fft(Mesh *m, const string &name, 
			  int yread, int ydest, int ysize, 
			  int xge, int xlt, BoutReal ***var);
  
  bool readgrid_3dvar_real(Mesh *m, const string &name, 
			   int yread, int ydest, int ysize, 
			   int xge, int xlt, BoutReal ***var);
};

class GridFromOptions : GridDataSource {
public:
  GridFromOptions(Options *opt = NULL) : options(opt) {}
  
  bool hasVar(const string &name);
  
  bool get(Mesh *m, int &ival,      const string &name); ///< Get an integer
  bool get(Mesh *m, BoutReal &rval, const string &name); ///< Get a BoutReal number
  bool get(Mesh *m, Field2D &var,   const string &name, BoutReal def=0.0);
  bool get(Mesh *m, Field3D &var,   const string &name, BoutReal def=0.0);
  
  bool get(Mesh *m, vector<int> &var,      const string &name, int len, int offset=0, GridDataSource::Direction dir = GridDataSource::X);
  bool get(Mesh *m, vector<BoutReal> &var, const string &name, int len, int offset=0, GridDataSource::Direction dir = GridDataSource::X);
  
private:
  Options *options;
};

/*
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
*/
#endif // __GRIDDATA_H__
