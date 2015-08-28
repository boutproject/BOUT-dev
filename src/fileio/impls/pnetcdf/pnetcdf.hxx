/*!
 * \file pnetcdf.hxx
 *
 * \brief Parallel NetCDF data format interface
 *
 * \author B.Dudson
 * \date   May 2012
 *
 * Records: In netCDF, the time dimension for each dimension must be
 * the same. Hence when a record is appended to a variable, the size
 * of all variables is increased. To work out which record to write to,
 * a map of variable names to record number is kept.
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

#ifndef PNCDF

#include "../emptyformat.hxx"
typedef EmptyFormat PncFormat;

#else

class PncFormat;

#ifndef __PNCFORMAT_H__
#define __PNCFORMAT_H__

#include "dataformat.hxx"

#include <map>
#include <string>

using std::string;
using std::map;

class PncFormat : public DataFormat {
public:
  PncFormat();
  PncFormat(const char *name);
  PncFormat(const string &name);
  ~PncFormat();

  bool openr(const string &name) { return openr(name.c_str()); }
  bool openr(const char *name);
  bool openr(const string &name, int mype) {return openr(name);}

  bool openw(const string &name, bool append=false) { return openw(name.c_str(), append); }
  bool openw(const char *name, bool append=false);
  bool openw(const string &name, int mype, bool append=false) {return openw(name, append);}

  bool is_valid() { return fname != NULL;}

  void close();

  void flush();

  const char* filename() { return fname; };

  const vector<int> getSize(const char *var);
  const vector<int> getSize(const string &var) { return getSize(var.c_str()); }

  // Set the origin for all subsequent calls
  bool setGlobalOrigin(int x = 0, int y = 0, int z = 0);
  bool setRecord(int t); // negative -> latest

  // Read / Write simple variables up to 3D

  bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool read(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0);
  bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool read(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0);

  bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
  bool write(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0);
  bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
  bool write(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0);

  // Read / Write record-based variables

  bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool read_rec(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0);
  bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool read_rec(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0);

  bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
  bool write_rec(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0);
  bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
  bool write_rec(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0);

  void setLowPrecision() { lowPrecision = true; }

private:

  char *fname; ///< Current file name

  /// ID of netCDF file
  int ncfile;
  bool valid; // True if the file is valid

  /// Dimensions
  int xDim, yDim, zDim, tDim;

  int *dimList; ///< List of dimensions (x,y,z)
  int recDimList[4]; ///< List of dimensions (t,x,y,z)

  bool appending;
  bool lowPrecision; ///< When writing, down-convert to floats

  int x0, y0, z0, t0; ///< Data origins (global offsets)

  map<string, int> rec_nr; // Record number for each variable (bit nasty)
  int default_rec;  // Starting record. Useful when appending to existing file
};

#endif // __PNCFORMAT_H__

#endif // PNCDF
