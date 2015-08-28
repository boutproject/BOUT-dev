/*!
 * \file ncxx4.hxx
 *
 * \brief netCDF-4 data format interface
 *
 * \author B.Dudson
 * \date   September 2012
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

#ifndef NCDF4

#include "../emptyformat.hxx"
typedef EmptyFormat Ncxx4;

#else

class Ncxx4;

#ifndef __NCFORMAT4_H__
#define __NCFORMAT4_H__

#include "dataformat.hxx"

#include <netcdf>

#include <map>
#include <string>
#include <vector>

class Ncxx4 : public DataFormat {
public:
  Ncxx4();
  Ncxx4(const char *name);
  Ncxx4(const std::string &name);
  ~Ncxx4();

  bool openr(const std::string &name);
  bool openr(const char *name);
  bool openw(const std::string &name, bool append=false);
  bool openw(const char *name, bool append=false);

  bool is_valid();

  void close();

  void flush();

  const char* filename() { return fname; };

  const vector<int> getSize(const char *var);
  const vector<int> getSize(const string &var);

  // Set the origin for all subsequent calls
  bool setGlobalOrigin(int x = 0, int y = 0, int z = 0);
  bool setLocalOrigin(int x = 0, int y = 0, int z = 0) { return setGlobalOrigin(x,y,z); }
  bool setRecord(int t); // negative -> latest

  // Read / Write simple variables up to 3D

  bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool read(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0);
  bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool read(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0);

  bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
  bool write(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0);
  bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
  bool write(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0);

  // Read / Write record-based variables

  bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool read_rec(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0);
  bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool read_rec(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0);

  bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
  bool write_rec(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0);
  bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
  bool write_rec(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0);

  void setLowPrecision() { lowPrecision = true; }

private:

  char *fname; ///< Current file name

  /// Pointer to netCDF file
  netCDF::NcFile *dataFile;

  /// Dimensions
  netCDF::NcDim xDim, yDim, zDim, tDim;
  const netCDF::NcDim **dimList; ///< List of dimensions (x,y,z)
  const netCDF::NcDim **recDimList; ///< List of dimensions (t,x,y,z)

  bool appending;
  bool lowPrecision; ///< When writing, down-convert to floats

  int x0, y0, z0, t0; ///< Data origins

  std::map<std::string, int> rec_nr; // Record number for each variable (bit nasty)
  int default_rec;  // Starting record. Useful when appending to existing file

  std::vector<netCDF::NcDim> getDimVec(int nd);
  std::vector<netCDF::NcDim> getRecDimVec(int nd);
};

#endif // __NCFORMAT4_H__

#endif // NCDF4
