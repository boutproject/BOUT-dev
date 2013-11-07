/*!
 * \file nc_capi.hxx
 *
 * \brief netCDF data format interface, using C API
 *
 * \author B.Dudson
 * \date   October 2013
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

#ifndef NCDFCAPI

#include "../emptyformat.hxx"
typedef EmptyFormat NcdfCapi;

#else

class NcdfCapi;

#ifndef __NCFORMATCAPI_H__
#define __NCFORMATCAPI_H__

#include "dataformat.hxx"

#include <netcdf.h>

#include <map>
#include <string>

class NcdfCapi : public DataFormat {
 public:
  NcdfCapi();
  NcdfCapi(const std::string &name);
  ~NcdfCapi();
  
  bool openr(const std::string &name);
  bool openw(const std::string &name, bool append=false);
  
  bool is_valid();
  
  void close();
  
  void flush();

  const vector<int> getSize(const string &var);
  
  // Set the origin for all subsequent calls
  bool setGlobalOrigin(int x = 0, int y = 0, int z = 0);
  bool setLocalOrigin(int x = 0, int y = 0, int z = 0) { return setGlobalOrigin(x,y,z); }
  bool setRecord(int t); // negative -> latest
  
  // Read / Write simple variables up to 3D

  bool read(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0);
  bool read(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0);

  bool write(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0);
  bool write(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0);

  // Read / Write record-based variables

  bool read_rec(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0);
  bool read_rec(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0);

  bool write_rec(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0);
  bool write_rec(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0);
  
  void setLowPrecision() { lowPrecision = true; }

 private:
  /// netCDF file ID
  int datafile;
  
  /// Dimensions
  int xDim, yDim, zDim, tDim;
  int *dimList;      ///< List of dimensions (x,y,z)
  int recDimList[4]; ///< List of dimensions (t,x,y,z)

  bool appending;
  bool lowPrecision; ///< When writing, down-convert to floats

  int x0, y0, z0, t0; ///< Data origins

  std::map<std::string, int> rec_nr; // Record number for each variable (bit nasty)
  int default_rec;  // Starting record. Useful when appending to existing file
  
};

#endif // __NCFORMATCAPI_H__

#endif // NCDFCAPI
