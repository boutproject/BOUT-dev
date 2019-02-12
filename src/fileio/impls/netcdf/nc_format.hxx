/*!
 * \file nc_format.hxx
 *
 * \brief netCDF data format interface
 *
 * \author B.Dudson
 * \date   April 2009
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

#ifndef NCDF

#include "../emptyformat.hxx"
typedef EmptyFormat NcFormat;

#else

class NcFormat;

#ifndef __NCFORMAT_H__
#define __NCFORMAT_H__

#include "dataformat.hxx"
#include "unused.hxx"

#include <netcdfcpp.h>

#include <map>
#include <string>

class NcFormat : public DataFormat {
 public:
  NcFormat(Mesh* mesh_in = nullptr);
  NcFormat(const char *name, Mesh* mesh_in = nullptr);
  NcFormat(const std::string &name, Mesh* mesh_in = nullptr)
    : NcFormat(name.c_str(), mesh_in) {}
  ~NcFormat();

  using DataFormat::openr;
  bool openr(const char *name) override;
  using DataFormat::openw;
  bool openw(const char *name, bool append=false) override;
  
  bool is_valid() override;
  
  void close() override;
  
  void flush() override;

  const char* filename() { return fname; };

  const std::vector<int> getSize(const char *var) override;
  const std::vector<int> getSize(const std::string &var) override;
  
  // Set the origin for all subsequent calls
  bool setGlobalOrigin(int x = 0, int y = 0, int z = 0) override;
  bool setLocalOrigin(int x = 0, int y = 0, int z = 0,
                      int UNUSED(offset_x) = 0, int UNUSED(offset_y) = 0, int UNUSED(offset_z) = 0) override {
    return setGlobalOrigin(x,y,z);
  }
  bool setRecord(int t) override; // negative -> latest
  
  // Add a variable to the file
  bool addVarInt(const std::string &name, bool repeat) override;
  bool addVarBoutReal(const std::string &name, bool repeat) override;
  bool addVarField2D(const std::string &name, bool repeat) override;
  bool addVarField3D(const std::string &name, bool repeat) override;

  // Read / Write simple variables up to 3D

  bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) override;

  bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) override;

  // Read / Write record-based variables

  bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read_rec(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read_rec(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) override;

  bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write_rec(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write_rec(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) override;
  
  void setLowPrecision() override { lowPrecision = true; }

  // Attributes

  void setAttribute(const std::string &varname, const std::string &attrname,
                    const std::string &text) override;
  void setAttribute(const std::string &varname, const std::string &attrname,
                    int value) override;
  void setAttribute(const std::string &varname, const std::string &attrname,
                    BoutReal value) override;
  bool getAttribute(const std::string &varname, const std::string &attrname, std::string &text) override;
  bool getAttribute(const std::string &varname, const std::string &attrname, int &value) override;
  bool getAttribute(const std::string &varname, const std::string &attrname, BoutReal &value) override;
  
 private:

  char *fname; ///< Current file name

  /// Pointer to netCDF file
  NcFile *dataFile;
  
  /// Dimensions
  NcDim *xDim, *yDim, *zDim, *tDim;
  const NcDim **dimList; ///< List of dimensions (x,y,z)
  const NcDim **recDimList; ///< List of dimensions (t,x,y,z)

  bool lowPrecision; ///< When writing, down-convert to floats

  int x0, y0, z0, t0; ///< Data origins

  std::map<std::string, int> rec_nr; // Record number for each variable (bit nasty)
  int default_rec;  // Starting record. Useful when appending to existing file

  void checkName(const char* name); ///< Check if a name contains invalid characters
  
};

#endif // __NCFORMAT_H__

#endif // NCDF
