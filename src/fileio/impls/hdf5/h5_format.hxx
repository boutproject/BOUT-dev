/*!
 * \file h5_format.hxx
 *
 * \brief HDF5 data format interface
 *
 * \author John Omotani
 * \date   2015
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

#ifndef HDF5

#include "../emptyformat.hxx"
using H5Format = EmptyFormat;

#else

class H5Format;

#ifndef __H5FORMAT_H__
#define __H5FORMAT_H__

#include "dataformat.hxx"

#include <hdf5.h>

#include <map>
#include <string>

class H5Format : public DataFormat {
 public:
  H5Format(bool parallel_in = false, Mesh* mesh_in = nullptr);
  H5Format(const char *name, bool parallel_in = false, Mesh* mesh_in = nullptr);
  H5Format(const std::string &name, bool parallel_in = false, Mesh* mesh_in = nullptr)
    : H5Format(name.c_str(), parallel_in, mesh_in) {}
  ~H5Format();

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
  bool setLocalOrigin(int x = 0, int y = 0, int z = 0, int offset_x = 0, int offset_y = 0, int offset_z = 0) override;
  bool setRecord(int t) override; // negative -> latest

  // Add a variable to the file
  bool addVarInt(const std::string &name, bool repeat) override;
  bool addVarIntVec(const std::string &name, bool repeat, size_t size) override;
  bool addVarBoutReal(const std::string &name, bool repeat) override;
  bool addVarField2D(const std::string &name, bool repeat) override;
  bool addVarField3D(const std::string &name, bool repeat) override;
  bool addVarFieldPerp(const std::string &name, bool repeat) override;
  
  // Read / Write simple variables up to 3D

  bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read_perp(BoutReal *var, const std::string &name, int lx = 1, int lz = 0) override;

  bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write_perp(BoutReal *var, const std::string &name, int lx = 0, int lz = 0) override;

  // Read / Write record-based variables

  bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read_rec(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read_rec(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) override;
  bool read_rec_perp(BoutReal *var, const std::string &name, int lx = 1, int lz = 0) override;

  bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write_rec(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write_rec(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) override;
  bool write_rec_perp(BoutReal *var, const std::string &name, int lx = 0, int lz = 0) override;
  
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
  
  hid_t dataFile;
  hid_t dataFile_plist;
  hid_t dataSet_plist;

  bool lowPrecision; ///< When writing, down-convert to floats
  bool parallel;

  int x0, y0, z0, t0; ///< Data origins for file access
  int x0_local, y0_local, z0_local; ///< Data origins for memory access
  
  hsize_t chunk_length;

  bool addVar(const std::string &name, bool repeat, hid_t write_hdf5_type, std::string datatype,
              int lx = 0, int ly = 0, int lz = 0);
  bool read(void *var, hid_t hdf5_type, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool write(void *var, hid_t mem_hdf5_type, const char *name, int lx = 0, int ly = 0, int lz = 0);
  bool read_rec(void *var, hid_t hdf5_type, const char *name, int lx = 1, int ly = 0, int lz = 0);
  bool write_rec(void *var, hid_t mem_hdf5_type, const char *name, int lx = 0, int ly = 0, int lz = 0);

  // Attributes

  void setAttribute(const hid_t &dataSet, const std::string &attrname,
                    const std::string &text);
  void setAttribute(const hid_t &dataSet, const std::string &attrname,
                    int value);
  void setAttribute(const hid_t &dataSet, const std::string &attrname,
                    BoutReal value);
  bool getAttribute(const hid_t &dataSet, const std::string &attrname, std::string &text);
  bool getAttribute(const hid_t &dataSet, const std::string &attrname, int &value);
  bool getAttribute(const hid_t &dataSet, const std::string &attrname, BoutReal &value);
};

#endif // __H5FORMAT_H__

#endif // HDF5
