/*!
 * \file dataformat.hxx
 *
 * \brief Generic interface for file formats e.g. netCDF, HDF5
 *
 * \author B.Dudson
 * \date   April 2009
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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

class DataFormat;

#ifndef __DATAFORMAT_H__
#define __DATAFORMAT_H__

#include "bout_types.hxx"
#include "unused.hxx"

#include <string>
#include <memory>
using std::string;

#include <vector>
using std::vector;

// Can't copy, to control access to file
class DataFormat {
 public:
  virtual ~DataFormat() { }
  // File opening routines
  virtual bool openr(const char *name) = 0;
  virtual bool openr(const string &name) {
    return openr(name.c_str());
  }
  virtual bool openr(const string &base, int mype);
  virtual bool openw(const char *name, bool append=false) = 0;
  virtual bool openw(const string &name, bool append=false) {
    return openw(name.c_str(), append);
  }
  virtual bool openw(const string &base, int mype, bool append=false);
  
  virtual bool is_valid() = 0;
  
  virtual void close() = 0;

  virtual void flush() = 0;

  virtual const vector<int> getSize(const char *var) = 0;
  virtual const vector<int> getSize(const string &var) = 0;

  // Set the origin for all subsequent calls
  virtual bool setGlobalOrigin(int x = 0, int y = 0, int z = 0) = 0; 
  virtual bool setLocalOrigin(int x = 0, int y = 0, int z = 0, int offset_x = 0, int offset_y = 0, int offset_z = 0);
  virtual bool setRecord(int t) = 0; // negative -> latest

  // Add a variable to the file
  virtual bool addVarInt(const string &name, bool repeat) = 0;
  virtual bool addVarBoutReal(const string &name, bool repeat) = 0;
  virtual bool addVarField2D(const string &name, bool repeat) = 0;
  virtual bool addVarField3D(const string &name, bool repeat) = 0;
  
  // Read / Write simple variables up to 3D

  virtual bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;

  virtual bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;

  // Read / Write record-based variables

  virtual bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read_rec(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read_rec(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;

  virtual bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write_rec(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write_rec(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;

  // Optional functions
  
  virtual void setLowPrecision() { }  // By default doesn't do anything

  // Attributes

  /// Sets a string attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist
  /// @param[in] attrname    Attribute name
  /// @param[in] text        A string attribute to attach to the variable
  virtual void setAttribute(const string &varname, const string &attrname,
                            const string &text) = 0;

  /// Sets an integer attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist
  /// @param[in] attrname    Attribute name
  /// @param[in] value       A string attribute to attach to the variable
  virtual void setAttribute(const string &varname, const string &attrname,
                            int value) = 0;
  /// Gets a string attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist
  /// @param[in] attrname    Attribute name
  ///
  /// Returns
  /// -------
  /// text                   A string attribute of the variable
  virtual bool getAttribute(const string &varname, const string &attrname, std::string &text) = 0;

  /// Sets an integer attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist
  /// @param[in] attrname    Attribute name
  ///
  /// Returns
  /// -------
  /// value                  An int attribute of the variable
  virtual bool getAttribute(const string &varname, const string &attrname, int &value) = 0;
};

// For backwards compatability. In formatfactory.cxx
std::unique_ptr<DataFormat> data_format(const char *filename = nullptr);

#endif // __DATAFORMAT_H__
