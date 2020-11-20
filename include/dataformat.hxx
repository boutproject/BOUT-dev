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
#include <vector>

class Mesh;
class Field;
class FieldPerp;

// Can't copy, to control access to file
class DataFormat {
 public:
  DataFormat(Mesh* mesh_in = nullptr);
  virtual ~DataFormat() = default;
  // File opening routines
  virtual bool openr(const char *name) = 0;
  virtual bool openr(const std::string &name) {
    return openr(name.c_str());
  }
  virtual bool openr(const std::string &base, int mype);
  virtual bool openw(const char *name, bool append=false) = 0;
  virtual bool openw(const std::string &name, bool append=false) {
    return openw(name.c_str(), append);
  }
  virtual bool openw(const std::string &base, int mype, bool append=false);
  
  virtual bool is_valid() = 0;
  
  virtual void close() = 0;

  virtual void flush() = 0;

  virtual const std::vector<int> getSize(const char *var) = 0;
  virtual const std::vector<int> getSize(const std::string &var) = 0;

  // Set the origin for all subsequent calls
  virtual bool setGlobalOrigin(int x = 0, int y = 0, int z = 0) = 0; 
  virtual bool setLocalOrigin(int x = 0, int y = 0, int z = 0, int offset_x = 0, int offset_y = 0, int offset_z = 0);
  virtual bool setRecord(int t) = 0; // negative -> latest

  // Add a variable to the file
  virtual bool addVarInt(const std::string &name, bool repeat) = 0;
  virtual bool addVarIntVec(const std::string &name, bool repeat, size_t size) = 0;
  virtual bool addVarBoutReal(const std::string &name, bool repeat) = 0;
  virtual bool addVarField2D(const std::string &name, bool repeat) = 0;
  virtual bool addVarField3D(const std::string &name, bool repeat) = 0;
  virtual bool addVarFieldPerp(const std::string &name, bool repeat) = 0;
  
  // Read / Write simple variables up to 3D

  virtual bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read_perp(BoutReal *var, const std::string &name, int lx = 1, int lz = 0) = 0;

  virtual bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write_perp(BoutReal *var, const std::string &name, int lx = 0, int lz = 0) = 0;

  // Read / Write record-based variables

  virtual bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read_rec(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read_rec(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read_rec_perp(BoutReal *var, const std::string &name, int lx = 1, int lz = 0) = 0;

  virtual bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write_rec(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write_rec(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write_rec_perp(BoutReal *var, const std::string &name, int lx = 0, int lz = 0) = 0;

  // Optional functions
  
  virtual void setLowPrecision() { }  // By default doesn't do anything

  // Attributes

  /// Sets a string attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist. If
  ///                        varname is the empty string "" then the attribute
  ///                        will be added to the file instead of to a
  ///                        variable.
  /// @param[in] attrname    Attribute name
  /// @param[in] text        A string attribute to attach to the variable
  virtual void setAttribute(const std::string &varname, const std::string &attrname,
                            const std::string &text) = 0;

  /// Sets an integer attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist. If
  ///                        varname is the empty string "" then the attribute
  ///                        will be added to the file instead of to a
  ///                        variable.
  /// @param[in] attrname    Attribute name
  /// @param[in] value       An int attribute to attach to the variable
  virtual void setAttribute(const std::string &varname, const std::string &attrname,
                            int value) = 0;

  /// Sets a BoutReal attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist. If
  ///                        varname is the empty string "" then the attribute
  ///                        will be added to the file instead of to a
  ///                        variable.
  /// @param[in] attrname    Attribute name
  /// @param[in] value       A BoutReal attribute to attach to the variable
  virtual void setAttribute(const std::string &varname, const std::string &attrname,
                            BoutReal value) = 0;

  /// Gets a string attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist. If
  ///                        varname is the empty string "" then get the
  ///                        attribute from the top-level of the file instead
  ///                        of from a variable.
  /// @param[in] attrname    Attribute name
  ///
  /// Returns
  /// -------
  /// text                   A string attribute of the variable
  virtual bool getAttribute(const std::string &varname, const std::string &attrname, std::string &text) = 0;

  /// Gets an integer attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist. If
  ///                        varname is the empty string "" then get the
  ///                        attribute from the top-level of the file instead
  ///                        of from a variable.
  /// @param[in] attrname    Attribute name
  ///
  /// Returns
  /// -------
  /// value                  An int attribute of the variable
  virtual bool getAttribute(const std::string &varname, const std::string &attrname, int &value) = 0;

  /// Gets a BoutReal attribute
  ///
  /// Inputs
  /// ------
  ///
  /// @param[in] varname     Variable name. The variable must already exist. If
  ///                        varname is the empty string "" then get the
  ///                        attribute from the top-level of the file instead
  ///                        of from a variable.
  /// @param[in] attrname    Attribute name
  ///
  /// Returns
  /// -------
  /// value                  A BoutReal attribute of the variable
  virtual bool getAttribute(const std::string &varname, const std::string &attrname, BoutReal &value) = 0;

  /// Write out the meta-data of a field as attributes of the variable
  void writeFieldAttributes(const std::string& name, const Field& f);
  /// Overload for FieldPerp so we can also write 'yindex'
  void writeFieldAttributes(const std::string& name, const FieldPerp& f);

  /// Read the attributes of a field
  void readFieldAttributes(const std::string& name, Field& f);
  /// Overload for FieldPerp so we can also read 'yindex'
  void readFieldAttributes(const std::string& name, FieldPerp& f);

 protected:
  Mesh* mesh;
};

// For backwards compatability. In formatfactory.cxx
std::unique_ptr<DataFormat> data_format(const char *filename = nullptr);

#endif // __DATAFORMAT_H__
