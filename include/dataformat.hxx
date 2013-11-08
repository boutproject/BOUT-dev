/*!
 * \file dataformat.hxx
 *
 * \brief Generic interface for file formats e.g. PDB, netCDF
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
#include <string>
using std::string;

#include <vector>
using std::vector;

// Can't copy, to control access to file
class DataFormat {
 public:
  virtual ~DataFormat() { }
  // File opening routines
  virtual bool openr(const string &name) = 0;
  virtual bool openr(const string &base, int mype);
  virtual bool openw(const string &name, bool append=false) = 0;
  virtual bool openw(const string &base, int mype, bool append=false);

  virtual bool is_valid() = 0;
  
  virtual void close() = 0;

  virtual void flush() = 0;

  virtual const vector<int> getSize(const string &var) = 0;

  // Set the origin for all subsequent calls
  virtual bool setGlobalOrigin(int x = 0, int y = 0, int z = 0) = 0; 
  virtual bool setLocalOrigin(int x = 0, int y = 0, int z = 0);
  virtual bool setRecord(int t) = 0; // negative -> latest
  
  // Read / Write simple variables up to 3D

  virtual bool read(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;

  virtual bool write(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;

  // Read / Write record-based variables

  virtual bool read_rec(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
  virtual bool read_rec(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;

  virtual bool write_rec(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
  virtual bool write_rec(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;

  // Optional functions
  
  virtual void setLowPrecision() { }  // By default doesn't do anything
};

// For backwards compatability. In formatfactory.cxx
DataFormat* data_format(const char *filename = NULL);

#endif // __DATAFORMAT_H__
