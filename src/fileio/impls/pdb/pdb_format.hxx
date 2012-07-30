/*!
 * \file pdb_format.h
 *
 * \brief PDB data format interface
 *
 * \author B.Dudson
 * \date   April 2009
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
 */

#ifndef PDBF

#include "../emptyformat.hxx"
typedef EmptyFormat PdbFormat;

#else

class PdbFormat;

#ifndef __PDBFORMAT_H__
#define __PDBFORMAT_H__

#include "dataformat.hxx"

#include "pdb.h"

class PdbFormat : public DataFormat {
 public:
  PdbFormat();
  PdbFormat(const char *name);
  PdbFormat(const string &name);
  ~PdbFormat();
  
  bool openr(const string &name);
  bool openr(const char *name);
  bool openw(const string &name, bool append=false);
  bool openw(const char *name, bool append=false);
  
  virtual bool is_valid();
  
  void close();
  
  void flush();
  
  const char* filename();

  const vector<int> getSize(const char *var);
  const vector<int> getSize(const string &var);

  // Set the origin for all subsequent calls
  bool setGlobalOrigin(int x = 0, int y = 0, int z = 0);
  bool setLocalOrigin(int x = 0, int y = 0, int z = 0) { return setGlobalOrigin(x,y,z); }
  
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
  PDBfile *fp;
  char *fname;
  bool appending;
  bool lowPrecision; ///< When writing, down-convert to floats
  
  int x0, y0, z0, t0; // Data origins

  int nrecs(const char *name); // Returns the number of records
};

#endif // __PDBFORMAT_H__

#endif // PDBF

