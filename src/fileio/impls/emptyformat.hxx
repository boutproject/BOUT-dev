/**************************************************************************
 * Empty format class for throwing errors
 * 
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
 **************************************************************************/

class EmptyFormat;

#ifndef __EMPTYSOLVER_H__
#define __EMPTYSOLVER_H__

#include <dataformat.hxx>
#include <boutexception.hxx>

class EmptyFormat {
  EmptyFormat() {throw BoutException("File format not enabled!");}
  
  bool openr(const string &name) {return false; }
  bool openw(const string &name, bool append) {return false; }
  
  bool is_valid() {return false;}
  
  void close() {}
  
  const vector<int> getSize(const char *var) {vector<int> tmp; return tmp;}
  const vector<int> getSize(const string &var) {vector<int> tmp; return tmp;}
  
  bool setOrigin(int x = 0, int y = 0, int z = 0) {return false;}
  bool setRecord(int t) {return false;}
  
  bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0)        {return false;}
  bool read(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0)      {return false;}
  bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0)   {return false;}
  bool read(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0) {return false;}
  
  bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) {return false;}
  bool write(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0) {return false;}
  bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) {return false;}
  bool write(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0) {return false;}
  
  bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) {return false;}
  bool read_rec(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0) {return false;}
  bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) {return false;}
  bool read_rec(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0) {return false;}

  bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) {return false;}
  bool write_rec(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0) {return false;}
  bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) {return false;}
  bool write_rec(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0) {return false;}
};

#endif // __EMPTYSOLVER_H__
