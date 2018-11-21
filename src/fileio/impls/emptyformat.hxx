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
#include <unused.hxx>

class EmptyFormat {
  EmptyFormat() {throw BoutException("File format not enabled!");}
  
  bool openr(const std::string &UNUSED(name)) {return false; }
  bool openw(const std::string &UNUSED(name), bool UNUSED(append)) {return false; }
  
  bool is_valid() {return false;}
  
  void close() {}
  
  const std::vector<int> getSize(const char *UNUSED(var)) {std::vector<int> tmp; return tmp;}
  const std::vector<int> getSize(const std::string &UNUSED(var)) {std::vector<int> tmp; return tmp;}
  
  bool setOrigin(int UNUSED(x) = 0, int UNUSED(y) = 0, int UNUSED(z) = 0) {return false;}
  bool setRecord(int UNUSED(t)) {return false;}
  
  bool read(int *UNUSED(var), const char *UNUSED(name), int UNUSED(lx) = 1,
            int UNUSED(ly) = 0, int UNUSED(lz) = 0)        {return false;}
  bool read(int *UNUSED(var), const std::string &UNUSED(name), int UNUSED(lx) = 1,
            int UNUSED(ly) = 0, int UNUSED(lz) = 0)      {return false;}
  bool read(BoutReal *UNUSED(var), const char *UNUSED(name), int UNUSED(lx) = 1,
            int UNUSED(ly) = 0, int UNUSED(lz) = 0)   {return false;}
  bool read(BoutReal *UNUSED(var), const std::string &UNUSED(name), int UNUSED(lx) = 1,
            int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  
  bool write(int *UNUSED(var), const char *UNUSED(name), int UNUSED(lx) = 0,
             int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  bool write(int *UNUSED(var), const std::string &UNUSED(name), int UNUSED(lx) = 0,
             int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  bool write(BoutReal *UNUSED(var), const char *UNUSED(name), int UNUSED(lx) = 0,
             int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  bool write(BoutReal *UNUSED(var), const std::string &UNUSED(name), int UNUSED(lx) = 0,
             int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  
  bool read_rec(int *UNUSED(var), const char *UNUSED(name), int UNUSED(lx) = 1,
                int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  bool read_rec(int *UNUSED(var), const std::string &UNUSED(name), int UNUSED(lx) = 1,
                int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  bool read_rec(BoutReal *UNUSED(var), const char *UNUSED(name), int UNUSED(lx) = 1,
                int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  bool read_rec(BoutReal *UNUSED(var), const std::string &UNUSED(name), int UNUSED(lx) = 1,
                int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}

  bool write_rec(int *UNUSED(var), const char *UNUSED(name), int UNUSED(lx) = 0,
                 int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  bool write_rec(int *UNUSED(var), const std::string &UNUSED(name), int UNUSED(lx) = 0,
                 int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  bool write_rec(BoutReal *UNUSED(var), const char *UNUSED(name), int UNUSED(lx) = 0,
                 int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
  bool write_rec(BoutReal *UNUSED(var), const std::string &UNUSED(name), int UNUSED(lx) = 0,
                 int UNUSED(ly) = 0, int UNUSED(lz) = 0) {return false;}
};

#endif // __EMPTYSOLVER_H__
