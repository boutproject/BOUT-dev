/*!
 * \file griddata.cxx
 * \brief Functions to read variables from a grid file
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

#include "mpi.h"

#include <utils.hxx>

#include <dataformat.hxx> // Abstract data class

#include <string.h>
#include <stdlib.h>
#include <cmath>

#include <fft.hxx> // For reading 3D variables
#include <dcomplex.hxx>

#include <bout_types.hxx>

#include <output.hxx>

#include <bout/griddata.hxx>

#include <msg_stack.hxx>

/*******************************************************************************
 * GridFile class
 * 
 * This is a fairly thin wrapper around the Datafile class. Only reading
 * routines are required though. Takes care of opening and closing file
 * then just passes read requests through.
 *******************************************************************************/

/// Default constructor
GridFile::GridFile() {
  file = NULL; ///< Signal that the file is invalid
  isOpen = false;
}

GridFile::~GridFile() {
  if(file != NULL)
    delete file;
}

/// Constructor passing a file format and name
GridFile::GridFile(DataFormat *format, const char *gridfilename) {
  file = NULL; ///< Signal that the file is invalid
  isOpen = false;
  setFile(format, gridfilename);
}

bool GridFile::isValid() {
  if(!file)
    return false;
  
  if(!isOpen) {
    // File not open, so need to open
    if(!file->openr(filename))
      return false;
  }
  
  bool valid = file->is_valid();
  
  if(!isOpen)
    file->close();

  return valid;
}

void GridFile::setFile(DataFormat *format, const char *gridfilename) {
  if(file != NULL)
    delete file;
  file = format;
  filename = string(gridfilename);
  isOpen = false;
}

bool GridFile::hasVar(const string &name) {
  if(file == NULL)
    return false;

  /// Try to get the size of the variable
  vector<int> s = getSize(name);
  
  /// Test if the variable has zero size
  return s.size() != 0;
}

/// Return the size of variable
vector<int> GridFile::getSize(const string &name) {
  vector<int> s;
  
  if(file == NULL)
    return s;
  
  if(!isOpen) {
    // File not open, so need to open
    if(!file->openr(filename))
      return s;
  }

  s = file->getSize(name);
  
  if(!isOpen)
    file->close(); // If was closed before then close again
  
  return s;
}

bool GridFile::setGlobalOrigin(int x, int y, int z) {
  if(file == NULL)
    return false;
  
  return file->setGlobalOrigin(x,y,z);
}

bool GridFile::fetch(int *var, const char *name, int lx, int ly, int lz) {
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

bool GridFile::fetch(int *var, const string &name, int lx, int ly, int lz) {
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

bool GridFile::fetch(BoutReal *var, const char *name, int lx, int ly, int lz) {
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

bool GridFile::fetch(BoutReal *var, const string &name, int lx, int ly, int lz) {
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

void GridFile::open(const char *name) {
  file->openr(filename);
  if(!file->is_valid()) {
    output << "\tERROR: Could not open file " << filename << endl;
  }
  isOpen = true;
}

void GridFile::close() {
  file->close();
  isOpen = false;
}

/*******************************************************************************
 * GridFromOptions class
 * 
 * 
 *******************************************************************************/

GridFromOptions::GridFromOptions(Options *opt) : options(opt), fieldmesh(NULL) {
  if(options == NULL)
    options = Options::getRoot()->getSection("mesh");
}

bool GridFromOptions::hasVar(const string &name) {
  return options->isSet(name);
}

vector<int> GridFromOptions::getSize(const string &name) {
  vector<int> v(1);
  v[0] = 1;
  return v;
}

bool GridFromOptions::fetch(int *var, const char *name, int lx, int ly, int lz) {
  // Get the value from options. Could improve to use FieldGenerator?
  int value;
  options->get(name, value, 0);
  if(lx < 1)
    lx = 1;
  if(ly < 1)
    ly = 1;
  if(lz < 1)
    lz = 1;
  for(int i=0;i<lx*ly*lz;i++)
    var[i] = value;
  
  return true;
}

bool GridFromOptions::fetch(int *var, const string &name, int lx, int ly, int lz) {
  return fetch(var, name.c_str(), lx, ly, lz);
}

bool GridFromOptions::fetch(BoutReal *var, const char *name, int lx, int ly, int lz) {
  BoutReal value;
  options->get(name, value, 0);
  
  if(lx < 1)
    lx = 1;
  if(ly < 1)
    ly = 1;
  if(lz < 1)
    lz = 1;
  
  for(int i=0;i<lx*ly*lz;i++)
    var[i] = value;
  
  return true;
}

bool GridFromOptions::fetch(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return fetch(var, name.c_str(), lx, ly, lz);
}

/*******************************************************************************
 * GridDataGroup class
 * 
 * 
 *******************************************************************************/

void GridDataGroup::add(GridDataSource *source) {
  if(source == NULL)
    return;
  
  source_list.push_front(source);
}

GridDataSource* GridDataGroup::findSource(const char *name) {

  int msg_point = msg_stack.push("Finding source");
  for(std::list<GridDataSource*>::iterator it = source_list.begin(); 
      it != source_list.end(); it++) {
    
    // Query this source
    if((*it) != NULL)
      if((*it)->hasVar(name)) {
	msg_stack.pop(msg_point);
	return *it;
      }
  }
  msg_stack.pop(msg_point);
  return NULL;
}
