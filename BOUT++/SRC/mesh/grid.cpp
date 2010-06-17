/*!
 * \file grid.cpp
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

#include "globals.h"
#include "mesh_topology.h"
#include "utils.h"

#include "dataformat.h" // Abstract data class

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "fft.h" // For reading 3D variables
#include "dcomplex.h"

#include "bout_types.h"

/*******************************************************************************
 * GridFile class
 * 
 * This is a fairly thin wrapper around the Datafile class. Only reading
 * routines are required though. Takes care of opening and closing file
 * then just passes read requests through.
 *******************************************************************************/

/// Default constructor
GridFile::GridFile()
{
  file = NULL; ///< Signal that the file is invalid
  isOpen = false;
}

/// Constructor passing a file format and name
GridFile::GridFile(DataFormat *format, const char *gridfilename)
{
  setFile(format, gridfilename);
}

void GridFile::setFile(DataFormat *format, const char *gridfilename)
{
  file = format;
  filename = string(gridfilename);
  isOpen = false;
}

bool GridFile::hasVar(const char *name)
{
  if(file == NULL)
    return false;

  /// Try to get the size of the variable
  vector<int> s = getSize(name);
  
  /// Test if the variable has zero size
  return s.size() != 0;
}

/// Return the size of variable
vector<int> GridFile::getSize(const char *name)
{
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

bool GridFile::setOrigin(int x, int y, int z)
{
  if(file == NULL)
    return false;
  
  return file->setOrigin(x,y,z);
}

bool GridFile::fetch(int *var, const char *name, int lx, int ly, int lz)
{
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

bool GridFile::fetch(int *var, const string &name, int lx, int ly, int lz)
{
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

bool GridFile::fetch(real *var, const char *name, int lx, int ly, int lz)
{
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

bool GridFile::fetch(real *var, const string &name, int lx, int ly, int lz)
{
  if(file == NULL)
    return false;
  if(!file->is_valid())
    return false;
  
  return file->read(var, name, lx, ly, lz);
}

void GridFile::open(const char *name)
{
  file->openr(filename);
  if(!file->is_valid()) {
    output << "\tERROR: Could not open file " << filename << endl;
  }
  isOpen = true;
}

void GridFile::close()
{
  file->close();
  isOpen = false;
}
