/**************************************************************************
 * Sets initial profiles
 *
 * ChangeLog
 * =========
 * 
 * 2011-02-12 Ben Dudson <bd512@york.ac.uk>
 *    * Changed to use new options system. For now the structure of the
 *      options is the same, but this could be modified more easily in future
 *
 * 2010-05-12 Ben Dudson <bd512@york.ac.uk>
 *    
 *    * Changed random numbers to use a hash of the parameters
 *      so that the phase doesn't vary with number of processors or grid size
 *      User can vary phase to give a different random sequence
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

#include <globals.hxx>
#include <initialprofiles.hxx>
#include <boutexception.hxx>
#include <field_factory.hxx>
#include <output.hxx>
#include <bout/constants.hxx>
#include <msg_stack.hxx>

#include <math.h>
#include <string.h>
#include <stdlib.h>

void initial_profile(const string &name, Field3D &var) {
  TRACE("initial_profile(string, Field3D)");
  
  CELL_LOC loc = CELL_DEFAULT;
  if (mesh->StaggerGrids) {
    loc = var.getLocation();
  }
  
  // Get the section for this specific variable 
  Options *varOpts = Options::getRoot()->getSection(name);
  
  // Use FieldFactory to generate values
    
  FieldFactory f(mesh);

  string function;
  VAROPTION(varOpts, function, "0.0");
  
  // Create a 3D variable
  var = f.create3D(function, varOpts, NULL, loc);

  // Optionally scale the variable
  BoutReal scale;
  VAROPTION(varOpts, scale, 1.0);
  var *= scale;
}

// For 2D variables almost identical, just no z dependence
void initial_profile(const string &name, Field2D &var) {
  
  CELL_LOC loc = var.getLocation();
  
  // Get the section for this variable
  Options *varOpts = Options::getRoot()->getSection(name);
  output << name;
  
  // Use FieldFactory to generate values
    
  FieldFactory f(mesh);

  string function;
  VAROPTION(varOpts, function, "0.0");

  var = f.create2D(function, varOpts, NULL, loc);

  // Optionally scale the variable
  BoutReal scale;
  VAROPTION(varOpts, scale, 1.0);
  var *= scale;
}

void initial_profile(const string &name, Vector2D &var) {
  if(var.covariant) {
    initial_profile(name + "_x", var.x);
    initial_profile(name + "_y", var.y);
    initial_profile(name + "_z", var.z);
  }else {
    initial_profile(name + "x", var.x);
    initial_profile(name + "y", var.y);
    initial_profile(name + "z", var.z);
  }
}

void initial_profile(const string &name, Vector3D &var) {
  if(var.covariant) {
    initial_profile(name + "_x", var.x);
    initial_profile(name + "_y", var.y);
    initial_profile(name + "_z", var.z);
  }else {
    initial_profile(name + "x", var.x);
    initial_profile(name + "y", var.y);
    initial_profile(name + "z", var.z);
  }
}

/**************************************************************************
 * Routines to generate profiles 
 **************************************************************************/

/// Generate a 3D field with a given Z oscillation
/*!
  @param[in] n      Mode number. Note that this is mode-number in the domain
  @param[in] phase  Phase shift in units of pi
*/
const Field3D genZMode(int n, BoutReal phase) {
  Field3D result;

  result.allocate();
  
  for(int jz=0;jz<mesh->LocalNz;jz++) {
    BoutReal val = sin(phase*PI +  TWOPI * ((BoutReal) jz)/ ((BoutReal) mesh->LocalNz) );
    for(int jx=0;jx<mesh->LocalNx;jx++)
      for(int jy=0;jy<mesh->LocalNy;jy++)
	result(jx,jy,jz) = val;
  }
  return result;
}
