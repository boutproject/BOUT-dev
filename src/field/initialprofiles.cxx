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

#include <bout/mesh.hxx>
#include <field2d.hxx>
#include <field3d.hxx>
#include <field_factory.hxx>
#include <globals.hxx>
#include <initialprofiles.hxx>
#include <msg_stack.hxx>

void initial_profile(const std::string& name, Field3D& var) {
  AUTO_TRACE();

  Mesh* localmesh = var.getMesh();

  Options* varOpts = Options::getRoot()->getSection(name);

  FieldFactory f(localmesh);

  std::string function;
  VAROPTION(varOpts, function, "0.0");

  var = f.create3D(function, varOpts, nullptr, var.getLocation());

  // Optionally scale the variable
  BoutReal scale;
  VAROPTION(varOpts, scale, 1.0);
  var *= scale;
}

void initial_profile(const std::string& name, Field2D& var) {
  AUTO_TRACE();

  Mesh* localmesh = var.getMesh();

  Options* varOpts = Options::getRoot()->getSection(name);

  FieldFactory f(localmesh);

  std::string function;
  VAROPTION(varOpts, function, "0.0");

  var = f.create2D(function, varOpts, nullptr, var.getLocation());

  // Optionally scale the variable
  BoutReal scale;
  VAROPTION(varOpts, scale, 1.0);
  var *= scale;
}

void initial_profile(const std::string& name, Vector2D& var) {
  AUTO_TRACE();

  if (var.covariant) {
    initial_profile(name + "_x", var.x);
    initial_profile(name + "_y", var.y);
    initial_profile(name + "_z", var.z);
  } else {
    initial_profile(name + "x", var.x);
    initial_profile(name + "y", var.y);
    initial_profile(name + "z", var.z);
  }
}

void initial_profile(const std::string& name, Vector3D& var) {
  AUTO_TRACE();

  if (var.covariant) {
    initial_profile(name + "_x", var.x);
    initial_profile(name + "_y", var.y);
    initial_profile(name + "_z", var.z);
  } else {
    initial_profile(name + "x", var.x);
    initial_profile(name + "y", var.y);
    initial_profile(name + "z", var.z);
  }
}
