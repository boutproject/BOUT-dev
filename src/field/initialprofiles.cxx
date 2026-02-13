/**************************************************************************
 * Sets initial profiles
 *
 **************************************************************************
 * Copyright 2010 - 2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/field_factory.hxx>
#include <bout/globals.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/mesh.hxx>

void initial_profile(const std::string& name, Field3D& var) {

  Mesh* localmesh = var.getMesh();

  Options* varOpts = Options::getRoot()->getSection(name);

  std::string function;
  VAROPTION(varOpts, function, "0.0");

  var = FieldFactory::get()->create3D(function, varOpts, localmesh, var.getLocation());

  // Optionally scale the variable
  BoutReal scale;
  VAROPTION(varOpts, scale, 1.0);
  var *= scale;
}

void initial_profile(const std::string& name, Field2D& var) {

  Mesh* localmesh = var.getMesh();

  Options* varOpts = Options::getRoot()->getSection(name);

  std::string function;
  VAROPTION(varOpts, function, "0.0");

  var = FieldFactory::get()->create2D(function, varOpts, localmesh, var.getLocation());

  // Optionally scale the variable
  BoutReal scale;
  VAROPTION(varOpts, scale, 1.0);
  var *= scale;
}

void initial_profile(const std::string& name, Vector2D& var) {

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
