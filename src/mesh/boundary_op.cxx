/***************************************************************************
 * Copyright 2018 B.D. Dudson, J.T. Omotani
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

#include <boundary_op.hxx>
#include <bout/constants.hxx>
#include <bout/mesh.hxx>

void BoundaryOp::apply_ddt(Field2D &f) {
  Field2D *dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < f.getNz(); z++)
      (*dt)(bndry->x, bndry->y, z) = 0.; // Set time derivative to zero
}

void BoundaryOp::apply_ddt(Field3D &f) {
  Field3D *dt = f.timeDeriv();
  for (bndry->first(); !bndry->isDone(); bndry->next())
    for (int z = 0; z < f.getNz(); z++)
      (*dt)(bndry->x, bndry->y, z) = 0.; // Set time derivative to zero
}
