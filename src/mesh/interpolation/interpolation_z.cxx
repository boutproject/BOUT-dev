/**************************************************************************
 * Copyright 2020 B.D.Dudson, P. Hill, J. Omotani, J. Parker
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

#include <bout/globals.hxx>
#include <bout/interpolation_z.hxx>
#include <bout/mesh.hxx>
#include <bout/region.hxx>

ZInterpolation::ZInterpolation(int y_offset, Mesh* mesh, const Region<Ind3D>& region_in)
    : localmesh(mesh == nullptr ? bout::globals::mesh : mesh),
      region(region_in.size() == 0 ? localmesh->getRegion3D("RGN_NOYZ")
                                   : region_in),
      y_offset(y_offset) {}
