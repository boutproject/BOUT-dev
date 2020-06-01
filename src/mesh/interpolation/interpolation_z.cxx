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

#include <bout/mesh.hxx>
#include <interpolation_z.hxx>

ZInterpolation::ZInterpolation(int y_offset, Mesh* mesh, Region<Ind3D>* region_ptr)
    : localmesh(mesh == nullptr ? bout::globals::mesh : mesh),
      region(region_ptr), y_offset(y_offset) {
  if (region == nullptr) {
    // Construct region that skips calculating interpolation in y-boundary regions that
    // should be filled by boundary conditions

    region = std::make_unique<Region<Ind3D>>(localmesh->getRegion3D("RGN_NOBNDRY"));

    const int ny = localmesh->LocalNy;
    const int nz = localmesh->LocalNz;
    auto mask_region = Region<Ind3D>(0, 0, 0, 0, 0, 0, ny, nz);
    if (y_offset > 0) {
      for (auto it = localmesh->iterateBndryUpperY(); not it.isDone(); it.next()) {
        mask_region.getUnion(Region<Ind3D>(it.ind, it.ind, localmesh->yend - y_offset + 1,
                                           localmesh->yend, localmesh->zstart,
                                           localmesh->zend, ny, nz));
      }
    }
    else if (y_offset < 0) {
      for (auto it = localmesh->iterateBndryLowerY(); not it.isDone(); it.next()) {
        mask_region.getUnion(Region<Ind3D>(it.ind, it.ind,
                                           localmesh->ystart, localmesh->ystart - y_offset - 1,
                                           localmesh->zstart, localmesh->zend,
                                           ny, nz));
      }
    }

    region->mask(mask_region);
  }
}
