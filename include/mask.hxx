/**************************************************************************
 * Copyright 2015 B.D.Dudson, P. Hill
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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
 **************************************************************************/

#ifndef __MASK_H__
#define __MASK_H__

#include <vector>

#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "globals.hxx"

/**
 * 3D array of bools to mask Field3Ds
 *
 * Wrapper around a 3D vector of bools to enable masking of
 * Field3Ds. Masking is not automatic, but can be accomplished by
 *
 *     // Create mask the size of mesh with all false values
 *     BoutMask mask(mesh, false);
 *     // Set an index to true to skip this index
 *     mask(3, 4, 5) = true;
 *     // Iterate over field
 *     for (auto index : field) {
 *       // Skip any indices which are set to true in the mask
 *       if (mask(index.x, index.y, index.z)) continue;
 *       ...
 *     }
 */
class BoutMask {
  // Typedef for internal container
  typedef std::vector<std::vector<std::vector<bool>>> vec3bool;

  // Dimensions of mask
  int nx;
  int ny;
  int nz;
  // Internal data
  vec3bool mask;
public:
  BoutMask(int nx, int ny, int nz, bool value=false) :
    nx(nx), ny(ny), nz(nz),
    mask(nx, std::vector<std::vector<bool>>
         (ny, std::vector<bool>
          (nz, value))) {}
  BoutMask(Mesh& mesh, bool value=false) :
    BoutMask(mesh.LocalNx, mesh.LocalNy, mesh.LocalNz, value) {}
  // Default constructor uses global mesh
  BoutMask() : BoutMask(*mesh) {}

  // Assignment from bool
  BoutMask& operator=(bool value) {
    mask = vec3bool(nx, std::vector<std::vector<bool>>
                    (ny, std::vector<bool>
                     (nz, value)));
    return *this;
  }

  // vector<bool> is weird nonsense and doesn't actually store bools.
  // Hence we need to return the vector<bool>::reference member type
  inline std::vector<bool>::reference operator()(int jx,int jy,int jz) {
#if CHECK > 2
    // Perform bounds checking
    if((jx < 0) || (jx >= nx) ||
       (jy < 0) || (jy >= ny) ||
       (jz < 0) || (jz >= nz))
      throw BoutException("BoutMask: (%d, %d, %d) operator out of bounds (%d, %d, %d)",
                          jx, jy, jz, nx, ny, nz);
#endif
    return mask[jx][jy][jz];
  }
  inline std::vector<bool>::const_reference operator()(int jx,int jy,int jz) const {
#if CHECK > 2
    // Perform bounds checking
    if((jx < 0) || (jx >= nx) ||
       (jy < 0) || (jy >= ny) ||
       (jz < 0) || (jz >= nz))
      throw BoutException("BoutMask: (%d, %d, %d) operator out of bounds (%d, %d, %d)",
                          jx, jy, jz, nx, ny, nz);
#endif
    return mask[jx][jy][jz];
  }

};




#endif //__MASK_H__
