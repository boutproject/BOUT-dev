/**************************************************************************
 * Flux-coordinate Independent parallel derivatives
 *
 **************************************************************************
 * Copyright 2014 B.D.Dudson, P. Hill
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

#ifndef __FCITRANSFORM_H__
#define __FCITRANSFORM_H__

#include <bout/paralleltransform.hxx>
#include <interpolation.hxx>
#include <mask.hxx>
#include <parallel_boundary_region.hxx>
#include <unused.hxx>

#include <vector>

/*!
 * Field line map - contains the coefficients for interpolation
 */
class FCIMap {
  /// Interpolation object
  Interpolation *interp;        // Cell centre
  Interpolation *interp_corner; // Cell corner at (x+1, z+1)

public:
  FCIMap() = delete;
  FCIMap(Mesh& mesh, int offset, BoundaryRegionPar* boundary, bool zperiodic);

  /// Direction of map
  const int offset;

  BoutMask boundary_mask;      /**< boundary mask - has the field line left the domain */
  BoutMask corner_boundary_mask; ///< If any of the integration area has left the domain
  
  Field3D y_prime;             /**< distance to intersection with boundary */

  Field3D interpolate(Field3D &f) const { return interp->interpolate(f); }

  Field3D integrate(Field3D &f) const;
};

/*!
 * Flux Coordinate Independent method for parallel derivatives
 */
class FCITransform : public ParallelTransform {
public:
  FCITransform() = delete;
  FCITransform(Mesh& mesh, bool zperiodic = true) : mesh(mesh), zperiodic(zperiodic) {
    auto forward_boundary = new BoundaryRegionPar("FCI_forward", BNDRY_PAR_FWD, +1, &mesh);
    auto backward_boundary = new BoundaryRegionPar("FCI_backward", BNDRY_PAR_BKWD, -1, &mesh);

    // Add the boundary region to the mesh's vector of parallel boundaries
    mesh.addBoundaryPar(forward_boundary);
    mesh.addBoundaryPar(backward_boundary);

    field_line_maps.reserve(mesh.ystart * 2);
    for (int offset = 1; offset < mesh.ystart + 1; ++offset) {
      field_line_maps.emplace_back(mesh, offset, forward_boundary, zperiodic);
      field_line_maps.emplace_back(mesh, -offset, backward_boundary, zperiodic);
    }
  }

  void calcYUpDown(Field3D &f) override;
  
  void integrateYUpDown(Field3D &f) override;
  
  const Field3D toFieldAligned(const Field3D &UNUSED(f)) override {
    throw BoutException("FCI method cannot transform into field aligned grid");
  }

  const Field3D fromFieldAligned(const Field3D &UNUSED(f)) override {
    throw BoutException("FCI method cannot transform into field aligned grid");
  }

  bool canToFromFieldAligned() override{
    return false;
  }
private:
  Mesh& mesh;

  /// FCI maps for field lines in +ve y
  std::vector<FCIMap> field_line_maps;

  /// Is the z-direction periodic?
  bool zperiodic;
};

#endif // __FCITRANSFORM_H__
