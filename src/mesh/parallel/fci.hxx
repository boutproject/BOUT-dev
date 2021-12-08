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
#include <interpolation_xz.hxx>
#include <mask.hxx>
#include <parallel_boundary_region.hxx>
#include <unused.hxx>

#include <memory>
#include <vector>


/// Field line map - contains the coefficients for interpolation
class FCIMap {
  /// Interpolation objects
  std::unique_ptr<XZInterpolation> interp;        // Cell centre
  std::unique_ptr<XZInterpolation> interp_corner; // Cell corner at (x+1, z+1)

public:
  FCIMap() = delete;
  FCIMap(Mesh& mesh, const Coordinates::FieldMetric& dy, Options& options, int offset,
         BoundaryRegionPar* inner_boundary, BoundaryRegionPar* outer_boundary,
         bool zperiodic);

  // The mesh this map was created on
  Mesh& map_mesh;

  /// Direction of map
  const int offset;

  /// boundary mask - has the field line left the domain
  BoutMask boundary_mask;
  /// If any of the integration area has left the domain
  BoutMask corner_boundary_mask;
  
  Field3D interpolate(Field3D& f) const {
    ASSERT1(&map_mesh == f.getMesh());
    return interp->interpolate(f);
  }

  Field3D integrate(Field3D &f) const;
};


/// Flux Coordinate Independent method for parallel derivatives
class FCITransform : public ParallelTransform {
public:
  FCITransform() = delete;
  FCITransform(Mesh& mesh, const Coordinates::FieldMetric& dy, bool zperiodic = true,
               Options* opt = nullptr)
      : ParallelTransform(mesh, opt) {

    // check the coordinate system used for the grid data source
    FCITransform::checkInputGrid();

    auto forward_boundary_xin =
        new BoundaryRegionPar("FCI_forward", BNDRY_PAR_FWD_XIN, +1, &mesh);
    auto backward_boundary_xin =
        new BoundaryRegionPar("FCI_backward", BNDRY_PAR_BKWD_XIN, -1, &mesh);
    auto forward_boundary_xout =
        new BoundaryRegionPar("FCI_forward", BNDRY_PAR_FWD_XOUT, +1, &mesh);
    auto backward_boundary_xout =
        new BoundaryRegionPar("FCI_backward", BNDRY_PAR_BKWD_XOUT, -1, &mesh);

    // Add the boundary region to the mesh's vector of parallel boundaries
    mesh.addBoundaryPar(forward_boundary_xin);
    mesh.addBoundaryPar(backward_boundary_xin);
    mesh.addBoundaryPar(forward_boundary_xout);
    mesh.addBoundaryPar(backward_boundary_xout);

    field_line_maps.reserve(mesh.ystart * 2);
    for (int offset = 1; offset < mesh.ystart + 1; ++offset) {
      field_line_maps.emplace_back(mesh, dy, options, offset, forward_boundary_xin,
                                   forward_boundary_xout, zperiodic);
      field_line_maps.emplace_back(mesh, dy, options, -offset, backward_boundary_xin,
                                   backward_boundary_xout, zperiodic);
    }
  }

  void calcParallelSlices(Field3D &f) override;
  
  void integrateParallelSlices(Field3D &f) override;

  Field3D toFieldAligned(const Field3D& UNUSED(f),
                         const std::string& UNUSED(region) = "RGN_ALL") override {
    throw BoutException("FCI method cannot transform into field aligned grid");
  }
  FieldPerp toFieldAligned(const FieldPerp& UNUSED(f),
                           const std::string& UNUSED(region) = "RGN_ALL") override {
    throw BoutException("FCI method cannot transform into field aligned grid");
  }

  Field3D fromFieldAligned(const Field3D& UNUSED(f),
                           const std::string& UNUSED(region) = "RGN_ALL") override {
    throw BoutException("FCI method cannot transform from field aligned grid");
  }
  FieldPerp fromFieldAligned(const FieldPerp& UNUSED(f),
                             const std::string& UNUSED(region) = "RGN_ALL") override {
    throw BoutException("FCI method cannot transform from field aligned grid");
  }

  bool canToFromFieldAligned() override { return false; }

  bool requiresTwistShift(bool UNUSED(twist_shift_enabled), MAYBE_UNUSED(YDirectionType ytype)) override {
    // No Field3Ds require twist-shift, because they cannot be field-aligned
    ASSERT1(ytype == YDirectionType::Standard);

    return false;
  }

protected:
  void checkInputGrid() override;


private:
  /// FCI maps for each of the parallel slices
  std::vector<FCIMap> field_line_maps;
};

#endif // __FCITRANSFORM_H__
