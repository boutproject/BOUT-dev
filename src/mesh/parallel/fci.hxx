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

#ifndef BOUT_FCITRANSFORM_H
#define BOUT_FCITRANSFORM_H

#include <bout/interpolation_xz.hxx>
#include <bout/mask.hxx>
#include <bout/parallel_boundary_region.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/unused.hxx>

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
         const std::shared_ptr<BoundaryRegionPar>& inner_boundary,
         const std::shared_ptr<BoundaryRegionPar>& outer_boundary, bool zperiodic);

  // The mesh this map was created on
  Mesh& map_mesh;

  /// Direction of map
  const int offset;

  /// region containing all points where the field line has not left the
  /// domain
  Region<Ind3D> region_no_boundary;
  /// If any of the integration area has left the domain
  BoutMask corner_boundary_mask;

  Field3D interpolate(Field3D& f) const {
    ASSERT1(&map_mesh == f.getMesh());
    return interp->interpolate(f);
  }

  Field3D integrate(Field3D& f) const;
};

/// Flux Coordinate Independent method for parallel derivatives
class FCITransform : public ParallelTransform {
public:
  FCITransform() = delete;
  FCITransform(Mesh& mesh, const Coordinates::FieldMetric& dy, bool zperiodic = true,
               Options* opt = nullptr)
    : ParallelTransform(mesh, opt), R{&mesh}, Z{&mesh} {

    // check the coordinate system used for the grid data source
    FCITransform::checkInputGrid();

    // Real-space coordinates of grid cells
    mesh.get(R, "R", 0.0, false);
    mesh.get(Z, "Z", 0.0, false);

    auto forward_boundary_xin =
        std::make_shared<BoundaryRegionPar>("FCI_forward", BNDRY_PAR_FWD_XIN, +1, &mesh);
    auto backward_boundary_xin = std::make_shared<BoundaryRegionPar>(
        "FCI_backward", BNDRY_PAR_BKWD_XIN, -1, &mesh);
    auto forward_boundary_xout =
        std::make_shared<BoundaryRegionPar>("FCI_forward", BNDRY_PAR_FWD_XOUT, +1, &mesh);
    auto backward_boundary_xout = std::make_shared<BoundaryRegionPar>(
        "FCI_backward", BNDRY_PAR_BKWD_XOUT, -1, &mesh);

    // Add the boundary region to the mesh's vector of parallel boundaries
    mesh.addBoundaryPar(forward_boundary_xin, BoundaryParType::xin_fwd);
    mesh.addBoundaryPar(backward_boundary_xin, BoundaryParType::xin_bwd);
    mesh.addBoundaryPar(forward_boundary_xout, BoundaryParType::xout_fwd);
    mesh.addBoundaryPar(backward_boundary_xout, BoundaryParType::xout_bwd);

    field_line_maps.reserve(mesh.ystart * 2);
    for (int offset = 1; offset < mesh.ystart + 1; ++offset) {
      field_line_maps.emplace_back(mesh, dy, options, offset, forward_boundary_xin,
                                   forward_boundary_xout, zperiodic);
      field_line_maps.emplace_back(mesh, dy, options, -offset, backward_boundary_xin,
                                   backward_boundary_xout, zperiodic);
    }
    ASSERT0(mesh.ystart == 1);
    std::shared_ptr<BoundaryRegionPar> bndries[]{
        forward_boundary_xin, forward_boundary_xout, backward_boundary_xin,
        backward_boundary_xout};
    for (auto& bndry : bndries) {
      for (const auto& bndry2 : bndries) {
        if (bndry->dir == bndry2->dir) {
          continue;
        }
        for (bndry->first(); !bndry->isDone(); bndry->next()) {
          if (bndry2->contains(*bndry)) {
            bndry->setValid(0);
          }
        }
      }
    }
  }

  void calcParallelSlices(Field3D& f) override;

  void integrateParallelSlices(Field3D& f) override;

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

  bool canToFromFieldAligned() const override { return false; }

  /// Save mesh variables to output
  /// If R and Z(x,y,z) coordinates are in the input then these are saved to output.
  void outputVars(Options& output_options) override;

  bool requiresTwistShift(bool UNUSED(twist_shift_enabled),
                          [[maybe_unused]] YDirectionType ytype) override {
    // No Field3Ds require twist-shift, because they cannot be field-aligned
    ASSERT1(ytype == YDirectionType::Standard);

    return false;
  }

protected:
  void checkInputGrid() override;

private:
  /// FCI maps for each of the parallel slices
  std::vector<FCIMap> field_line_maps;

  /// Real-space coordinates of grid points
  Field3D R, Z;
};

#endif // BOUT_FCITRANSFORM_H
