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

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/coordinates.hxx"
#include "bout/region.hxx"
#include <bout/interpolation_xz.hxx>
#include <bout/mask.hxx>
#include <bout/parallel_boundary_region.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/unused.hxx>

#include <memory>
#include <string>
#include <vector>

class BoundaryRegionPar;
class FieldPerp;
class Field2D;
class Field3D;
class Options;

/// Field line map - contains the coefficients for interpolation
class FCIMap {
  /// Interpolation objects
  std::unique_ptr<XZInterpolation> interp;        // Cell centre
  std::unique_ptr<XZInterpolation> interp_corner; // Cell corner at (x+1, z+1)

  // The mesh this map was created on
  Mesh* map_mesh;

  /// Direction of map
  int offset_;

  /// region containing all points where the field line has not left the
  /// domain
  Region<Ind3D> region_no_boundary;
  /// If any of the integration area has left the domain
  BoutMask corner_boundary_mask;

public:
  FCIMap() = delete;
  FCIMap(Mesh& mesh, const Coordinates::FieldMetric& dy, Options& options, int offset,
         const std::shared_ptr<BoundaryRegionPar>& inner_boundary,
         const std::shared_ptr<BoundaryRegionPar>& outer_boundary, bool zperiodic);

  /// Direction of map
  int offset() const { return offset_; }

  Field3D interpolate(Field3D& f) const {
    ASSERT1(map_mesh == f.getMesh());
    return interp->interpolate(f);
  }

  Field3D integrate(Field3D& f) const;
};

/// Flux Coordinate Independent method for parallel derivatives
class FCITransform : public ParallelTransform {
public:
  FCITransform() = delete;
  FCITransform(Mesh& mesh, const Coordinates::FieldMetric& dy, bool zperiodic = true,
               Options* opt = nullptr);

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
