/**************************************************************************
 * Shifted metric parallel derivatives, implemented with interpolation
 * instead of FFT
 *
 **************************************************************************
 * Copyright 2018 B.D.Dudson, P. Hill, J. Omotani, J. Parker
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

#ifndef __SHIFTEDINTERP_H__
#define __SHIFTEDINTERP_H__

#include <bout/interpolation_z.hxx>
#include <bout/paralleltransform.hxx>
#include <cstddef>

/*!
 * Shifted metric method
 * Each Y location is shifted in Z with respect to its neighbours
 * so that the grid is orthogonal in X-Z, but requires interpolation
 * to calculate the values of points along field-lines.
 *
 * In this implementation the interpolation is done using ZInterpolation objects
 */
class ShiftedMetricInterp : public ParallelTransform {
public:
  ShiftedMetricInterp() = delete;
  ShiftedMetricInterp(Mesh& mesh, CELL_LOC location_in, Field2D zShift_in,
                      BoutReal zlength_in, Options* opt = nullptr);

  /*!
   * Calculates the yup() and ydown() fields of f
   * by interpolating f through a toroidal shift angle
   */
  void calcParallelSlices(Field3D& f) override;

  /*!
   * Uses interpolation of f through a toroidal shift angle to align the grid
   * points with the y coordinate (along magnetic field usually).
   *
   * Note that the returned field will no longer be orthogonal in X-Z, and the
   * metric tensor will need to be changed if X derivatives are used.
   */
  Field3D toFieldAligned(const Field3D& f,
                         const std::string& region = "RGN_ALL") override;
  FieldPerp toFieldAligned(const FieldPerp& UNUSED(f),
                           const std::string& UNUSED(region) = "RGN_ALL") override {
    throw BoutException("Not implemented yet");
  }

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  Field3D fromFieldAligned(const Field3D& f,
                           const std::string& region = "RGN_ALL") override;
  FieldPerp fromFieldAligned(const FieldPerp& UNUSED(f),
                             const std::string& UNUSED(region) = "RGN_ALL") override {
    throw BoutException("Not implemented yet");
  }

  bool canToFromFieldAligned() const override { return true; }

  std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYApproximation(int i, int j, int k, int y_offset) override {
#if CHECK > 0
    if (y_offset == 0) {
      throw BoutException("Cannot get ShiftedMetricInterp parallel slice at zero offset");
    }
#endif
    const auto index = static_cast<std::size_t>(std::abs(y_offset) - 1)
                       + ((y_offset > 0) ? yup_index : ydown_index);

    return parallel_slice_interpolators[index]->getWeightsForYApproximation(i, j, k,
                                                                            y_offset);
  }

  bool requiresTwistShift(bool twist_shift_enabled, YDirectionType ytype) override {
    // Twist-shift only if field-aligned
    if (ytype == YDirectionType::Aligned and not twist_shift_enabled) {
      throw BoutException("'twistshift = true' is required to communicate field-aligned "
                          "Field3Ds when using ShiftedMetric.");
    }
    return ytype == YDirectionType::Aligned;
  }

protected:
  void checkInputGrid() override;

private:
  CELL_LOC location{CELL_CENTRE};

  /// This is the shift in toroidal angle (z) which takes a point from
  /// X-Z orthogonal to field-aligned along Y.
  Field2D zShift;

  /// Length of the z-domain in radians
  BoutReal zlength{0.};

  /// Cache of interpolators for the parallel slices. Slices are stored
  /// in the following order:
  ///     {+1, ..., +n, -1, ..., -n}
  /// parallel_slice_interpolator[i] stores interpolator for slice i+1
  /// parallel_slice_interpolator[n + i] stores offset -(i+1)
  /// where i goes from 0 to (n-1), with n the number of y guard cells
  std::vector<std::unique_ptr<ZInterpolation>> parallel_slice_interpolators;

  /// ZInterpolation objects for shifting to and from field-aligned coordinates
  std::unique_ptr<ZInterpolation> interp_to_aligned, interp_from_aligned;

  static constexpr std::size_t yup_index = 0;
  const std::size_t ydown_index;
};

#endif // __SHIFTEDINTERP_H__
