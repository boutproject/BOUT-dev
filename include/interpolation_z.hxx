/**************************************************************************
 * Copyright 2020 P. Hill, J.T. Omotani, J.T. Parker
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

#ifndef __INTERP_Z_H__
#define __INTERP_Z_H__

#include "bout/generic_factory.hxx"
#include "bout/paralleltransform.hxx"
#include "bout/region.hxx"

class ZInterpolation {
protected:
  Mesh* localmesh{nullptr};

  Region<Ind3D> region;

public:
  explicit ZInterpolation(int y_offset = 0, Mesh* mesh = nullptr,
                          Region<Ind3D> region_in = {});
  virtual ~ZInterpolation() = default;

  virtual void calcWeights(const Field3D& delta_z) = 0;

  virtual Field3D interpolate(const Field3D& f,
                              const std::string& region_str = "DEFAULT") const = 0;
  virtual Field3D interpolate(const Field3D& f, const Field3D& delta_z,
                              const std::string& region_str = "DEFAULT") = 0;

  void setRegion(Region<Ind3D> new_region) {
    region = new_region;
  }

  virtual std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYUpApproximation(int i, int j, int k) const {
    return getWeightsForYApproximation(i, j, k, 1);
  }
  virtual std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYDownApproximation(int i, int j, int k) const {
    return getWeightsForYApproximation(i, j, k, -1);
  }
  virtual std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYApproximation(int UNUSED(i), int UNUSED(j), int UNUSED(k),
                              int UNUSED(yoffset)) const {
    throw BoutException(
        "ZInterpolation::getWeightsForYApproximation not implemented in this subclass");
  }

  // Interpolate using the field at (x,y+y_offset,z), rather than (x,y,z)
  const int y_offset;
};

class ZInterpolationFactory
    : public Factory<ZInterpolation, ZInterpolationFactory, int, Mesh*, Region<Ind3D>> {
public:
  static constexpr auto type_name = "ZInterpolation";
  static constexpr auto section_name = "zinterpolation";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "hermitespline";

  using Factory::create;
  ReturnType create(Options* options, int y_offset = 0, Mesh* mesh = nullptr,
                    Region<Ind3D> region_in = {}) {
    return Factory::create(options, y_offset, mesh, region_in);
  }
  ReturnType create(int y_offset = 0, Mesh* mesh = nullptr,
                    Region<Ind3D> region_in = {}) {
    return Factory::create(getType(nullptr), y_offset, mesh, region_in);
  }

  static void ensureRegistered();
};

template <class DerivedType>
using RegisterZInterpolation = ZInterpolationFactory::RegisterInFactory<DerivedType>;

class ZHermiteSpline : public ZInterpolation {
public:
  explicit ZHermiteSpline(int y_offset = 0, Mesh* mesh = nullptr,
                          Region<Ind3D> region_in = {});

  void calcWeights(const Field3D& delta_z) override;

  // Use precalculated weights
  Field3D interpolate(const Field3D& f,
                      const std::string& region_str = "DEFAULT") const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D& f, const Field3D& delta_z,
                      const std::string& region_str = "DEFAULT") override;
  std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYApproximation(int i, int j, int k, int yoffset) const override;

private:
  const std::string fz_region;

  Array<Ind3D> k_corner; // z-index of left grid point

  // Basis functions for cubic Hermite spline interpolation
  //    see http://en.wikipedia.org/wiki/Cubic_Hermite_spline
  // The h00 and h01 basis functions are applied to the function itself
  // and the h10 and h11 basis functions are applied to its derivative
  // along the interpolation direction.

  Field3D h00;
  Field3D h01;
  Field3D h10;
  Field3D h11;
};

#endif // __INTERP_Z_H__
