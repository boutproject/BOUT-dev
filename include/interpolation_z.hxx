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
#include "mask.hxx"

class ZInterpolation {
protected:
  Mesh* localmesh{nullptr};

  // 3D vector of points to skip (true -> skip this point)
  BoutMask skip_mask;

public:
  ZInterpolation(int y_offset = 0, Mesh* localmeshIn = nullptr)
      : localmesh(localmeshIn == nullptr ? bout::globals::mesh : localmeshIn),
        skip_mask(*localmesh, false), y_offset(y_offset) {}
  ZInterpolation(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : ZInterpolation(y_offset, mesh) {
    skip_mask = mask;
  }
  virtual ~ZInterpolation() = default;

  virtual void calcWeights(const Field3D &delta_z,
                           const std::string& region = "RGN_NOBNDRY") = 0;
  virtual void calcWeights(const Field3D &delta_z, const BoutMask &mask,
                           const std::string& region = "RGN_NOBNDRY") = 0;

  virtual Field3D interpolate(const Field3D &f,
                              const std::string& region = "RGN_NOBNDRY") const = 0;
  virtual Field3D interpolate(const Field3D &f, const Field3D &delta_z,
                              const std::string& region = "RGN_NOBNDRY") = 0;
  virtual Field3D interpolate(const Field3D &f, const Field3D &delta_z,
                              const BoutMask &mask,
                              const std::string& region = "RGN_NOBNDRY") = 0;

  void setMask(const BoutMask &mask) { skip_mask = mask; }

  /// Interpolate using the field at (x,y+y_offset,z), rather than (x,y,z)
  void setYOffset(int offset) { y_offset = offset; }

  virtual std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYUpApproximation(int i, int j, int k) {
    return getWeightsForYApproximation(i, j, k, 1);
  }
  virtual std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYDownApproximation(int i, int j, int k) {
    return getWeightsForYApproximation(i, j, k, -1);
  }
  virtual std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYApproximation(int UNUSED(i), int UNUSED(j), int UNUSED(k),
                              int UNUSED(yoffset)) {
    throw BoutException(
        "ZInterpolation::getWeightsForYApproximation not implemented in this subclass");
  }

protected:
  // Interpolate using the field at (x,y+y_offset,z), rather than (x,y,z)
  int y_offset;
};

class ZInterpolationFactory
    : public Factory<ZInterpolation, ZInterpolationFactory,
                             std::function<std::unique_ptr<ZInterpolation>(Mesh*)>> {
public:
  static constexpr auto type_name = "ZInterpolation";
  static constexpr auto section_name = "zinterpolation";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "hermitespline";

  using Factory::create;
  ReturnType create(Mesh* mesh = nullptr) {
    return Factory::create(getType(nullptr), mesh);
  }

  static void ensureRegistered();
};

template <class DerivedType>
class RegisterZInterpolation {
public:
  RegisterZInterpolation(const std::string& name) {
    ZInterpolationFactory::getInstance().add(
        name, [](Mesh* mesh) -> std::unique_ptr<ZInterpolation> {
          return std::make_unique<DerivedType>(mesh);
        });
  }
};

class ZHermiteSpline : public ZInterpolation {
public:
  ZHermiteSpline(Mesh *mesh = nullptr) : ZHermiteSpline(0, mesh) {}
  ZHermiteSpline(int y_offset = 0, Mesh *mesh = nullptr);
  ZHermiteSpline(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : ZHermiteSpline(y_offset, mesh) {
    skip_mask = mask;
  }

  /// Callback function for ZInterpolationFactory
  static ZInterpolation *CreateZHermiteSpline(Mesh *mesh) {
    return new ZHermiteSpline(mesh);
  }

  void calcWeights(const Field3D& delta_z,
                   const std::string& region = "RGN_NOBNDRY") override;
  void calcWeights(const Field3D& delta_z, const BoutMask &mask,
                   const std::string& region = "RGN_NOBNDRY") override;

  // Use precalculated weights
  Field3D interpolate(const Field3D &f,
                      const std::string& region = "RGN_NOBNDRY") const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D &f, const Field3D &delta_z,
                      const std::string& region = "RGN_NOBNDRY") override;
  Field3D interpolate(const Field3D &f, const Field3D &delta_z, const BoutMask &mask,
                      const std::string& region = "RGN_NOBNDRY") override;
  std::vector<ParallelTransform::PositionsAndWeights>
    getWeightsForYApproximation(int i, int j, int k, int yoffset);

private:
  Tensor<int> k_corner; // z-index of left grid point

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
