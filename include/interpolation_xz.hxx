/**************************************************************************
 * Copyright 2010-2020 B.D.Dudson, S.Farley, P. Hill, J.T. Omotani, J.T. Parker,
 * M.V.Umansky, X.Q.Xu
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

#ifndef __INTERP_XZ_H__
#define __INTERP_XZ_H__

#include "mask.hxx"

/// Interpolate a field onto a perturbed set of points
const Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                          const Field3D &delta_z);

const Field3D interpolate(const Field2D &f, const Field3D &delta_x,
                          const Field3D &delta_z);
const Field3D interpolate(const Field2D &f, const Field3D &delta_x);

class XZInterpolation {
protected:
  Options& options;

  Mesh* localmesh{nullptr};

  // 3D vector of points to skip (true -> skip this point)
  BoutMask skip_mask;

public:
  XZInterpolation(int y_offset = 0, Mesh* localmeshIn = nullptr, Options* opt = nullptr)
      : options(opt == nullptr ? Options::root()["xzinterpolation"] : *opt),
        localmesh(localmeshIn == nullptr ? bout::globals::mesh : localmeshIn),
        skip_mask(*localmesh, false), y_offset(y_offset) {}
  XZInterpolation(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr,
                  Options* opt = nullptr)
      : XZInterpolation(y_offset, mesh, opt) {
    skip_mask = mask;
  }
  virtual ~XZInterpolation() = default;

  virtual void calcWeights(const Field3D &delta_x, const Field3D &delta_z) = 0;
  virtual void calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                           const BoutMask &mask) = 0;

  virtual Field3D interpolate(const Field3D &f) const = 0;
  virtual Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                              const Field3D &delta_z) = 0;
  virtual Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                              const Field3D &delta_z, const BoutMask &mask) = 0;

  void setMask(const BoutMask &mask) { skip_mask = mask; }

  // Interpolate using the field at (x,y+y_offset,z), rather than (x,y,z)
  int y_offset;
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
        "XZInterpolation::getWeightsForYApproximation not implemented in this subclass");
  }
};

class XZHermiteSpline : public XZInterpolation {
protected:
  /// This is protected rather than private so that it can be
  /// extended and used by HermiteSplineMonotonic

  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  // Basis functions for cubic Hermite spline interpolation
  //    see http://en.wikipedia.org/wiki/Cubic_Hermite_spline
  // The h00 and h01 basis functions are applied to the function itself
  // and the h10 and h11 basis functions are applied to its derivative
  // along the interpolation direction.

  Field3D h00_x;
  Field3D h01_x;
  Field3D h10_x;
  Field3D h11_x;
  Field3D h00_z;
  Field3D h01_z;
  Field3D h10_z;
  Field3D h11_z;

public:
  XZHermiteSpline(Mesh *mesh = nullptr, Options* opt = nullptr)
      : XZHermiteSpline(0, mesh, opt) {}
  XZHermiteSpline(int y_offset = 0, Mesh *mesh = nullptr, Options* opt = nullptr);
  XZHermiteSpline(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr,
                  Options* opt = nullptr)
      : XZHermiteSpline(y_offset, mesh, opt) {
    skip_mask = mask;
  }

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z) override;
  void calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                   const BoutMask &mask) override;

  // Use precalculated weights
  Field3D interpolate(const Field3D &f) const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                      const Field3D &delta_z) override;
  Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z,
                      const BoutMask &mask) override;
  std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYApproximation(int i, int j, int k, int yoffset) override;
};


/// Monotonic Hermite spline interpolator
///
/// Similar to XZHermiteSpline, so uses most of the same code.
/// Forces the interpolated result to be in the range of the
/// neighbouring cell values. This prevents unphysical overshoots,
/// but also degrades accuracy near maxima and minima.
/// Perhaps should only impose near boundaries, since that is where
/// problems most obviously occur.
class XZMonotonicHermiteSpline : public XZHermiteSpline {
public:
  XZMonotonicHermiteSpline(Mesh *mesh = nullptr, Options* opt = nullptr)
      : XZHermiteSpline(0, mesh, opt) {}
  XZMonotonicHermiteSpline(int y_offset = 0, Mesh *mesh = nullptr, Options* opt = nullptr)
      : XZHermiteSpline(y_offset, mesh, opt) {}
  XZMonotonicHermiteSpline(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr,
      Options* opt = nullptr)
      : XZHermiteSpline(mask, y_offset, mesh, opt) {}

  using XZHermiteSpline::interpolate;
  /// Interpolate using precalculated weights.
  /// This function is called by the other interpolate functions
  /// in the base class XZHermiteSpline.
  Field3D interpolate(const Field3D &f) const override;
};

class XZLagrange4pt : public XZInterpolation {
  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  Field3D t_x, t_z;

public:
  XZLagrange4pt(Mesh *mesh = nullptr, Options* opt = nullptr)
      : XZLagrange4pt(0, mesh, opt) {}
  XZLagrange4pt(int y_offset = 0, Mesh *mesh = nullptr, Options* opt = nullptr);
  XZLagrange4pt(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr,
                Options* opt = nullptr)
      : XZLagrange4pt(y_offset, mesh, opt) {
    skip_mask = mask;
  }

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z) override;
  void calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                   const BoutMask &mask) override;

  // Use precalculated weights
  Field3D interpolate(const Field3D &f) const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                      const Field3D &delta_z) override;
  Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z,
                      const BoutMask &mask) override;
  BoutReal lagrange_4pt(BoutReal v2m, BoutReal vm, BoutReal vp, BoutReal v2p,
                        BoutReal offset) const;
  BoutReal lagrange_4pt(const BoutReal v[], BoutReal offset) const;
};

class XZBilinear : public XZInterpolation {
  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  Field3D w0, w1, w2, w3;

public:
  XZBilinear(Mesh *mesh = nullptr, Options* opt = nullptr) : XZBilinear(0, mesh, opt) {}
  XZBilinear(int y_offset = 0, Mesh *mesh = nullptr, Options* opt = nullptr);
  XZBilinear(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr,
             Options* opt = nullptr)
      : XZBilinear(y_offset, mesh, opt) {
    skip_mask = mask;
  }

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z) override;
  void calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                   const BoutMask &mask) override;

  // Use precalculated weights
  Field3D interpolate(const Field3D &f) const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                      const Field3D &delta_z) override;
  Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z,
                      const BoutMask &mask) override;
};

class XZInterpolationFactory
    : public Factory<XZInterpolation, XZInterpolationFactory,
                     std::function<std::unique_ptr<XZInterpolation>(Mesh*, Options*)>> {
public:
  static constexpr auto type_name = "XZInterpolation";
  static constexpr auto section_name = "xzinterpolation";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "hermitespline";

  using Factory::create;
  ReturnType create(Mesh* mesh = nullptr, Options* opt = nullptr) {
    return Factory::create(getType(nullptr), mesh, opt);
  }

  static void ensureRegistered();
};

template <class DerivedType>
class RegisterXZInterpolation {
public:
  RegisterXZInterpolation(const std::string& name) {
    XZInterpolationFactory::getInstance().add(
        name, [](Mesh* mesh, Options* opt) -> std::unique_ptr<XZInterpolation> {
          return std::make_unique<DerivedType>(mesh, opt);
        });
  }
};

#endif // __INTERP_XZ_H__
