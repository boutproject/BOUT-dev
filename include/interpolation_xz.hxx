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

#define USE_NEW_WEIGHTS 1
#if BOUT_HAS_PETSC
#define HS_USE_PETSC 1
#endif

#ifdef HS_USE_PETSC
#include "bout/petsclib.hxx"
#endif

class Options;

/// Interpolate a field onto a perturbed set of points
const Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                          const Field3D &delta_z);

const Field3D interpolate(const Field2D &f, const Field3D &delta_x,
                          const Field3D &delta_z);
const Field3D interpolate(const Field2D &f, const Field3D &delta_x);

class XZInterpolation {
public:
  int y_offset;

protected:
  Mesh* localmesh{nullptr};

  int region_id{-1};

public:
  XZInterpolation(int y_offset = 0, Mesh* localmeshIn = nullptr)
      : y_offset(y_offset),
        localmesh(localmeshIn == nullptr ? bout::globals::mesh : localmeshIn) {}
  XZInterpolation(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : XZInterpolation(y_offset, mesh) {
    setMask(mask);
  }
  XZInterpolation(const std::string& region_name, int y_offset = 0, Mesh* mesh = nullptr)
    : y_offset(y_offset), localmesh(mesh), region_id(localmesh->getRegionID(region_name)) {}
  XZInterpolation(const Region<Ind3D>& region, int y_offset = 0,
                  Mesh* mesh = nullptr)
      : y_offset(y_offset), localmesh(mesh){
    setRegion(region);
  }
  virtual ~XZInterpolation() = default;

  void setMask(const BoutMask& mask) {
    setRegion(regionFromMask(mask, localmesh));
  }
  void setRegion(const std::string& region_name) {
    this->region_id = localmesh->getRegionID(region_name);
  }
  void setRegion(const Region<Ind3D>& region) {
    std::string name;
    int i=0;
    do {
      name = fmt::format("unsec_reg_xz_interp_{:d}",i++);
    } while (localmesh->hasRegion3D(name));
    localmesh->addRegion(name, region);
    this->region_id = localmesh->getRegionID(name);
  }
  const Region<Ind3D>& getRegion() const {
    ASSERT2(region_id != -1);
    return localmesh->getRegion(region_id);
  }
  const Region<Ind3D>& getRegion(const std::string& region) const {
    if (region_id == -1) {
      return localmesh->getRegion(region);
    }
    if (region == "" or region == "RGN_ALL"){
      return getRegion();
    }
    return localmesh->getRegion(localmesh->getCommonRegion(localmesh->getRegionID(region), region_id));
  }
  virtual void calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                           const std::string& region = "RGN_NOBNDRY") = 0;
  virtual void calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                           const BoutMask& mask,
                           const std::string& region = "RGN_NOBNDRY") = 0;

  virtual Field3D interpolate(const Field3D& f,
                              const std::string& region = "RGN_NOBNDRY") const = 0;
  virtual Field3D interpolate(const Field3D& f, const Field3D& delta_x,
                              const Field3D& delta_z,
                              const std::string& region = "RGN_NOBNDRY") = 0;
  virtual Field3D interpolate(const Field3D& f, const Field3D& delta_x,
                              const Field3D& delta_z, const BoutMask& mask,
                              const std::string& region = "RGN_NOBNDRY") = 0;

  // Interpolate using the field at (x,y+y_offset,z), rather than (x,y,z)
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

  Tensor<SpecificInd<IND_TYPE::IND_3D>> i_corner; // index of bottom-left grid point
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

  std::vector<Field3D> newWeights;

#if HS_USE_PETSC
  PetscLib* petsclib;
  bool isInit{false};
  Mat petscWeights;
  Vec rhs, result;
#endif

public:
  XZHermiteSpline(Mesh *mesh = nullptr)
      : XZHermiteSpline(0, mesh) {}
  XZHermiteSpline(int y_offset = 0, Mesh *mesh = nullptr);
  XZHermiteSpline(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : XZHermiteSpline(y_offset, mesh) {
    setRegion(regionFromMask(mask, localmesh));
  }
  ~XZHermiteSpline() {
#if HS_USE_PETSC
    if (isInit) {
      MatDestroy(&petscWeights);
      VecDestroy(&rhs);
      VecDestroy(&result);
      isInit = false;
      delete petsclib;
    }
#endif
  }

  void calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                   const std::string& region = "RGN_NOBNDRY") override;
  void calcWeights(const Field3D& delta_x, const Field3D& delta_z, const BoutMask& mask,
                   const std::string& region = "RGN_NOBNDRY") override;

  // Use precalculated weights
  Field3D interpolate(const Field3D& f,
                      const std::string& region = "RGN_NOBNDRY") const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D& f, const Field3D& delta_x, const Field3D& delta_z,
                      const std::string& region = "RGN_NOBNDRY") override;
  Field3D interpolate(const Field3D& f, const Field3D& delta_x, const Field3D& delta_z,
                      const BoutMask& mask,
                      const std::string& region = "RGN_NOBNDRY") override;
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
  XZMonotonicHermiteSpline(Mesh *mesh = nullptr)
      : XZHermiteSpline(0, mesh) {}
  XZMonotonicHermiteSpline(int y_offset = 0, Mesh *mesh = nullptr)
      : XZHermiteSpline(y_offset, mesh) {}
  XZMonotonicHermiteSpline(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : XZHermiteSpline(mask, y_offset, mesh) {}

  using XZHermiteSpline::interpolate;
  /// Interpolate using precalculated weights.
  /// This function is called by the other interpolate functions
  /// in the base class XZHermiteSpline.
  Field3D interpolate(const Field3D& f,
                      const std::string& region = "RGN_NOBNDRY") const override;
};

class XZLagrange4pt : public XZInterpolation {
  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  Field3D t_x, t_z;

public:
  XZLagrange4pt(Mesh *mesh = nullptr)
      : XZLagrange4pt(0, mesh) {}
  XZLagrange4pt(int y_offset = 0, Mesh *mesh = nullptr);
  XZLagrange4pt(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : XZLagrange4pt(y_offset, mesh) {
    setRegion(regionFromMask(mask, localmesh));
  }

  void calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                   const std::string& region = "RGN_NOBNDRY") override;
  void calcWeights(const Field3D& delta_x, const Field3D& delta_z, const BoutMask& mask,
                   const std::string& region = "RGN_NOBNDRY") override;

  // Use precalculated weights
  Field3D interpolate(const Field3D& f,
                      const std::string& region = "RGN_NOBNDRY") const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D& f, const Field3D& delta_x, const Field3D& delta_z,
                      const std::string& region = "RGN_NOBNDRY") override;
  Field3D interpolate(const Field3D& f, const Field3D& delta_x, const Field3D& delta_z,
                      const BoutMask& mask,
                      const std::string& region = "RGN_NOBNDRY") override;
  BoutReal lagrange_4pt(BoutReal v2m, BoutReal vm, BoutReal vp, BoutReal v2p,
                        BoutReal offset) const;
  BoutReal lagrange_4pt(const BoutReal v[], BoutReal offset) const;
};

class XZBilinear : public XZInterpolation {
  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  Field3D w0, w1, w2, w3;

public:
  XZBilinear(Mesh *mesh = nullptr) : XZBilinear(0, mesh) {}
  XZBilinear(int y_offset = 0, Mesh *mesh = nullptr);
  XZBilinear(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : XZBilinear(y_offset, mesh) {
    setRegion(regionFromMask(mask, localmesh));
  }

  void calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                   const std::string& region = "RGN_NOBNDRY") override;
  void calcWeights(const Field3D& delta_x, const Field3D& delta_z, const BoutMask& mask,
                   const std::string& region = "RGN_NOBNDRY") override;

  // Use precalculated weights
  Field3D interpolate(const Field3D& f,
                      const std::string& region = "RGN_NOBNDRY") const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D& f, const Field3D& delta_x, const Field3D& delta_z,
                      const std::string& region = "RGN_NOBNDRY") override;
  Field3D interpolate(const Field3D& f, const Field3D& delta_x, const Field3D& delta_z,
                      const BoutMask& mask,
                      const std::string& region = "RGN_NOBNDRY") override;
};

class XZInterpolationFactory
    : public Factory<XZInterpolation, XZInterpolationFactory, Mesh*> {
public:
  static constexpr auto type_name = "XZInterpolation";
  static constexpr auto section_name = "xzinterpolation";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "hermitespline";

  ReturnType create(Options* options = nullptr, Mesh* mesh = nullptr) const {
    return Factory::create(getType(options), mesh);
  }
  ReturnType create(const std::string& type, MAYBE_UNUSED(Options* options)) const {
    return Factory::create(type, nullptr);
  }

  static void ensureRegistered();
};

template <class DerivedType>
using RegisterXZInterpolation = XZInterpolationFactory::RegisterInFactory<DerivedType>;

#endif // __INTERP_XZ_H__
