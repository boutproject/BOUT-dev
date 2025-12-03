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

#ifndef BOUT_INTERP_XZ_H
#define BOUT_INTERP_XZ_H

#include <bout/bout_types.hxx>
#include <bout/generic_factory.hxx>
#include <bout/mask.hxx>

class Options;

/// Interpolate a field onto a perturbed set of points
Field3D interpolate(const Field3D& f, const Field3D& delta_x, const Field3D& delta_z);

Field3D interpolate(const Field2D& f, const Field3D& delta_x, const Field3D& delta_z);
Field3D interpolate(const Field2D& f, const Field3D& delta_x);

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
  XZInterpolation(const BoutMask& mask, int y_offset = 0, Mesh* mesh = nullptr)
      : XZInterpolation(y_offset, mesh) {
    setMask(mask);
  }
  XZInterpolation(const std::string& region_name, int y_offset = 0, Mesh* mesh = nullptr)
      : y_offset(y_offset), localmesh(mesh),
        region_id(localmesh->getRegionID(region_name)) {}
  XZInterpolation(const Region<Ind3D>& region, int y_offset = 0, Mesh* mesh = nullptr)
      : y_offset(y_offset), localmesh(mesh) {
    setRegion(region);
  }
  virtual ~XZInterpolation() = default;

  void setMask(const BoutMask& mask) { setRegion(regionFromMask(mask, localmesh)); }
  void setRegion(const std::string& region_name) {
    this->region_id = localmesh->getRegionID(region_name);
  }
  void setRegion(const std::unique_ptr<Region<Ind3D>> region) { setRegion(*region); }
  void setRegion(const Region<Ind3D>& region) {
    std::string name;
    int i = 0;
    do {
      name = fmt::format("unsec_reg_xz_interp_{:d}", i++);
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
    if (region == "" or region == "RGN_ALL") {
      return getRegion();
    }
    return localmesh->getRegion(
        localmesh->getCommonRegion(localmesh->getRegionID(region), region_id));
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

class XZInterpolationFactory
    : public Factory<XZInterpolation, XZInterpolationFactory, Mesh*, Options*> {
public:
  static constexpr auto type_name = "XZInterpolation";
  static constexpr auto section_name = "xzinterpolation";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "hermitespline";

  ReturnType create(Options* options = nullptr, Mesh* mesh = nullptr) const {
    return Factory::create(getType(options), mesh, options);
  }
  ReturnType create(const std::string& type, Options* options) const {
    return Factory::create(type, nullptr, options);
  }

  static void ensureRegistered();
};

template <class DerivedType>
using RegisterXZInterpolation = XZInterpolationFactory::RegisterInFactory<DerivedType>;

#endif // BOUT_INTERP_XZ_H
