/**************************************************************************
 * Describes coordinate systems
 *
 **************************************************************************
 * Copyright 2014-2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#ifndef BOUT_COORDINATES_H
#define BOUT_COORDINATES_H

#include "bout/assert.hxx"
#include "bout/field_data.hxx"
#include <bout/bout_types.hxx>
#include <bout/build_defines.hxx>
#include <bout/christoffel_symbols.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/g_values.hxx>
#include <bout/metric_tensor.hxx>
#include <bout/paralleltransform.hxx>
#include <optional>

#include <array>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

class Mesh;
class YBoundary;
struct MetricNormaliser;

/*!
 * Represents a coordinate system, and associated operators
 */
class Coordinates {
public:
  using FieldMetric = bout::FieldMetric;

  /// Standard constructor from input
  Coordinates(Mesh* mesh, Options* options = nullptr);

  /// Constructor interpolating from another Coordinates object
  /// By default attempts to read staggered Coordinates from grid data source,
  /// interpolating from CELL_CENTRE if not present. Set
  /// force_interpolate_from_centre argument to true to always interpolate
  /// (useful if CELL_CENTRE Coordinates have been changed, so reading from file
  /// would not be correct).
  Coordinates(Mesh* mesh, Options* options, CELL_LOC loc, const Coordinates* coords_in,
              bool force_interpolate_from_centre = false);

  /// A constructor useful for testing purposes. To use it, inherit
  /// from Coordinates. If \p calculate_geometry is true (default),
  /// calculate the non-uniform variables, Christoffel symbols
  Coordinates(Mesh* mesh, FieldMetric dx, FieldMetric dy, FieldMetric dz, FieldMetric J,
              FieldMetric Bxy, FieldMetric g11, FieldMetric g22, FieldMetric g33,
              FieldMetric g12, FieldMetric g13, FieldMetric g23, FieldMetric g_11,
              FieldMetric g_22, FieldMetric g_33, FieldMetric g_12, FieldMetric g_13,
              FieldMetric g_23, FieldMetric ShiftTorsion, FieldMetric IntShiftTorsion);

  /// Add variables to \p output_options, for post-processing
  void outputVars(Options& output_options);

  ///< Mesh spacing in x, y and z
  const FieldMetric& dx() const { return dx_; }
  const FieldMetric& dy() const { return dy_; }
  const FieldMetric& dz() const { return dz_; }

  const BoutReal& dx(int x, int y, int z) const { return dx_(x, y, z); }
  const BoutReal& dy(int x, int y, int z) const { return dy_(x, y, z); }
  const BoutReal& dz(int x, int y, int z) const { return dz_(x, y, z); }

#if BOUT_USE_METRIC_3D
  const BoutReal* dx(int x, int y) const { return dx_(x, y); }
  const BoutReal* dy(int x, int y) const { return dy_(x, y); }
  const BoutReal* dz(int x, int y) const { return dz_(x, y); }
#else
  const BoutReal& dx(int x, int y) const { return dx_(x, y); }
  const BoutReal& dy(int x, int y) const { return dy_(x, y); }
  const BoutReal& dz(int x, int y) const { return dz_(x, y); }
#endif

  const BoutReal& IntShiftTorsion(int x, int y, int z) const {
    return IntShiftTorsion_(x, y, z);
  }

#if not(BOUT_USE_METRIC_3D)
  const BoutReal& IntShiftTorsion(int x, int y) const { return IntShiftTorsion_(x, y); }
#endif

  const BoutReal& J(int x, int y, int z) const { return J()(x, y, z); }

#if not(BOUT_USE_METRIC_3D)
  const BoutReal& J(int x, int y) const { return J()(x, y); }
#endif

  void setDx(FieldMetric dx, bool communicate = true);
  void setDy(FieldMetric dy, bool communicate = true);
  void setDz(FieldMetric dz, bool communicate = true);

  void setD1_dx(FieldMetric d1_dx) { d1_dx_ = std::move(d1_dx); }
  void setD1_dy(FieldMetric d1_dy) { d1_dy_ = std::move(d1_dy); }
  void setD1_dz(FieldMetric d1_dz) { d1_dz_ = std::move(d1_dz); }

  /// Length of the Z domain. Used for FFTs
  const Field2D& zlength() const;

  const BoutReal& zlength(int x, int y) const { return zlength()(x, y); }

  /// True if corrections for non-uniform mesh spacing should be included in operators
  bool non_uniform() const { return non_uniform_; }
  void setNon_uniform(bool non_uniform) { non_uniform_ = non_uniform; }

  /// 2nd-order correction for non-uniform meshes d/di(1/dx), d/di(1/dy) and d/di(1/dz)
  const FieldMetric& d1_dx() const { return d1_dx_; }
  const FieldMetric& d1_dy() const { return d1_dy_; }
  const FieldMetric& d1_dz() const { return d1_dz_; }

#if BOUT_USE_METRIC_3D
  const BoutReal& d1_dx(int x, int y, int z) const { return d1_dx_(x, y, z); }
  const BoutReal& d1_dy(int x, int y, int z) const { return d1_dy_(x, y, z); }
  const BoutReal& d1_dz(int x, int y, int z) const { return d1_dz_(x, y, z); }
#else
  const BoutReal& d1_dx(int x, int y) const { return d1_dx_(x, y); }
  const BoutReal& d1_dy(int x, int y) const { return d1_dy_(x, y); }
  const BoutReal& d1_dz(int x, int y) const { return d1_dz_(x, y); }
#endif

  /// Covariant metric tensor
  const FieldMetric& g_11() const { return covariantMetricTensor.g11(); }
  const FieldMetric& g_22() const { return covariantMetricTensor.g22(); }
  const FieldMetric& g_33() const { return covariantMetricTensor.g33(); }
  const FieldMetric& g_12() const { return covariantMetricTensor.g12(); }
  const FieldMetric& g_13() const { return covariantMetricTensor.g13(); }
  const FieldMetric& g_23() const { return covariantMetricTensor.g23(); }

  /// get g_22 at the cell faces;
  const FieldMetric& g_22_ylow() const;
  const FieldMetric& g_22_yhigh() const;
  // Cell Areas
  const FieldMetric& cell_area_xlow() const {
    if (_cell_area_xlow.has_value()) {
      return *_cell_area_xlow;
    }
    _compute_cell_area_x();
    ASSERT2(_cell_area_xlow.has_value());
    return *_cell_area_xlow;
  }
  const FieldMetric& cell_area_xhigh() const {
    if (_cell_area_xhigh.has_value()) {
      return *_cell_area_xhigh;
    }
    _compute_cell_area_x();
    ASSERT2(_cell_area_xhigh.has_value());
    return *_cell_area_xhigh;
  }
  const FieldMetric& cell_area_ylow() const {
    if (_cell_area_ylow.has_value()) {
      return *_cell_area_ylow;
    }
    _compute_cell_area_y();
    ASSERT2(_cell_area_ylow.has_value());
    return *_cell_area_ylow;
  }
  const FieldMetric& cell_area_yhigh() const {
    if (_cell_area_yhigh.has_value()) {
      return *_cell_area_yhigh;
    }
    _compute_cell_area_y();
    ASSERT2(_cell_area_yhigh.has_value());
    return *_cell_area_yhigh;
  }
  const FieldMetric& cell_area_zlow() const {
    if (_cell_area_zlow.has_value()) {
      return *_cell_area_zlow;
    }
    _compute_cell_area_z();
    ASSERT2(_cell_area_zlow.has_value());
    return *_cell_area_zlow;
  }
  const FieldMetric& cell_area_zhigh() const {
    if (_cell_area_zhigh.has_value()) {
      return *_cell_area_zhigh;
    }
    _compute_cell_area_z();
    ASSERT2(_cell_area_zhigh.has_value());
    return *_cell_area_zhigh;
  }
  FieldMetric& cell_area_xlow() {
    if (_cell_area_xlow.has_value()) {
      return *_cell_area_xlow;
    }
    _compute_cell_area_x();
    ASSERT2(_cell_area_xlow.has_value());
    return *_cell_area_xlow;
  }
  FieldMetric& cell_area_xhigh() {
    if (_cell_area_xhigh.has_value()) {
      return *_cell_area_xhigh;
    }
    _compute_cell_area_x();
    ASSERT2(_cell_area_xhigh.has_value());
    return *_cell_area_xhigh;
  }
  FieldMetric& cell_area_ylow() {
    if (_cell_area_ylow.has_value()) {
      return *_cell_area_ylow;
    }
    _compute_cell_area_y();
    ASSERT2(_cell_area_ylow.has_value());
    return *_cell_area_ylow;
  }
  FieldMetric& cell_area_yhigh() {
    if (_cell_area_yhigh.has_value()) {
      return *_cell_area_yhigh;
    }
    _compute_cell_area_y();
    ASSERT2(_cell_area_yhigh.has_value());
    return *_cell_area_yhigh;
  }
  FieldMetric& cell_area_zlow() {
    if (_cell_area_zlow.has_value()) {
      return *_cell_area_zlow;
    }
    _compute_cell_area_z();
    ASSERT2(_cell_area_zlow.has_value());
    return *_cell_area_zlow;
  }
  FieldMetric& cell_area_zhigh() {
    if (_cell_area_zhigh.has_value()) {
      return *_cell_area_zhigh;
    }
    _compute_cell_area_z();
    ASSERT2(_cell_area_zhigh.has_value());
    return *_cell_area_zhigh;
  }
  // Cell Volume
  const FieldMetric& cell_volume() const {
    if (_cell_volume.has_value()) {
      return *_cell_volume;
    }
    _compute_cell_volume();
    ASSERT2(_cell_volume.has_value());
    return *_cell_volume;
  }
  FieldMetric& cell_volume() {
    if (_cell_volume.has_value()) {
      return *_cell_volume;
    }
    _compute_cell_volume();
    ASSERT2(_cell_volume.has_value());
    return *_cell_volume;
  }

private:
  mutable std::optional<FieldMetric> _g_22_ylow, _g_22_yhigh;
  mutable std::optional<FieldMetric> _cell_area_xlow, _cell_area_xhigh;
  mutable std::optional<FieldMetric> _cell_area_ylow, _cell_area_yhigh;
  mutable std::optional<FieldMetric> _cell_area_zlow, _cell_area_zhigh;
  mutable std::optional<FieldMetric> _cell_volume;
  void _compute_cell_area_x() const;
  void _compute_cell_area_y() const;
  void _compute_cell_area_z() const;
  void _compute_cell_volume() const;

public:
  /// Contravariant metric tensor (g^{ij})
  const FieldMetric& g11() const { return contravariantMetricTensor.g11(); }
  const FieldMetric& g22() const { return contravariantMetricTensor.g22(); }
  const FieldMetric& g33() const { return contravariantMetricTensor.g33(); }
  const FieldMetric& g12() const { return contravariantMetricTensor.g12(); }
  const FieldMetric& g13() const { return contravariantMetricTensor.g13(); }
  const FieldMetric& g23() const { return contravariantMetricTensor.g23(); }

  /// Covariant metric tensor
  const BoutReal& g_11(int x, int y, int z) const {
    return covariantMetricTensor.g11(x, y, z);
  }
  const BoutReal& g_22(int x, int y, int z) const {
    return covariantMetricTensor.g22(x, y, z);
  }
  const BoutReal& g_33(int x, int y, int z) const {
    return covariantMetricTensor.g33(x, y, z);
  }
  const BoutReal& g_12(int x, int y, int z) const {
    return covariantMetricTensor.g12(x, y, z);
  }
  const BoutReal& g_13(int x, int y, int z) const {
    return covariantMetricTensor.g13(x, y, z);
  }
  const BoutReal& g_23(int x, int y, int z) const {
    return covariantMetricTensor.g23(x, y, z);
  }

#if not(BOUT_USE_METRIC_3D)
  const BoutReal& g_11(int x, int y) const { return covariantMetricTensor.g11(x, y); }
  const BoutReal& g_22(int x, int y) const { return covariantMetricTensor.g22(x, y); }
  const BoutReal& g_33(int x, int y) const { return covariantMetricTensor.g33(x, y); }
  const BoutReal& g_12(int x, int y) const { return covariantMetricTensor.g12(x, y); }
  const BoutReal& g_13(int x, int y) const { return covariantMetricTensor.g13(x, y); }
  const BoutReal& g_23(int x, int y) const { return covariantMetricTensor.g23(x, y); }
#endif

  /// Contravariant metric tensor (g^{ij})
  const BoutReal& g11(int x, int y, int z) const {
    return contravariantMetricTensor.g11(x, y, z);
  }
  const BoutReal& g22(int x, int y, int z) const {
    return contravariantMetricTensor.g22(x, y, z);
  }
  const BoutReal& g33(int x, int y, int z) const {
    return contravariantMetricTensor.g33(x, y, z);
  }
  const BoutReal& g12(int x, int y, int z) const {
    return contravariantMetricTensor.g12(x, y, z);
  }
  const BoutReal& g13(int x, int y, int z) const {
    return contravariantMetricTensor.g13(x, y, z);
  }
  const BoutReal& g23(int x, int y, int z) const {
    return contravariantMetricTensor.g23(x, y, z);
  }

#if not(BOUT_USE_METRIC_3D)
  const BoutReal& g11(int x, int y) const { return contravariantMetricTensor.g11(x, y); }
  const BoutReal& g22(int x, int y) const { return contravariantMetricTensor.g22(x, y); }
  const BoutReal& g33(int x, int y) const { return contravariantMetricTensor.g33(x, y); }
  const BoutReal& g12(int x, int y) const { return contravariantMetricTensor.g12(x, y); }
  const BoutReal& g13(int x, int y) const { return contravariantMetricTensor.g13(x, y); }
  const BoutReal& g23(int x, int y) const { return contravariantMetricTensor.g23(x, y); }
#endif

  const ContravariantMetricTensor& getContravariantMetricTensor() const {
    return contravariantMetricTensor;
  }

  const CovariantMetricTensor& getCovariantMetricTensor() const {
    return covariantMetricTensor;
  }

  void setContravariantMetricTensor(const ContravariantMetricTensor& metric_tensor,
                                    const std::string& region = "RGN_ALL",
                                    bool recalculate_staggered = true,
                                    bool force_interpolate_from_centre = false);

  void setCovariantMetricTensor(const CovariantMetricTensor& metric_tensor,
                                const std::string& region = "RGN_ALL",
                                bool recalculate_staggered = true,
                                bool force_interpolate_from_centre = false);

  void setMetricTensor(const ContravariantMetricTensor& contravariant_metric_tensor,
                       const CovariantMetricTensor& covariant_metric_tensor);

  void setMetricTensorJB(const ContravariantMetricTensor& contravariant_metric_tensor,
                         const CovariantMetricTensor& covariant_metric_tensor,
                         const FieldMetric& J, const FieldMetric& Bxy);

  void communicateMetricTensor();

  void communicateDz();

  void normaliseMetric(const MetricNormaliser& norm);

  ///< Coordinate system Jacobian, so volume of cell is J*dx*dy*dz
  const FieldMetric& J() const;

  ///< Magnitude of B = nabla z times nabla x
  const FieldMetric& Bxy() const { return Bxy_; }

  void setJ(const FieldMetric& J, bool communicate = true);

  void setBxy(FieldMetric Bxy, bool communicate = true);

  /// d pitch angle / dx. Needed for vector differentials (Curl)
  const FieldMetric& ShiftTorsion() const { return ShiftTorsion_; }

  ///< Integrated shear (I in BOUT notation)
  const FieldMetric& IntShiftTorsion() const { return IntShiftTorsion_; }

  void setIntShiftTorsion(FieldMetric IntShiftTorsion) {
    IntShiftTorsion_ = std::move(IntShiftTorsion);
  }

  ///////////////////////////////////////////////////////////
  // Parallel transforms
  ///////////////////////////////////////////////////////////

  /// Set the parallel (y) transform for this mesh.
  /// Mostly useful for tests.
  void setParallelTransform(std::unique_ptr<ParallelTransform> pt) {
    transform = std::move(pt);
  }

  bool hasParallelTransform() const { return transform != nullptr; }
  /// Return the parallel transform
  ParallelTransform& getParallelTransform() const {
    ASSERT1(hasParallelTransform());
    return *transform;
  }

  /// Christoffel symbol of the second kind (connection coefficients)
  const FieldMetric& G1_11() const { return christoffel_symbols().G1_11(); }
  const FieldMetric& G1_22() const { return christoffel_symbols().G1_22(); }
  const FieldMetric& G1_33() const { return christoffel_symbols().G1_33(); }
  const FieldMetric& G1_12() const { return christoffel_symbols().G1_12(); }
  const FieldMetric& G1_13() const { return christoffel_symbols().G1_13(); }
  const FieldMetric& G1_23() const { return christoffel_symbols().G1_23(); }
  const FieldMetric& G2_11() const { return christoffel_symbols().G2_11(); }
  const FieldMetric& G2_22() const { return christoffel_symbols().G2_22(); }
  const FieldMetric& G2_33() const { return christoffel_symbols().G2_33(); }
  const FieldMetric& G2_12() const { return christoffel_symbols().G2_12(); }
  const FieldMetric& G2_13() const { return christoffel_symbols().G2_13(); }
  const FieldMetric& G2_23() const { return christoffel_symbols().G2_23(); }
  const FieldMetric& G3_11() const { return christoffel_symbols().G3_11(); }
  const FieldMetric& G3_22() const { return christoffel_symbols().G3_22(); }
  const FieldMetric& G3_33() const { return christoffel_symbols().G3_33(); }
  const FieldMetric& G3_12() const { return christoffel_symbols().G3_12(); }
  const FieldMetric& G3_13() const { return christoffel_symbols().G3_13(); }
  const FieldMetric& G3_23() const { return christoffel_symbols().G3_23(); }

  const FieldMetric& G1() const { return g_values().G1(); }
  const FieldMetric& G2() const { return g_values().G2(); }
  const FieldMetric& G3() const { return g_values().G3(); }

  const BoutReal& G1(int x, int y, int z) const { return G1()(x, y, z); }
  const BoutReal& G2(int x, int y, int z) const { return G2()(x, y, z); }
  const BoutReal& G3(int x, int y, int z) const { return G3()(x, y, z); }

#if not(BOUT_USE_METRIC_3D)
  const BoutReal& G1(int x, int y) const { return G1()(x, y); }
  const BoutReal& G2(int x, int y) const { return G2()(x, y); }
  const BoutReal& G3(int x, int y) const { return G3()(x, y); }
#endif

  const FieldMetric& Grad2_par2_DDY_invSg(CELL_LOC outloc,
                                          const std::string& method) const;

  const FieldMetric& invSg() const;

  const ChristoffelSymbols& christoffel_symbols() const;

  GValues& g_values() const;

  void recalculateAndReset(bool recalculate_staggered,
                           bool force_interpolate_from_centre);

  FieldMetric recalculateJacobian() const;

  friend std::shared_ptr<YBoundary> getYBoundary(Coordinates* coords, YBndryType type);

private:
  std::shared_ptr<YBoundary> makeYBoundary(YBndryType type) const;
  int nz; // Size of mesh in Z. This is mesh->ngz-1
  Mesh* localmesh;
  Options* localoptions{nullptr};
  CELL_LOC location;

  /// True if corrections for non-uniform mesh spacing should be included in operators
  bool non_uniform_{};

  FieldMetric dx_, dy_, dz_; ///< Mesh spacing in x, y and z

  /// 2nd-order correction for non-uniform meshes d/di(1/dx), d/di(1/dy) and d/di(1/dz)
  FieldMetric d1_dx_, d1_dy_, d1_dz_;

  /// d pitch angle / dx. Needed for vector differentials (Curl)
  FieldMetric ShiftTorsion_;

  ///< Integrated shear (I in BOUT notation)
  FieldMetric IntShiftTorsion_;

  /// Handles calculation of yup and ydown
  std::unique_ptr<ParallelTransform> transform{nullptr};

  /// Cache variable for `zlength`. Invalidated when
  /// `Coordinates::recalculateAndReset` is called
  mutable std::unique_ptr<Field2D> zlength_cache{nullptr};

  /// Cache variable for Grad2_par2
  mutable std::map<std::string, std::unique_ptr<FieldMetric>> Grad2_par2_DDY_invSgCache;
  mutable std::unique_ptr<FieldMetric> invSgCache{nullptr};

  ContravariantMetricTensor contravariantMetricTensor;
  CovariantMetricTensor covariantMetricTensor;

  /// Christoffel symbol of the second kind (connection coefficients)
  mutable std::unique_ptr<ChristoffelSymbols> christoffel_symbols_cache{nullptr};

  /// `g_values` needs renaming, when we know what the name should be
  mutable std::unique_ptr<GValues> g_values_cache{nullptr};

  mutable std::unique_ptr<FieldMetric> jacobian_cache{nullptr};

  FieldMetric Bxy_; ///< Magnitude of B = nabla z times nabla x

  /// Set the parallel (y) transform from the options file.
  /// Used in the constructor to create the transform object.
  void setParallelTransform(Options* options);

  // check that covariant tensors are positive (if expected) and finite (always)
  void checkCovariant();
  // check that contravariant tensors are positive (if expected) and finite (always)
  void checkContravariant();
  void invalidateMetricCaches();
  void invalidateJacobianCaches();
  void invalidateCellGeometryCaches();
  void invalidateAccessorCache() const;

  mutable std::array<std::shared_ptr<YBoundary>, 3> ybndrys;

  FieldMetric recalculateBxy() const;

  /// Non-uniform meshes. Need to use DDX, DDY
  void correctionForNonUniformMeshes(bool force_interpolate_from_centre);

  void interpolateFromCoordinates(Options* options, const Coordinates* coords_in);

  /// Read quantities with given suffix from `Mesh`
  void readFromMesh(Options* options, const std::string& suffix);

  /// Read parallel slices of metric components from `Mesh`
  void readParallelMetricComponents();

protected:
  /// For testing purposes only; inherit and make this public
  void splitBxyParallelSlices();
};

namespace bout {
std::string parallelSliceFieldName(std::string_view field, int offset);
}

/// Represents a way to normalise the coordinate system
/// If a component returns nothing, no normalisation is performed.
/// Coordinate values are divided by the respective component from
/// MetricNormaliser, with the exception of the contravariant metric
/// tensor, which is multiplied by the normalisation factor.
struct MetricNormaliser {
  std::optional<BoutReal> g = std::nullopt;
  std::optional<BoutReal> g11 = std::nullopt;
  std::optional<BoutReal> g22 = std::nullopt;
  std::optional<BoutReal> g33 = std::nullopt;
  std::optional<BoutReal> g12 = std::nullopt;
  std::optional<BoutReal> g13 = std::nullopt;
  std::optional<BoutReal> g23 = std::nullopt;
  std::optional<BoutReal> dx = std::nullopt;
  std::optional<BoutReal> dy = std::nullopt;
  std::optional<BoutReal> dz = std::nullopt;
  std::optional<BoutReal> J = std::nullopt;
  std::optional<BoutReal> Bxy = std::nullopt;
};

#endif // BOUT_COORDINATES_H
