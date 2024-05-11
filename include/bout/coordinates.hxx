/**************************************************************************
 * Describes coordinate systems
 *
 * ChangeLog
 * =========
 * 
 * 2014-11-10 Ben Dudson <bd512@york.ac.uk>
 *    * Created by separating metric from Mesh
 *
 * 
 **************************************************************************
 * Copyright 2014 B.D.Dudson
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

#ifndef BOUT_COORDINATES_H
#define BOUT_COORDINATES_H

#include "christoffel_symbols.hxx"
#include "g_values.hxx"
#include "bout/metric_tensor.hxx"
#include "bout/paralleltransform.hxx"

class Mesh;

/*!
 * Represents a coordinate system, and associated operators
 */
class Coordinates {
public:
#if BOUT_USE_METRIC_3D
  using FieldMetric = Field3D;
#else
  using FieldMetric = Field2D;
#endif

  /// Constructor interpolating from another Coordinates object
  /// By default attempts to read staggered Coordinates from grid data source,
  /// interpolating from CELL_CENTRE if not present. Set
  /// force_interpolate_from_centre argument to true to always interpolate
  /// (useful if CELL_CENTRE Coordinates have been changed, so reading from file
  /// would not be correct).
  explicit Coordinates(Mesh* mesh, Options* options = nullptr, CELL_LOC loc = CELL_CENTRE,
                       const Coordinates* coords_in = nullptr,
                       bool force_interpolate_from_centre = false);

  /// A constructor useful for testing purposes. To use it, inherit
  /// from Coordinates. If \p calculate_geometry is true (default),
  /// calculate the non-uniform variables, Christoffel symbols
  Coordinates(Mesh* mesh, FieldMetric dx, FieldMetric dy, FieldMetric dz,
              [[maybe_unused]] const FieldMetric& J, FieldMetric Bxy,
              const FieldMetric& g11, const FieldMetric& g22, const FieldMetric& g33,
              const FieldMetric& g12, const FieldMetric& g13, const FieldMetric& g23,
              const FieldMetric& g_11, const FieldMetric& g_22, const FieldMetric& g_33,
              const FieldMetric& g_12, const FieldMetric& g_13, const FieldMetric& g_23,
              FieldMetric ShiftTorsion, FieldMetric IntShiftTorsion);

  /// Add variables to \p output_options, for post-processing
  void outputVars(Options& output_options);

  ///< Mesh spacing in x, y and z
  const FieldMetric& dx() const { return dx_; }
  const FieldMetric& dy() const { return dy_; }
  const FieldMetric& dz() const { return dz_; }

  void setDx(FieldMetric dx) { dx_ = std::move(dx); }
  void setDy(FieldMetric dy) { dy_ = std::move(dy); }
  void setDz(FieldMetric dz) { dz_ = std::move(dz); }

#if BOUT_USE_METRIC_3D
  void setDy(BoutReal value, int x, int y, int z) { dy_(x, y, z) = value; }
#else
  void setDy(BoutReal value, int x, int y) { dy_(x, y) = value; }
#endif

  void setD1_dx(FieldMetric d1_dx) { d1_dx_ = std::move(d1_dx); }
  void setD1_dy(FieldMetric d1_dy) { d1_dy_ = std::move(d1_dy); }
  void setD1_dz(FieldMetric d1_dz) { d1_dz_ = std::move(d1_dz); }

  /// Length of the Z domain. Used for FFTs
  const Field2D& zlength() const;

  /// True if corrections for non-uniform mesh spacing should be included in operators
  bool non_uniform() const { return non_uniform_; }
  void setNon_uniform(bool non_uniform) { non_uniform_ = non_uniform; }

  /// 2nd-order correction for non-uniform meshes d/di(1/dx), d/di(1/dy) and d/di(1/dz)
  const FieldMetric& d1_dx() const { return d1_dx_; }
  const FieldMetric& d1_dy() const { return d1_dy_; }
  const FieldMetric& d1_dz() const { return d1_dz_; }

  /// Covariant metric tensor
  const MetricTensor::FieldMetric& g_11() const { return covariantMetricTensor.g11(); }
  const MetricTensor::FieldMetric& g_22() const { return covariantMetricTensor.g22(); }
  const MetricTensor::FieldMetric& g_33() const { return covariantMetricTensor.g33(); }
  const MetricTensor::FieldMetric& g_12() const { return covariantMetricTensor.g12(); }
  const MetricTensor::FieldMetric& g_13() const { return covariantMetricTensor.g13(); }
  const MetricTensor::FieldMetric& g_23() const { return covariantMetricTensor.g23(); }

  /// Contravariant metric tensor (g^{ij})
  const MetricTensor::FieldMetric& g11() const { return contravariantMetricTensor.g11(); }
  const MetricTensor::FieldMetric& g22() const { return contravariantMetricTensor.g22(); }
  const MetricTensor::FieldMetric& g33() const { return contravariantMetricTensor.g33(); }
  const MetricTensor::FieldMetric& g12() const { return contravariantMetricTensor.g12(); }
  const MetricTensor::FieldMetric& g13() const { return contravariantMetricTensor.g13(); }
  const MetricTensor::FieldMetric& g23() const { return contravariantMetricTensor.g23(); }

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

  void setMetricTensor(ContravariantMetricTensor contravariant_metric_tensor,
                                      CovariantMetricTensor covariant_metric_tensor);

  ///< Coordinate system Jacobian, so volume of cell is J*dx*dy*dz
  const FieldMetric& J() const;

  ///< Magnitude of B = nabla z times nabla x
  const FieldMetric& Bxy() const { return Bxy_; }

  void setJ(const FieldMetric& J);

#if BOUT_USE_METRIC_3D
  void setJ(BoutReal value, int x, int y, int z);
#else
  void setJ(BoutReal value, int x, int y);
#endif

  void setBxy(FieldMetric Bxy);

  /// d pitch angle / dx. Needed for vector differentials (Curl)
  const FieldMetric& ShiftTorsion() const { return ShiftTorsion_; }

  ///< Integrated shear (I in BOUT notation)
  const FieldMetric& IntShiftTorsion() const { return IntShiftTorsion_; }

  void setIntShiftTorsion(FieldMetric IntShiftTorsion) {
    IntShiftTorsion_ = std::move(IntShiftTorsion);
  }

  int communicateAndCheckMeshSpacing() const;

  ///////////////////////////////////////////////////////////
  // Parallel transforms
  ///////////////////////////////////////////////////////////

  /// Set the parallel (y) transform for this mesh.
  /// Mostly useful for tests.
  void setParallelTransform(std::unique_ptr<ParallelTransform> pt) {
    transform = std::move(pt);
  }

  /// Return the parallel transform
  ParallelTransform& getParallelTransform() {
    ASSERT1(transform != nullptr);
    return *transform;
  }

  ///////////////////////////////////////////////////////////
  // Operators
  ///////////////////////////////////////////////////////////

  FieldMetric DDX(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT",
                  const std::string& region = "RGN_NOBNDRY") const;

  FieldMetric DDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT",
                  const std::string& region = "RGN_NOBNDRY") const;

  FieldMetric DDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT",
                  const std::string& region = "RGN_NOBNDRY") const;

  Field3D DDX(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field3D DDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field3D DDZ(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  /// Gradient along magnetic field  b.Grad(f)
  FieldMetric Grad_par(const Field2D& var, CELL_LOC outloc = CELL_DEFAULT,
                       const std::string& method = "DEFAULT");

  Field3D Grad_par(const Field3D& var, CELL_LOC outloc = CELL_DEFAULT,
                   const std::string& method = "DEFAULT");

  /// Advection along magnetic field V*b.Grad(f)
  FieldMetric Vpar_Grad_par(const Field2D& v, const Field2D& f,
                            CELL_LOC outloc = CELL_DEFAULT,
                            const std::string& method = "DEFAULT");

  Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                        CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");

  /// Divergence along magnetic field  Div(b*f) = B.Grad(f/B)
  FieldMetric Div_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                      const std::string& method = "DEFAULT");

  Field3D Div_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT");

  // Second derivative along magnetic field
  FieldMetric Grad2_par2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT");

  Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                     const std::string& method = "DEFAULT");
  // Perpendicular Laplacian operator, using only X-Z derivatives
  // NOTE: This might be better bundled with the Laplacian inversion code
  // since it makes use of the same coefficients and FFT routines
  FieldMetric Delp2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
  Field3D Delp2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
  FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);

  // Full parallel Laplacian operator on scalar field
  // Laplace_par(f) = Div( b (b dot Grad(f)) )
  FieldMetric Laplace_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT);
  Field3D Laplace_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT);

  // Full Laplacian operator on scalar field
  FieldMetric Laplace(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                      const std::string& dfdy_boundary_conditions = "free_o3",
                      const std::string& dfdy_dy_region = "");
  Field3D Laplace(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& dfdy_boundary_conditions = "free_o3",
                  const std::string& dfdy_dy_region = "");

  // Full perpendicular Laplacian, in form of inverse of Laplacian operator in LaplaceXY
  // solver
  Field2D Laplace_perpXY(const Field2D& A, const Field2D& f) const;

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

  void setG1(const FieldMetric& G1) const { g_values().setG1(G1); }
  void setG2(const FieldMetric& G2) const { g_values().setG2(G2); }
  void setG3(const FieldMetric& G3) const { g_values().setG3(G3); }

  const FieldMetric& Grad2_par2_DDY_invSg(CELL_LOC outloc,
                                          const std::string& method) const;

  const FieldMetric& invSg() const;

  ChristoffelSymbols& christoffel_symbols() const;

  GValues& g_values() const;

  void recalculateAndReset(bool recalculate_staggered,
                           bool force_interpolate_from_centre);

  FieldMetric recalculateJacobian() const;

  template <typename T, typename... Ts>
  void communicate(T& t, Ts... ts) const;

private:
  int nz; // Size of mesh in Z. This is mesh->ngz-1
  Mesh* localmesh;
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

  void applyToContravariantMetricTensor(
      const std::function<const FieldMetric(const FieldMetric)>& function);

  void applyToCovariantMetricTensor(
      const std::function<const FieldMetric(const FieldMetric)>& function);

  void applyToChristoffelSymbols(
      const std::function<const FieldMetric(const FieldMetric)>& function) const;

  mutable std::unique_ptr<FieldMetric> jacobian_cache{nullptr};

  FieldMetric Bxy_; ///< Magnitude of B = nabla z times nabla x

  void invalidateAndRecalculateCachedVariables();

  /// Set the parallel (y) transform from the options file.
  /// Used in the constructor to create the transform object.
  void setParallelTransform(Options* options);

  // check that covariant tensors are positive (if expected) and finite (always)
  void checkCovariant();
  // check that contravariant tensors are positive (if expected) and finite (always)
  void checkContravariant();

  FieldMetric getAtLocOrUnaligned(Mesh* mesh, const std::string& name,
                                  BoutReal default_value = 0.,
                                  const std::string& suffix = "",
                                  CELL_LOC cell_location = CELL_CENTRE);

  FieldMetric getUnaligned(const std::string& name, BoutReal default_value);

  void communicateChristoffelSymbolTerms() const;
  void extrapolateChristoffelSymbols();

  void communicateGValues() const;
  void extrapolateGValues();

  FieldMetric recalculateBxy() const;

  /// Non-uniform meshes. Need to use DDX, DDY
  void correctionForNonUniformMeshes(bool force_interpolate_from_centre);

  void interpolateFieldsFromOtherCoordinates(Options* options,
                                             const Coordinates* coords_in);

  void setBoundaryCells(Options* options, const std::string& suffix);

  FieldMetric getDzFromOptionsFile(Mesh* mesh, const std::string& suffix) const;

  void fixZShiftGuards(Field2D& zShift) const;

  Field2D interpolateAndExtrapolate(const Field2D& f, CELL_LOC location,
                                    bool extrapolate_x, bool extrapolate_y,
                                    bool no_extra_interpolate,
                                    ParallelTransform* UNUSED_pt,
                                    const std::string& region) const;

#if BOUT_USE_METRIC_3D
  Field3D interpolateAndExtrapolate(const Field3D& f_, CELL_LOC location,
                                    bool extrapolate_x, bool extrapolate_y,
                                    bool no_extra_interpolate,
                                    ParallelTransform* pt_) const;

#endif // BOUT_USE_METRIC_3D
};

/*
/// Standard coordinate system for tokamak simulations
class TokamakCoordinates : public Coordinates {
public:
  TokamakCoordinates(Mesh *mesh) : Coordinates(mesh) {
    
  }
private:
  
};
*/

#endif // BOUT_COORDINATES_H
