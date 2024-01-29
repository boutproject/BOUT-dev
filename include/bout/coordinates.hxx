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
#include "differential_operators.hxx"
#include "bout/metricTensor.hxx"
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
  Coordinates(Mesh* mesh, FieldMetric dx, FieldMetric dy, FieldMetric dz, FieldMetric J,
              FieldMetric Bxy, const FieldMetric& g11, const FieldMetric& g22,
              const FieldMetric& g33, const FieldMetric& g12, const FieldMetric& g13,
              const FieldMetric& g23, const FieldMetric& g_11, const FieldMetric& g_22,
              const FieldMetric& g_33, const FieldMetric& g_12, const FieldMetric& g_13,
              const FieldMetric& g_23, FieldMetric ShiftTorsion,
              FieldMetric IntShiftTorsion);

  /// Add variables to \p output_options, for post-processing
  void outputVars(Options& output_options);

  ///< Mesh spacing in x, y and z
  const FieldMetric& dx() const;
  const FieldMetric& dy() const;
  const FieldMetric& dz() const;

  void setDx(FieldMetric dx);
  void setDy(FieldMetric dy);
  void setDz(FieldMetric dz);

  void setDy(BoutReal value, int x, int y);

  void setD1_dx(FieldMetric d1_dx);
  void setD1_dy(FieldMetric d1_dy);
  void setD1_dz(FieldMetric d1_dz);

  /// Length of the Z domain. Used for FFTs
  const Field2D& zlength() const;

  /// True if corrections for non-uniform mesh spacing should be included in operators
  bool non_uniform() const;
  void setNon_uniform(bool non_uniform);

  /// 2nd-order correction for non-uniform meshes d/di(1/dx), d/di(1/dy) and d/di(1/dz)
  const FieldMetric& d1_dx() const;
  const FieldMetric& d1_dy() const;
  const FieldMetric& d1_dz() const;

  /// Covariant metric tensor
  const FieldMetric& g_11() const;
  const FieldMetric& g_22() const;
  const FieldMetric& g_33() const;
  const FieldMetric& g_12() const;
  const FieldMetric& g_13() const;
  const FieldMetric& g_23() const;

  /// Contravariant metric tensor (g^{ij})
  const FieldMetric& g11() const;
  const FieldMetric& g22() const;
  const FieldMetric& g33() const;
  const FieldMetric& g12() const;
  const FieldMetric& g13() const;
  const FieldMetric& g23() const;

  const MetricTensor& getContravariantMetricTensor() const;
  const MetricTensor& getCovariantMetricTensor() const;

  void setContravariantMetricTensor(const MetricTensor& metric_tensor,
                                    const std::string& region = "RGN_ALL",
                                    bool recalculate_staggered = true,
                                    bool force_interpolate_from_centre = false);

  void setCovariantMetricTensor(const MetricTensor& metric_tensor,
                                const std::string& region = "RGN_ALL",
                                bool recalculate_staggered = true,
                                bool force_interpolate_from_centre = false);

  ///< Coordinate system Jacobian, so volume of cell is J*dx*dy*dz
  const FieldMetric& J() const;

  ///< Magnitude of B = nabla z times nabla x
  const FieldMetric& Bxy() const;

  void setJ(FieldMetric J);
  void setJ(BoutReal value, int x, int y);

  void setBxy(FieldMetric Bxy);

  /// d pitch angle / dx. Needed for vector differentials (Curl)
  const FieldMetric& ShiftTorsion() const;

  ///< Integrated shear (I in BOUT notation)
  const FieldMetric& IntShiftTorsion() const;
  void setIntShiftTorsion(FieldMetric IntShiftTorsion);

  /// Calculate differential geometry quantities from the metric tensor
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
    ASSERT1(transform != nullptr)
    return *transform;
  }

  ///////////////////////////////////////////////////////////
  // Operators
  ///////////////////////////////////////////////////////////

  Field2D DDX(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field2D DDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field2D DDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
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
  Field2D Grad_par(const Field2D& var, CELL_LOC outloc = CELL_DEFAULT,
                   const std::string& method = "DEFAULT");

  Field3D Grad_par(const Field3D& var, CELL_LOC outloc = CELL_DEFAULT,
                   const std::string& method = "DEFAULT");

  /// Advection along magnetic field V*b.Grad(f)
  Field2D Vpar_Grad_par(const Field2D& v, const Field2D& f,
                        CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");

  Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                        CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");

  /// Divergence along magnetic field  Div(b*f) = B.Grad(f/B)
  Field2D Div_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT");

  Field3D Div_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT");

  // Second derivative along magnetic field
  Field2D Grad2_par2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                     const std::string& method = "DEFAULT");

  Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                     const std::string& method = "DEFAULT");
  // Perpendicular Laplacian operator, using only X-Z derivatives
  // NOTE: This might be better bundled with the Laplacian inversion code
  // since it makes use of the same coefficients and FFT routines
  Field2D Delp2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
  Field3D Delp2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
  FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);

  // Full parallel Laplacian operator on scalar field
  // Laplace_par(f) = Div( b (b dot Grad(f)) )
  Field2D Laplace_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT);
  Field3D Laplace_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT);

  // Full Laplacian operator on scalar field
  Field2D Laplace(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& dfdy_boundary_conditions = "free_o3",
                  const std::string& dfdy_dy_region = "");
  Field3D Laplace(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& dfdy_boundary_conditions = "free_o3",
                  const std::string& dfdy_dy_region = "");

  // Full perpendicular Laplacian, in form of inverse of Laplacian operator in LaplaceXY
  // solver
  Field2D Laplace_perpXY(const Field2D& A, const Field2D& f) const;

  /// Christoffel symbol of the second kind (connection coefficients)
  const FieldMetric& G1_11() const;
  const FieldMetric& G1_22() const;
  const FieldMetric& G1_33() const;
  const FieldMetric& G1_12() const;
  const FieldMetric& G1_13() const;
  const FieldMetric& G1_23() const;

  const FieldMetric& G2_11() const;
  const FieldMetric& G2_22() const;
  const FieldMetric& G2_33() const;
  const FieldMetric& G2_12() const;
  const FieldMetric& G2_13() const;
  const FieldMetric& G2_23() const;

  const FieldMetric& G3_11() const;
  const FieldMetric& G3_22() const;
  const FieldMetric& G3_33() const;
  const FieldMetric& G3_12() const;
  const FieldMetric& G3_13() const;
  const FieldMetric& G3_23() const;

  const FieldMetric& G1() const;
  const FieldMetric& G2() const;
  const FieldMetric& G3() const;

  void setG1(const FieldMetric& G1);
  void setG2(const FieldMetric& G2);
  void setG3(const FieldMetric& G3);

  const FieldMetric& Grad2_par2_DDY_invSg(CELL_LOC outloc,
                                          const std::string& method) const;

  const FieldMetric& invSg() const;

  ChristoffelSymbols& christoffel_symbols() const;

  void recalculateAndReset(bool recalculate_staggered,
                           bool force_interpolate_from_centre);

  FieldMetric recalculateJacobian() const;

private:
  int nz; // Size of mesh in Z. This is mesh->ngz-1
  Mesh* localmesh;
  CELL_LOC location;

  DifferentialOperators* differential_operators;

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

  MetricTensor contravariantMetricTensor;
  MetricTensor covariantMetricTensor;

  /// Christoffel symbol of the second kind (connection coefficients)
  mutable std::unique_ptr<ChristoffelSymbols> christoffel_symbols_cache{nullptr};

  void applyToContravariantMetricTensor(
      const std::function<const FieldMetric(const FieldMetric)>& function);

  void applyToCovariantMetricTensor(
      const std::function<const FieldMetric(const FieldMetric)>& function);

  void applyToChristoffelSymbols(
      const std::function<const FieldMetric(const FieldMetric)>& function) const;

  FieldMetric jacobian_cache;

  FieldMetric Bxy_; ///< Magnitude of B = nabla z times nabla x

  FieldMetric G1_, G2_, G3_;

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

  FieldMetric recalculateBxy() const;

  /// Non-uniform meshes. Need to use DDX, DDY
  void correctionForNonUniformMeshes(bool force_interpolate_from_centre);

  void interpolateFieldsFromOtherCoordinates(Options* options,
                                             const Coordinates* coords_in);

  void setBoundaryCells(Options* options, const std::string& suffix);

  FieldMetric getDzFromOptionsFile(Mesh* mesh, const std::string& suffix) const;
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
