/**************************************************************************
 * Describes coordinate geometry
 *
 * ChangeLog
 * =========
 * 
 * 2021-11-23 Tom Chapman <tpc522@york.ac.uk>
 *    * Extracted from Coordinates
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

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include "coordinates.hxx"
#include "differential_operators.hxx"
#include "metricTensor.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/paralleltransform.hxx"
#include "bout/utils.hxx"
#include <bout/bout_types.hxx>

/*!
 * Represents the geometry a coordinate system, and associated operators
 *
 * This is a container for a collection of metric tensor components
 */
class Geometry {

public:
#if BOUT_USE_METRIC_3D
  using FieldMetric = Field3D;
#else
  using FieldMetric = Field2D;
#endif

  //  Geometry();

  Geometry(
      //      Coordinates& coordinates,
      FieldMetric J, FieldMetric Bxy, FieldMetric g11, FieldMetric g22, FieldMetric g33,
      FieldMetric g12, FieldMetric g13, FieldMetric g23, FieldMetric g_11,
      FieldMetric g_22, FieldMetric g_33, FieldMetric g_12, FieldMetric g_13,
      FieldMetric g_23, DifferentialOperators differential_operators);

  Geometry(Mesh* mesh, DifferentialOperators differential_operators);

  /// Christoffel symbol of the second kind (connection coefficients)
  FieldMetric G1_11, G1_22, G1_33, G1_12, G1_13, G1_23;
  FieldMetric G2_11, G2_22, G2_33, G2_12, G2_13, G2_23;
  FieldMetric G3_11, G3_22, G3_33, G3_12, G3_13, G3_23;

  FieldMetric G1, G2, G3;

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
  //  MetricTensor& getCovariantMetricTensor();

  ///< Coordinate system Jacobian, so volume of cell is J*dx*dy*dz
  const FieldMetric& J() const;

  ///< Magnitude of B = nabla z times nabla x
  const FieldMetric& Bxy() const;

  void setContravariantMetricTensor(MetricTensor metric_tensor, CELL_LOC cell_location,
                                    const std::string& region = "RGN_ALL");

  void setCovariantMetricTensor(MetricTensor metric_tensor, CELL_LOC cell_location,
                                const std::string& region = "RGN_ALL");

  void setJ(FieldMetric J);
  void setJ(BoutReal value, int x, int y);

  void setBxy(FieldMetric Bxy);

  void calcCovariant(CELL_LOC cell_location, const std::string& region = "RGN_ALL");

  /// Invert covariant metric to get contravariant components
  void calcContravariant(CELL_LOC cell_location, const std::string& region = "RGN_ALL");

  //  void jacobian(bool extrapolate_x, bool extrapolate_y); ///< Calculate J and Bxy

  void CalculateChristoffelSymbols(); /// Calculate Christoffel symbol terms

  // check that covariant tensors are positive (if expected) and finite (always)
  void checkCovariant(int ystart);

  // check that contravariant tensors are positive (if expected) and finite (always)
  void checkContravariant(int ystart);

  FieldMetric recalculateJacobian();
  FieldMetric recalculateBxy();

  void applyToContravariantMetricTensor(
      std::function<const FieldMetric(const FieldMetric)> function);

  void applyToCovariantMetricTensor(
      std::function<const FieldMetric(const FieldMetric)> function);

private:
  //  Coordinates& coordinates;

  //  int calculateGeometry(FieldMetric& dx, FieldMetric& dy, FieldMetric& dz,
  //                        bool recalculate_staggered, bool force_interpolate_from_centre);

  MetricTensor contravariantMetricTensor;
  MetricTensor covariantMetricTensor;

  FieldMetric this_J;
  FieldMetric this_Bxy; ///< Magnitude of B = nabla z times nabla x

  DifferentialOperators differential_operators;

  //  template <typename T, typename... Ts>
  //  void communicate(T& t, Ts... ts);
};

#endif // __GEOMETRY_H__
