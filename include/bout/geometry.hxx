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

#include "differential_operators.hxx"
#include "metricTensor.hxx"

using FieldMetric = MetricTensor::FieldMetric;

/*!
 * Represents the geometry a coordinate system, and associated operators
 *
 * This is a container for a collection of metric tensor components
 */
class Geometry {

public:
  Geometry(const FieldMetric& J, const FieldMetric& Bxy, const FieldMetric& g11,
           const FieldMetric& g22, const FieldMetric& g33, const FieldMetric& g12,
           const FieldMetric& g13, const FieldMetric& g23, const FieldMetric& g_11,
           const FieldMetric& g_22, const FieldMetric& g_33, const FieldMetric& g_12,
           const FieldMetric& g_13, const FieldMetric& g_23,
           DifferentialOperators* differential_operators);

  Geometry(Mesh* mesh, DifferentialOperators* differential_operators);

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

  ///< Coordinate system Jacobian, so volume of cell is J*dx*dy*dz
  const FieldMetric& J() const;

  ///< Magnitude of B = nabla z times nabla x
  const FieldMetric& Bxy() const;

  void setContravariantMetricTensor(MetricTensor metric_tensor, CELL_LOC cell_location,
                                    const std::string& region = "RGN_ALL");

  void setCovariantMetricTensor(MetricTensor metric_tensor, CELL_LOC cell_location,
                                const std::string& region = "RGN_ALL");

  void setG1_11(FieldMetric G1_11);
  void setG1_22(FieldMetric G1_22);
  void setG1_33(FieldMetric G1_33);
  void setG1_12(FieldMetric G1_12);
  void setG1_13(FieldMetric G1_13);
  void setG1_23(FieldMetric G1_23);

  void setG2_11(FieldMetric G2_11);
  void setG2_22(FieldMetric G2_22);
  void setG2_33(FieldMetric G2_33);
  void setG2_12(FieldMetric G2_12);
  void setG2_13(FieldMetric G2_13);
  void setG2_23(FieldMetric G2_23);

  void setG3_11(FieldMetric G3_11);
  void setG3_22(FieldMetric G3_22);
  void setG3_33(FieldMetric G3_33);
  void setG3_12(FieldMetric G3_12);
  void setG3_13(FieldMetric G3_13);
  void setG3_23(FieldMetric G3_23);

  void setG3(FieldMetric G3);
  void setG1(FieldMetric G1);
  void setG2(FieldMetric G2);

  void setJ(FieldMetric J);
  void setJ(BoutReal value, int x, int y);

  void setBxy(FieldMetric Bxy);

  void calcCovariant(CELL_LOC cell_location, const std::string& region = "RGN_ALL");

  /// Invert covariant metric to get contravariant components
  void calcContravariant(CELL_LOC cell_location, const std::string& region = "RGN_ALL");

  /// Calculate Christoffel symbol terms
  void CalculateChristoffelSymbols(FieldMetric& dx, FieldMetric& dy);

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
  /// Christoffel symbol of the second kind (connection coefficients)
  FieldMetric G1_11_, G1_22_, G1_33_, G1_12_, G1_13_, G1_23_;
  FieldMetric G2_11_, G2_22_, G2_33_, G2_12_, G2_13_, G2_23_;
  FieldMetric G3_11_, G3_22_, G3_33_, G3_12_, G3_13_, G3_23_;

  FieldMetric G1_, G2_, G3_;

  MetricTensor contravariantMetricTensor;
  MetricTensor covariantMetricTensor;

  FieldMetric this_J;
  FieldMetric this_Bxy; ///< Magnitude of B = nabla z times nabla x

  DifferentialOperators* differential_operators;
};

#endif // __GEOMETRY_H__
