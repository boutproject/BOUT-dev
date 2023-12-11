
#ifndef BOUT_CHRISTOFFELSYMBOLS_HXX
#define BOUT_CHRISTOFFELSYMBOLS_HXX

#include "differential_operators.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/metricTensor.hxx"
#include <bout/bout_types.hxx>

using FieldMetric = MetricTensor::FieldMetric;

class ChristoffelSymbols {

public:
  ChristoffelSymbols(FieldMetric G1_11, FieldMetric G1_22, FieldMetric G1_33,
                     FieldMetric G1_12, FieldMetric G1_13, FieldMetric G1_23,
                     FieldMetric G2_11, FieldMetric G2_22, FieldMetric G2_33,
                     FieldMetric G2_12, FieldMetric G2_13, FieldMetric G2_23,
                     FieldMetric G3_11, FieldMetric G3_22, FieldMetric G3_33,
                     FieldMetric G3_12, FieldMetric G3_13, FieldMetric G3_23,
                     FieldMetric G1, FieldMetric G2, FieldMetric G3,
                     DifferentialOperators* differential_operators);

  //  ChristoffelSymbols(BoutReal g11, BoutReal g22, BoutReal g33, BoutReal g12, BoutReal g13,
  //                     BoutReal g23, Mesh* mesh);

  explicit ChristoffelSymbols(DifferentialOperators* differential_operators);

  ChristoffelSymbols(Mesh* mesh, DifferentialOperators* differential_operators);

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

  void setG1(FieldMetric& G1);
  void setG2(FieldMetric& G2);
  void setG3(FieldMetric& G3);

  void setChristoffelSymbols(
      const FieldMetric& G1_11, const FieldMetric& G1_22, const FieldMetric& G1_33,
      const FieldMetric& G1_12, const FieldMetric& G1_13, const FieldMetric& G1_23,
      const FieldMetric& G2_11, const FieldMetric& G2_22, const FieldMetric& G2_33,
      const FieldMetric& G2_12, const FieldMetric& G2_13, const FieldMetric& G2_23,
      const FieldMetric& G3_11, const FieldMetric& G3_22, const FieldMetric& G3_33,
      const FieldMetric& G3_12, const FieldMetric& G3_13, const FieldMetric& G3_23,
      const FieldMetric& G1, const FieldMetric& G2, const FieldMetric& G3);

  //  void Allocate();
  //
  //  void setLocation(CELL_LOC location);
  //
  //  // Transforms the ChristoffelSymbols by applying the given function to every component
  //  void map(std::function<const Field2D(const FieldMetric)> function);
  //
  //  ChristoffelSymbols
  //  applyToComponents(std::function<const FieldMetric(const FieldMetric)> function) const;

  void CalculateChristoffelSymbols(const FieldMetric& dx, const FieldMetric& dy,
                                   const MetricTensor& contravariantMetricTensor,
                                   const MetricTensor& covariantMetricTensor);

  // Transforms the ChristoffelSymbols by applying the given function to every element
  void map(const std::function<const FieldMetric(const FieldMetric)>& function);

private:
  FieldMetric G1_11_, G1_22_, G1_33_, G1_12_, G1_13_, G1_23_;
  FieldMetric G2_11_, G2_22_, G2_33_, G2_12_, G2_13_, G2_23_;
  FieldMetric G3_11_, G3_22_, G3_33_, G3_12_, G3_13_, G3_23_;
  FieldMetric G1_, G2_, G3_;

  DifferentialOperators* differential_operators;

  ChristoffelSymbols applyToComponents(
      const std::function<const FieldMetric(const FieldMetric)>& function) const;
};

#endif //BOUT_CHRISTOFFELSYMBOLS_HXX
