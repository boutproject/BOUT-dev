
#ifndef BOUT_CHRISTOFFELSYMBOLS_HXX
#define BOUT_CHRISTOFFELSYMBOLS_HXX

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/metricTensor.hxx"
#include <bout/bout_types.hxx>

using FieldMetric = MetricTensor::FieldMetric;

class Coordinates;

class ChristoffelSymbols {

public:
  ChristoffelSymbols(FieldMetric G1_11, FieldMetric G1_22, FieldMetric G1_33,
                     FieldMetric G1_12, FieldMetric G1_13, FieldMetric G1_23,
                     FieldMetric G2_11, FieldMetric G2_22, FieldMetric G2_33,
                     FieldMetric G2_12, FieldMetric G2_13, FieldMetric G2_23,
                     FieldMetric G3_11, FieldMetric G3_22, FieldMetric G3_33,
                     FieldMetric G3_12, FieldMetric G3_13, FieldMetric G3_23);

  explicit ChristoffelSymbols(const Coordinates& coordinates);

  //  ChristoffelSymbols(BoutReal g11, BoutReal g22, BoutReal g33, BoutReal g12, BoutReal g13,
  //                     BoutReal g23, Mesh* mesh);

  const FieldMetric& G1_11() const { return G1_11_; }
  const FieldMetric& G1_22() const { return G1_22_; }
  const FieldMetric& G1_33() const { return G1_33_; }
  const FieldMetric& G1_12() const { return G1_12_; }
  const FieldMetric& G1_13() const { return G1_13_; }
  const FieldMetric& G1_23() const { return G1_23_; }

  const FieldMetric& G2_11() const { return G2_11_; }
  const FieldMetric& G2_22() const { return G2_22_; }
  const FieldMetric& G2_33() const { return G2_33_; }
  const FieldMetric& G2_12() const { return G2_12_; }
  const FieldMetric& G2_13() const { return G2_13_; }
  const FieldMetric& G2_23() const { return G2_23_; }

  const FieldMetric& G3_11() const { return G3_11_; }
  const FieldMetric& G3_22() const { return G3_22_; }
  const FieldMetric& G3_33() const { return G3_33_; }
  const FieldMetric& G3_12() const { return G3_12_; }
  const FieldMetric& G3_13() const { return G3_13_; }
  const FieldMetric& G3_23() const { return G3_23_; }

  void setChristoffelSymbols(const FieldMetric& G1_11, const FieldMetric& G1_22,
                             const FieldMetric& G1_33, const FieldMetric& G1_12,
                             const FieldMetric& G1_13, const FieldMetric& G1_23,
                             const FieldMetric& G2_11, const FieldMetric& G2_22,
                             const FieldMetric& G2_33, const FieldMetric& G2_12,
                             const FieldMetric& G2_13, const FieldMetric& G2_23,
                             const FieldMetric& G3_11, const FieldMetric& G3_22,
                             const FieldMetric& G3_33, const FieldMetric& G3_12,
                             const FieldMetric& G3_13, const FieldMetric& G3_23) {
    G1_11_ = G1_11;
    G1_22_ = G1_22;
    G1_33_ = G1_33;
    G1_12_ = G1_12;
    G1_13_ = G1_13;
    G1_23_ = G1_23;
    G2_11_ = G2_11;
    G2_22_ = G2_22;
    G2_33_ = G2_33;
    G2_12_ = G2_12;
    G2_13_ = G2_13;
    G2_23_ = G2_23;
    G3_11_ = G3_11;
    G3_22_ = G3_22;
    G3_33_ = G3_33;
    G3_12_ = G3_12;
    G3_13_ = G3_13;
    G3_23_ = G3_23;
  }

  //  void Allocate();
  //
  //  void setLocation(CELL_LOC location);

  // Transforms the ChristoffelSymbols by applying the given function to every element
  void map(const std::function<const FieldMetric(const FieldMetric)>& function);

private:
  FieldMetric G1_11_, G1_22_, G1_33_, G1_12_, G1_13_, G1_23_;
  FieldMetric G2_11_, G2_22_, G2_33_, G2_12_, G2_13_, G2_23_;
  FieldMetric G3_11_, G3_22_, G3_33_, G3_12_, G3_13_, G3_23_;

  ChristoffelSymbols applyToComponents(
      const std::function<const FieldMetric(const FieldMetric)>& function) const;
};

#endif //BOUT_CHRISTOFFELSYMBOLS_HXX
