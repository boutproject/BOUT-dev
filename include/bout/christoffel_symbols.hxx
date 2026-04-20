#ifndef BOUT_CHRISTOFFELSYMBOLS_HXX
#define BOUT_CHRISTOFFELSYMBOLS_HXX

#include <bout/metric_tensor.hxx>

class Coordinates;

class ChristoffelSymbols {

public:
  explicit ChristoffelSymbols(Coordinates& coordinates);

  const bout::FieldMetric& G1_11() const { return G1_11_m; }
  const bout::FieldMetric& G1_22() const { return G1_22_m; }
  const bout::FieldMetric& G1_33() const { return G1_33_m; }
  const bout::FieldMetric& G1_12() const { return G1_12_m; }
  const bout::FieldMetric& G1_13() const { return G1_13_m; }
  const bout::FieldMetric& G1_23() const { return G1_23_m; }

  const bout::FieldMetric& G2_11() const { return G2_11_m; }
  const bout::FieldMetric& G2_22() const { return G2_22_m; }
  const bout::FieldMetric& G2_33() const { return G2_33_m; }
  const bout::FieldMetric& G2_12() const { return G2_12_m; }
  const bout::FieldMetric& G2_13() const { return G2_13_m; }
  const bout::FieldMetric& G2_23() const { return G2_23_m; }

  const bout::FieldMetric& G3_11() const { return G3_11_m; }
  const bout::FieldMetric& G3_22() const { return G3_22_m; }
  const bout::FieldMetric& G3_33() const { return G3_33_m; }
  const bout::FieldMetric& G3_12() const { return G3_12_m; }
  const bout::FieldMetric& G3_13() const { return G3_13_m; }
  const bout::FieldMetric& G3_23() const { return G3_23_m; }

  // Transforms the ChristoffelSymbols by applying the given function to every element
  template <class F>
  void map(F function) {
    G1_11_m = function(G1_11_m);
    G1_22_m = function(G1_22_m);
    G1_33_m = function(G1_33_m);
    G1_12_m = function(G1_12_m);
    G1_13_m = function(G1_13_m);
    G1_23_m = function(G1_23_m);
    G2_11_m = function(G2_11_m);
    G2_22_m = function(G2_22_m);
    G2_33_m = function(G2_33_m);
    G2_12_m = function(G2_12_m);
    G2_13_m = function(G2_13_m);
    G2_23_m = function(G2_23_m);
    G3_11_m = function(G3_11_m);
    G3_22_m = function(G3_22_m);
    G3_33_m = function(G3_33_m);
    G3_12_m = function(G3_12_m);
    G3_13_m = function(G3_13_m);
    G3_23_m = function(G3_23_m);
  }

private:
  bout::FieldMetric G1_11_m, G1_22_m, G1_33_m, G1_12_m, G1_13_m, G1_23_m;
  bout::FieldMetric G2_11_m, G2_22_m, G2_33_m, G2_12_m, G2_13_m, G2_23_m;
  bout::FieldMetric G3_11_m, G3_22_m, G3_33_m, G3_12_m, G3_13_m, G3_23_m;
};

#endif //BOUT_CHRISTOFFELSYMBOLS_HXX
