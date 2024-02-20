
#include "bout/g_values.hxx"

GValues::GValues(FieldMetric G1, FieldMetric G2, FieldMetric G3)
    : G1_(std::move(G1)), G2_(std::move(G2)), G3_(std::move(G3)){};

GValues::GValues() = default;

const FieldMetric& GValues::G1() const { return G1_; }
const FieldMetric& GValues::G2() const { return G2_; }
const FieldMetric& GValues::G3() const { return G3_; }

void GValues::setG1(const FieldMetric& G1) {}
void GValues::setG2(const FieldMetric& G2) {}
void GValues::setG3(const FieldMetric& G3) {}
