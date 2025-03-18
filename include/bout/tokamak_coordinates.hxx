
#ifndef BOUT_TOKAMAK_COORDINATES_HXX
#define BOUT_TOKAMAK_COORDINATES_HXX

#include "bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/metric_tensor.hxx"

using FieldMetric = MetricTensor::FieldMetric;

struct TokamakOptions {
    Field2D Rxy;
    Field2D Bpxy;
    Field2D Btxy;
    Field2D Bxy;
    Field2D hthe;
    FieldMetric I;
    FieldMetric dx;

    TokamakOptions(Mesh &mesh);

    void normalise(BoutReal Lbar, BoutReal Bbar, BoutReal ShearFactor) {
        Rxy /= Lbar;
        Bpxy /= Bbar;
        Btxy /= Bbar;
        Bxy /= Bbar;
        hthe /= Lbar;
        I *= Lbar * Lbar * Bbar * ShearFactor;
        dx /= Lbar * Lbar * Bbar;
    }
};

BoutReal get_sign_of_bp(const Field2D& Bpxy);

void set_tokamak_coordinates_on_mesh(TokamakOptions &tokamak_options, Mesh &mesh, BoutReal Lbar,
                                     BoutReal Bbar, BoutReal ShearFactor = 0.0);

#endif //BOUT_TOKAMAK_COORDINATES_HXX
