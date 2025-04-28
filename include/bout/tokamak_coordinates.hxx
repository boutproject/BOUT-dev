
#ifndef BOUT_TOKAMAK_COORDINATES_HXX
#define BOUT_TOKAMAK_COORDINATES_HXX

#include "bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/metric_tensor.hxx"

using FieldMetric = MetricTensor::FieldMetric;

namespace bout {

    struct Coordinates3D {

        Field3D x;
        Field3D y;
        Field3D z;

        Coordinates3D(Field3D x, Field3D y, Field3D z) : x(x), y(y), z(z) {}
    };

    struct TokamakOptions {
        Field2D Rxy;
        Field2D Zxy;
        Field2D Bpxy;
        Field2D Btxy;
        Field2D Bxy;
        Field2D hthe;
        FieldMetric I;
        FieldMetric dx;

        TokamakOptions(Mesh &mesh);

        void normalise(BoutReal Lbar, BoutReal Bbar, BoutReal ShearFactor);

        Coordinates3D CylindricalCoordinatesToCartesian();
    };

    BoutReal get_sign_of_bp(const Field2D &Bpxy);

    void set_tokamak_coordinates_on_mesh(TokamakOptions &tokamak_options, Mesh &mesh, BoutReal Lbar,
                                         BoutReal Bbar, BoutReal ShearFactor = 0.0);

}

#endif //BOUT_TOKAMAK_COORDINATES_HXX
