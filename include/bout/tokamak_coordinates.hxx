
#ifndef BOUT_TOKAMAK_COORDINATES_HXX
#define BOUT_TOKAMAK_COORDINATES_HXX

#include "bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/christoffel_symbols.hxx"
#include "bout/coordinates.hxx"
#include "bout/field2d.hxx"
#include "bout/utils.hxx"
#include "bout/vector2d.hxx"

struct TokamakOptions {
    Field2D Rxy;
    Field2D Bpxy;
    Field2D Btxy;
    Field2D Bxy;
    Field2D hthe;
    FieldMetric ShearFactor;
    FieldMetric dx;

    TokamakOptions(Mesh& mesh) {
        mesh.get(Rxy, "Rxy");
        //    mesh->get(Zxy, "Zxy");
        mesh.get(Bpxy, "Bpxy");
        mesh.get(Btxy, "Btxy");
        mesh.get(Bxy, "Bxy");
        mesh.get(hthe, "hthe");
        mesh.get(ShearFactor, "sinty");
        mesh.get(dx, "dpsi");
    }

    void normalise(BoutReal Lbar, BoutReal Bbar, BoutReal ShearFactor) {
        Rxy /= Lbar;
        Bpxy /= Bbar;
        Btxy /= Bbar;
        Bxy /= Bbar;
        hthe /= Lbar;
        ShearFactor *= Lbar * Lbar * Bbar * ShearFactor;
        dx /= Lbar * Lbar * Bbar;
    }
};

BoutReal get_sign_of_bp(Field2D Bpxy) {
    if (min(Bpxy, true) < 0.0) {
        return -1.0;
    }
    return 1.0;
}

void set_tokamak_coordinates_on_mesh(TokamakOptions& tokamak_options, Mesh& mesh, const bool noshear, BoutReal Lbar,
                                     BoutReal Bbar, BoutReal ShearFactor = 1.0) {

    tokamak_options.normalise(Lbar, Bbar, ShearFactor);

    const BoutReal sign_of_bp = get_sign_of_bp(tokamak_options.Bpxy);

    if (noshear) {
        ShearFactor = 0.0;
    }

    auto *coord = mesh.getCoordinates();

    const auto g11 = SQ(tokamak_options.Rxy * tokamak_options.Bpxy);
    const auto g22 = 1.0 / SQ(tokamak_options.hthe);
    const auto g33 = SQ(ShearFactor) * g11 + SQ(tokamak_options.Bxy) / g11;
    const auto g12 = 0.0;
    const auto g13 = -ShearFactor * g11;
    const auto g23 = -sign_of_bp * tokamak_options.Btxy /
                     (tokamak_options.hthe * tokamak_options.Bpxy * tokamak_options.Rxy);

    const auto g_11 = 1.0 / g11 + SQ(ShearFactor * tokamak_options.Rxy);
    const auto g_22 = SQ(tokamak_options.Bxy * tokamak_options.hthe / tokamak_options.Bpxy);
    const auto g_33 = tokamak_options.Rxy * tokamak_options.Rxy;
    const auto g_12 =
            sign_of_bp * tokamak_options.Btxy * tokamak_options.hthe * ShearFactor * tokamak_options.Rxy /
            tokamak_options.Bpxy;
    const auto g_13 = ShearFactor * tokamak_options.Rxy * tokamak_options.Rxy;
    const auto g_23 = sign_of_bp * tokamak_options.Btxy * tokamak_options.hthe * tokamak_options.Rxy /
                      tokamak_options.Bpxy;

    coord->setMetricTensor(ContravariantMetricTensor(g11, g22, g33, g12, g13, g23),
                           CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));


    coord->setJ(tokamak_options.hthe / tokamak_options.Bpxy);
    coord->setBxy(tokamak_options.Bxy);
    coord->setDx(tokamak_options.dx);
}

#endif //BOUT_TOKAMAK_COORDINATES_HXX
