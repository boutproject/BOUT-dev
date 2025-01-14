
#ifndef BOUT_TOKAMAK_COORDINATES_HXX
#define BOUT_TOKAMAK_COORDINATES_HXX

#include "bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/christoffel_symbols.hxx"
#include "bout/coordinates.hxx"
#include "bout/field2d.hxx"
#include "bout/utils.hxx"
#include "bout/vector2d.hxx"

struct TokamakCoordinates {
    Field2D Rxy;
    Field2D Bpxy;
    Field2D Btxy;
    Field2D Bxy;
    Field2D hthe;
    FieldMetric ShearFactor;
    FieldMetric dx;

    TokamakCoordinates(Mesh &mesh) {
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

void set_tokamak_coordinates_on_mesh(TokamakCoordinates& tokamakCoordinates, Mesh& mesh, const bool noshear, 
                                BoutReal Lbar, BoutReal Bbar, BoutReal ShearFactor = 1.0) {

    tokamakCoordinates.normalise(Lbar, Bbar, ShearFactor);

    const BoutReal sign_of_bp = get_sign_of_bp(tokamakCoordinates.Bpxy);

    if (noshear) {
        ShearFactor = 0.0;
    }

    auto *coord = mesh.getCoordinates();

    const auto g11 = SQ(tokamakCoordinates.Rxy * tokamakCoordinates.Bpxy);
    const auto g22 = 1.0 / SQ(tokamakCoordinates.hthe);
    const auto g33 = SQ(ShearFactor) * g11 + SQ(tokamakCoordinates.Bxy) / g11;
    const auto g12 = 0.0;
    const auto g13 = -ShearFactor * g11;
    const auto g23 = -sign_of_bp * tokamakCoordinates.Btxy /
                     (tokamakCoordinates.hthe * tokamakCoordinates.Bpxy * tokamakCoordinates.Rxy);

    const auto g_11 = 1.0 / g11 + SQ(ShearFactor * tokamakCoordinates.Rxy);
    const auto g_22 = SQ(tokamakCoordinates.Bxy * tokamakCoordinates.hthe / tokamakCoordinates.Bpxy);
    const auto g_33 = tokamakCoordinates.Rxy * tokamakCoordinates.Rxy;
    const auto g_12 =
            sign_of_bp * tokamakCoordinates.Btxy * tokamakCoordinates.hthe * ShearFactor * tokamakCoordinates.Rxy /
            tokamakCoordinates.Bpxy;
    const auto g_13 = ShearFactor * tokamakCoordinates.Rxy * tokamakCoordinates.Rxy;
    const auto g_23 = sign_of_bp * tokamakCoordinates.Btxy * tokamakCoordinates.hthe * tokamakCoordinates.Rxy /
                      tokamakCoordinates.Bpxy;

    coord->setMetricTensor(ContravariantMetricTensor(g11, g22, g33, g12, g13, g23),
                           CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));


    coord->setJ(tokamakCoordinates.hthe / tokamakCoordinates.Bpxy);
    coord->setBxy(tokamakCoordinates.Bxy);
    coord->setDx(tokamakCoordinates.dx);
}

#endif //BOUT_TOKAMAK_COORDINATES_HXX
