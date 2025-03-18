
#include "bout/tokamak_coordinates.hxx"
#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/utils.hxx"

BoutReal get_sign_of_bp(const Field2D& Bpxy) {
    if (min(Bpxy, true) < 0.0) {
        return -1.0;
    }
    return 1.0;
}

TokamakOptions::TokamakOptions(Mesh &mesh) {
    mesh.get(Rxy, "Rxy");
    //    mesh->get(Zxy, "Zxy");
    mesh.get(Bpxy, "Bpxy");
    mesh.get(Btxy, "Btxy");
    mesh.get(Bxy, "Bxy");
    mesh.get(hthe, "hthe");
    mesh.get(I, "sinty");
    if (mesh.get(dx, "dpsi")) {
        dx = mesh.getCoordinates()->dx();
    }
}

void set_tokamak_coordinates_on_mesh(TokamakOptions &tokamak_options, Mesh &mesh, BoutReal Lbar,
                                     BoutReal
                                     Bbar, BoutReal
                                     ShearFactor) {

    tokamak_options.normalise(Lbar, Bbar, ShearFactor);

    const BoutReal sign_of_bp = get_sign_of_bp(tokamak_options.Bpxy);

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