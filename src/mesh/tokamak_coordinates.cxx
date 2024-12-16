#include <bout/tokamak_coordinates.hxx>


void TokamakCoordinates::normalise(BoutReal Lbar, BoutReal Bbar, BoutReal ShearFactor) {

    Rxy_m /= Lbar;
    Bpxy_m /= Bbar;
    Btxy_m /= Bbar;
    Bxy_m /= Bbar;
    hthe_m /= Lbar;
    ShearFactor_m *= Lbar * Lbar * Bbar * ShearFactor;
    dx_m /= Lbar * Lbar * Bbar;
  }

  TokamakCoordinates::TokamakCoordinates(Mesh& mesh) : mesh_m(mesh) {

    mesh.get(Rxy_m, "Rxy");
    //    mesh->get(Zxy, "Zxy");
    mesh.get(Bpxy_m, "Bpxy");
    mesh.get(Btxy_m, "Btxy");
    mesh.get(Bxy_m, "Bxy");
    mesh.get(hthe_m, "hthe");
    mesh.get(ShearFactor_m, "sinty");
    mesh.get(dx_m, "dpsi");
  }

  BoutReal TokamakCoordinates::get_sign_of_bp() {
    if (min(Bpxy_m, true) < 0.0) {
      return -1.0;
    }
    return 1.0;
  }

  Coordinates* TokamakCoordinates::make_coordinates(const bool noshear, BoutReal Lbar,
                                                    BoutReal Bbar, BoutReal ShearFactor) {

    normalise(Lbar, Bbar, ShearFactor);

    const BoutReal sign_of_bp = get_sign_of_bp();

    if (noshear) {
      ShearFactor_m = 0.0;
    }

    auto* coord = mesh_m.getCoordinates();

    const auto g11 = SQ(Rxy_m * Bpxy_m);
    const auto g22 = 1.0 / SQ(hthe_m);
    const auto g33 = SQ(ShearFactor_m) * g11 + SQ(Bxy_m) / g11;
    const auto g12 = 0.0;
    const auto g13 = -ShearFactor_m * g11;
    const auto g23 = -sign_of_bp * Btxy_m / (hthe_m * Bpxy_m * Rxy_m);

    const auto g_11 = 1.0 / g11 + SQ(ShearFactor_m * Rxy_m);
    const auto g_22 = SQ(Bxy_m * hthe_m / Bpxy_m);
    const auto g_33 = Rxy_m * Rxy_m;
    const auto g_12 = sign_of_bp * Btxy_m * hthe_m * ShearFactor_m * Rxy_m / Bpxy_m;
    const auto g_13 = ShearFactor_m * Rxy_m * Rxy_m;
    const auto g_23 = sign_of_bp * Btxy_m * hthe_m * Rxy_m / Bpxy_m;

    coord->setMetricTensor(ContravariantMetricTensor(g11, g22, g33, g12, g13, g23),
                           CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));

    coord->setJ(hthe_m / Bpxy_m);
    coord->setBxy(Bxy_m);
    coord->setDx(dx_m);

    return coord;
  }
