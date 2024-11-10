
#ifndef BOUT_TOKAMAK_COORDINATES_FACTORY_HXX
#define BOUT_TOKAMAK_COORDINATES_FACTORY_HXX

#include "bout.hxx"

class TokamakCoordinatesFactory {

private:
  
  Mesh& mesh_m;
  Field2D Rxy_m;
  Field2D Bpxy_m;
  Field2D Btxy_m;
  Field2D Bxy_m;
  Field2D hthe_m;
  FieldMetric ShearFactor_m;
  FieldMetric dx_m;
  Vector2D b0xcv_m; // Curvature term


public:

  TokamakCoordinatesFactory(Mesh& mesh)
      : mesh_m(mesh) {

    mesh.get(Rxy_m, "Rxy");
    //    mesh->get(Zxy, "Zxy");
    mesh.get(Bpxy_m, "Bpxy");
    mesh.get(Btxy_m, "Btxy");
    mesh.get(Bxy_m, "Bxy");
    mesh.get(hthe_m, "hthe");
    mesh.get(ShearFactor_m, "sinty");
    mesh.get(dx_m, "dpsi");
  }

  //  void setShearFactor(const bool noshear = false) {
  //
  //    if (shifted_metric_method) {
  //      // No integrated shear in metric
  //      ShearFactor_m = 0.0;
  //    }
  //  }

  BoutReal get_sign_of_bp() {
    if (min(Bpxy_m, true) < 0.0) {
      return -1.0;
    }
    return 1.0;
  }

  Coordinates* make_tokamak_coordinates(const bool noshear, const bool include_curvature) {

    BoutReal sign_of_bp = get_sign_of_bp();

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

    if (!include_curvature) {
      b0xcv_m = 0.0;
    }

    if (noshear) {
      if (include_curvature) {
        b0xcv_m.z += ShearFactor_m * b0xcv_m.x;
      }
      ShearFactor_m = 0.0;
    }
    //      if (ShiftXderivs and not!mesh->IncIntShear) {
    //        // Dimits style, using local coordinate system
    //        if (include_curvature) {
    //          b0xcv.z += I * b0xcv.x;
    //        }
    //        I = 0.0; // I disappears from metric
    //      }

    return coord;
  }

  void normalise(BoutReal Lbar, BoutReal Bbar) {

    Rxy_m /= Lbar;
    Bpxy_m /= Bbar;
    Btxy_m / Bbar;
    Bxy_m / Bbar;
    hthe_m / Lbar;
    ShearFactor_m *= Lbar * Lbar * Bbar;
    dx_m /= Lbar * Lbar * Bbar;

    // Normalise curvature term
    b0xcv_m.x /= Bbar;
    b0xcv_m.y *= Lbar * Lbar;
    b0xcv_m.z *= Lbar * Lbar;
  }

  const Field2D& get_Rxy() const { return Rxy_m; }
  const Field2D& get_Bpxy() const { return Bpxy_m; }
  const Field2D& get_Btxy() const { return Btxy_m; }
  const Field2D& get_Bxy() const { return Bxy_m; }
  const Field2D& get_hthe() const { return hthe_m; }
  const FieldMetric& get_ShearFactor() const { return ShearFactor_m; }

  void set_ShearFactor(FieldMetric& shearFactor) { ShearFactor_m = shearFactor; }
};

#endif //BOUT_TOKAMAK_COORDINATES_FACTORY_HXX
