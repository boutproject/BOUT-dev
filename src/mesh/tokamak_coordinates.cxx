#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
#include <bout/field2d.hxx>
#include <bout/mesh.hxx>
#include <bout/metric_tensor.hxx>
#include <bout/tokamak_coordinates.hxx>
#include <bout/utils.hxx>

namespace bout {
TokamakCoordinates set_tokamak_coordinates(Mesh& mesh, BoutReal Lbar, BoutReal Bbar,
                                           bool no_shear, BoutReal shear_factor) {
  FieldMetric Rxy;
  mesh.get(Rxy, "Rxy"); // [m]
  Rxy /= Lbar;

  FieldMetric Zxy;
  mesh.get(Zxy, "Zxy"); // [m]
  Zxy /= Lbar;

  FieldMetric Bpxy;
  mesh.get(Bpxy, "Bpxy"); // [T]
  Bpxy /= Bbar;

  FieldMetric Btxy;
  mesh.get(Btxy, "Btxy"); // [T]
  Btxy /= Bbar;

  FieldMetric Bxy;
  mesh.get(Bxy, "Bxy"); // [T]
  Bxy /= Bbar;

  FieldMetric hthe;
  mesh.get(hthe, "hthe"); // [m / radian]
  hthe /= Lbar;

  bout::FieldMetric I;
  if (no_shear) {
    I = 0.0;
  } else {
    mesh.get(I, "sinty"); // [m^-2 T^-1]
  }
  const auto I_unnormalised = I;
  I *= Lbar * Lbar * Bbar * shear_factor;

  bout::FieldMetric dx;
  if (mesh.get(dx, "dpsi") != 0) {
    dx = mesh.getCoordinates()->dx();
  }
  dx /= Lbar * Lbar * Bbar;

  const BoutReal sign_of_bp = min(Bpxy, true) < 0.0 ? -1.0 : 1.0;

  auto* coord = mesh.getCoordinates();

  const FieldMetric g11 = SQ(Rxy * Bpxy);
  const FieldMetric g22 = 1.0 / SQ(hthe);
  const FieldMetric g33 = SQ(I) * g11 + SQ(Bxy) / g11;
  const FieldMetric g12 = 0.0;
  const FieldMetric g13 = -I * g11;
  const FieldMetric g23 = -sign_of_bp * Btxy / (hthe * Bpxy * Rxy);

  const FieldMetric g_11 = 1.0 / g11 + SQ(I * Rxy);
  const FieldMetric g_22 = SQ(Bxy * hthe / Bpxy);
  const FieldMetric g_33 = Rxy * Rxy;
  const FieldMetric g_12 = sign_of_bp * Btxy * hthe * I * Rxy / Bpxy;
  const FieldMetric g_13 = I * Rxy * Rxy;
  const FieldMetric g_23 = sign_of_bp * Btxy * hthe * Rxy / Bpxy;

  coord->setMetricTensor(ContravariantMetricTensor(g11, g22, g33, g12, g13, g23),
                         CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));

  coord->setJ(FieldMetric{hthe / Bpxy});
  coord->setBxy(Bxy);
  coord->setDx(dx);

  return {Rxy, Zxy, Bpxy, Btxy, Bxy, hthe, I, I_unnormalised};
}

MetricNormaliser TokamakOrFCIMetricNormaliser(const Mesh* mesh, BoutReal Bnorm,
                                              BoutReal rho_s0) {
  if (mesh->isFci()) {
    return {.g{SQ(rho_s0)}, .J{rho_s0 * rho_s0 * rho_s0}, .Bxy{Bnorm}};
  }
  return {.g11{SQ(Bnorm * rho_s0)},
          .g11_mul{true},
          .g22{SQ(rho_s0)},
          .g33{SQ(rho_s0)},
          .g12{Bnorm},
          .g12_mul{true},
          .g13{Bnorm},
          .g13_mul{true},
          .g23{SQ(rho_s0)},
          .dx{rho_s0 * rho_s0 * Bnorm},
          .J{Bnorm / rho_s0},
          .J_mul{true},
          .Bxy{Bnorm}};
}

} // namespace bout
