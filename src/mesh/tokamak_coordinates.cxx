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
  Field2D Rxy;
  mesh.get(Rxy, "Rxy"); // [m]
  Rxy /= Lbar;

  Field2D Zxy;
  mesh.get(Zxy, "Zxy"); // [m]
  Zxy /= Lbar;

  Field2D Bpxy;
  mesh.get(Bpxy, "Bpxy"); // [T]
  Bpxy /= Bbar;

  Field2D Btxy;
  mesh.get(Btxy, "Btxy"); // [T]
  Btxy /= Bbar;

  Field2D Bxy;
  mesh.get(Bxy, "Bxy"); // [T]
  Bxy /= Bbar;

  Field2D hthe;
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

  const auto g11 = SQ(Rxy * Bpxy);
  const auto g22 = 1.0 / SQ(hthe);
  const auto g33 = SQ(I) * g11 + SQ(Bxy) / g11;
  const auto g12 = 0.0;
  const auto g13 = -I * g11;
  const auto g23 = -sign_of_bp * Btxy / (hthe * Bpxy * Rxy);

  const auto g_11 = 1.0 / g11 + SQ(I * Rxy);
  const auto g_22 = SQ(Bxy * hthe / Bpxy);
  const auto g_33 = Rxy * Rxy;
  const auto g_12 = sign_of_bp * Btxy * hthe * I * Rxy / Bpxy;
  const auto g_13 = I * Rxy * Rxy;
  const auto g_23 = sign_of_bp * Btxy * hthe * Rxy / Bpxy;

  coord->setMetricTensor(ContravariantMetricTensor(g11, g22, g33, g12, g13, g23),
                         CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));

  coord->setJ(hthe / Bpxy);
  coord->setBxy(Bxy);
  coord->setDx(dx);

  return {Rxy, Zxy, Bpxy, Btxy, Bxy, hthe, I, I_unnormalised};
}
} // namespace bout
