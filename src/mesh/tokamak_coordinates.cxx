#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
#include <bout/field2d.hxx>
#include <bout/mesh.hxx>
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

  Coordinates::FieldMetric I;
  if (no_shear) {
    I = 0.0;
  } else {
    mesh.get(I, "sinty"); // [m^-2 T^-1]
  }
  const auto I_unnormalised = I;
  I *= Lbar * Lbar * Bbar * shear_factor;

  Coordinates::FieldMetric dx;
  if (mesh.get(dx, "dpsi") != 0) {
    dx = mesh.getCoordinates()->dx;
  }
  dx /= Lbar * Lbar * Bbar;

  const BoutReal sign_of_bp = min(Bpxy, true) < 0.0 ? -1.0 : 1.0;

  auto* coords = mesh.getCoordinates();

  coords->Bxy = Bxy;
  coords->dx = dx;

  coords->g11 = SQ(Rxy * Bpxy);
  coords->g22 = 1.0 / SQ(hthe);
  coords->g33 = SQ(I) * coords->g11 + SQ(Bxy) / coords->g11;
  coords->g12 = 0.0;
  coords->g13 = -I * coords->g11;
  coords->g23 = -sign_of_bp * Btxy / (hthe * Bpxy * Rxy);

  coords->J = hthe / Bpxy;

  coords->g_11 = 1.0 / coords->g11 + SQ(I * Rxy);
  coords->g_22 = SQ(Bxy * hthe / Bpxy);
  coords->g_33 = Rxy * Rxy;
  coords->g_12 = sign_of_bp * Btxy * hthe * I * Rxy / Bpxy;
  coords->g_13 = I * Rxy * Rxy;
  coords->g_23 = sign_of_bp * Btxy * hthe * Rxy / Bpxy;

  coords->geometry();

  return {Rxy, Zxy, Bpxy, Btxy, Bxy, hthe, I, I_unnormalised};
}
} // namespace bout
