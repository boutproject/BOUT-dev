
#ifndef BOUT_TOKAMAK_COORDINATES_FACTORY_HXX
#define BOUT_TOKAMAK_COORDINATES_FACTORY_HXX

#include "bout.hxx"

class TokamakCoordinatesFactory {

private:
  
  Mesh& mesh_m;
  const Field2D& Rxy_m;
  const Field2D& Bpxy_m;
  const Field2D& Btxy_m;
  const FieldMetric& B_m;
  const FieldMetric& hthe_m;

public:

  TokamakCoordinatesFactory(Mesh& mesh, const Field2D& Rxy, const Field2D& Bpxy, const Field2D& Btxy, const FieldMetric& B, const FieldMetric& hthe)
      : mesh_m(mesh), Rxy_m(Rxy), Bpxy_m(Bpxy), Btxy_m(Btxy), B_m(B), hthe_m(hthe) {}

  Coordinates* make_tokamak_coordinates(const FieldMetric& I,
                                        const BoutReal& sign_of_bp = 1.0) const
  {
    auto* coord = mesh_m.getCoordinates();

    const auto g11 = SQ(Rxy() * Bpxy());
    const auto g22 = 1.0 / SQ(hthe());
    const auto g33 = SQ(I) * g11 + SQ(B()) / g11;
    const auto g12 = 0.0;
    const auto g13 = -I * g11;
    const auto g23 = -sign_of_bp * Btxy() / (hthe() * Bpxy() * Rxy());

    const auto g_11 = 1.0 / g11 + SQ(I * Rxy());
    const auto g_22 = SQ(B() * hthe() / Bpxy());
    const auto g_33 = Rxy() * Rxy();
    const auto g_12 = sign_of_bp * Btxy() * hthe() * I * Rxy() / Bpxy();
    const auto g_13 = I * Rxy() * Rxy();
    const auto g_23 = sign_of_bp * Btxy() * hthe() * Rxy() / Bpxy();

    coord->setMetricTensor(ContravariantMetricTensor(g11, g22, g33, g12, g13, g23),
                           CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));

    coord->setJ(hthe() / Bpxy());
    coord->setBxy(B());

    return coord;
  }

  const Field2D& Rxy() const { return Rxy_m; }
  const Field2D& Bpxy() const { return Bpxy_m; }
  const Field2D& Btxy() const { return Btxy_m; }
  const FieldMetric& B() const { return B_m; }
  const FieldMetric& hthe() const { return hthe_m; }

};

#endif //BOUT_TOKAMAK_COORDINATES_FACTORY_HXX
