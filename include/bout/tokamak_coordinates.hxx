
#ifndef BOUT_TOKAMAK_COORDINATES_HXX
#define BOUT_TOKAMAK_COORDINATES_HXX

#include "bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/christoffel_symbols.hxx"
#include "bout/coordinates.hxx"
#include "bout/field2d.hxx"
#include "bout/utils.hxx"
#include "bout/vector2d.hxx"


class TokamakCoordinates {

private:

  Mesh& mesh_m;
  Field2D Rxy_m;
  Field2D Bpxy_m;
  Field2D Btxy_m;
  Field2D Bxy_m;
  Field2D hthe_m;
  FieldMetric ShearFactor_m;
  FieldMetric dx_m;

  void normalise(BoutReal Lbar, BoutReal Bbar, BoutReal ShearFactor);

  BoutReal get_sign_of_bp();


public:

  TokamakCoordinates(Mesh& mesh);

  Coordinates* make_coordinates(const bool noshear, BoutReal Lbar, BoutReal Bbar,
                                        BoutReal ShearFactor = 1.0);

  const Field2D& Rxy() const { return Rxy_m; }
  const Field2D& Bpxy() const { return Bpxy_m; }
  const Field2D& Btxy() const { return Btxy_m; }
  const Field2D& Bxy() const { return Bxy_m; }
  const Field2D& hthe() const { return hthe_m; }
  const FieldMetric& ShearFactor() const { return ShearFactor_m; }
};

#endif //BOUT_TOKAMAK_COORDINATES_HXX
