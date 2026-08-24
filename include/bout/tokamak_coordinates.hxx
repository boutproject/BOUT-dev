#ifndef BOUT_TOKAMAK_COORDINATES_HXX
#define BOUT_TOKAMAK_COORDINATES_HXX

#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
#include <bout/field2d.hxx>
#include <bout/metric_tensor.hxx>

class Mesh;

namespace bout {
/// Variables used in the BOUT++ tokamak coordinate system.
///
/// When created by `set_tokamak_coordinates`, all variables except
/// `I_unnormalised` are normalised.
struct TokamakCoordinates {
  /// Major radius
  FieldMetric Rxy;
  /// Vertical height
  FieldMetric Zxy;
  /// Poloidal magnetic field
  FieldMetric Bpxy;
  /// Toroidal magnetic field
  FieldMetric Btxy;
  /// Total magnetic field
  FieldMetric Bxy;
  /// Poloidal arc length
  FieldMetric hthe;
  /// Integrated shear (normalised)
  FieldMetric I;
  /// Unnormalised integrated shear
  FieldMetric I_unnormalised;
};

/// Read, normalise, calculate, and set the metric components for a BOUT++
/// tokamak coordinate system.
///
/// Sets the cell-centre `Coordinates` on \p mesh
///
/// @param Lbar         Spatial normalisation
/// @param Bbar         Magnetic normalisation
/// @param no_shear     Don't use integrated shear in coordinate system
/// @param shear_factor Integrated shear normalisation
TokamakCoordinates set_tokamak_coordinates(Mesh& mesh, BoutReal Lbar = 1.0,
                                           BoutReal Bbar = 1.0, bool no_shear = false,
                                           BoutReal shear_factor = 1.0);

class TokamakMetricNormaliser : public MetricNormaliser {
public:
  TokamakMetricNormaliser(BoutReal Bnorm, BoutReal rho_s0)
      : Bnorm(Bnorm), rho_s0(rho_s0) {};
  std::optional<BoutReal> g11() const override;
  std::optional<BoutReal> g22() const override;
  std::optional<BoutReal> g33() const override;
  std::optional<BoutReal> g12() const override;
  std::optional<BoutReal> g13() const override;
  std::optional<BoutReal> g23() const override;
  std::optional<BoutReal> dx() const override;
  std::optional<BoutReal> J() const override;
  std::optional<BoutReal> Bxy() const override;

private:
  BoutReal Bnorm;
  BoutReal rho_s0;
};

class FCIMetricNormaliser : public MetricNormaliser {
public:
  FCIMetricNormaliser(BoutReal Bnorm, BoutReal rho_s0) : Bnorm(Bnorm), rho_s0(rho_s0) {};
  std::optional<BoutReal> g() const override;
  std::optional<BoutReal> J() const override;
  std::optional<BoutReal> Bxy() const override;

private:
  BoutReal Bnorm;
  BoutReal rho_s0;
};

std::unique_ptr<MetricNormaliser>
TokamakOrFCIMetricNormaliser(const Mesh* mesh, BoutReal Bnorm, BoutReal rho_s0);

} // namespace bout

#endif //BOUT_TOKAMAK_COORDINATES_HXX
