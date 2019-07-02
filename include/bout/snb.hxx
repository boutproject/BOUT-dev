#include "../options.hxx"
#include "../invert_parderiv.hxx"

#include <memory>

namespace bout {

/// Calculate heat flux using the Shurtz-Nicolai-Busquet (SNB) model
///
/// Useful references:
///
///   Braginskii equations by R.Fitzpatrick: http://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node35.html
/// 
///   J.P.Brodrick et al 2017: https://doi.org/10.1063/1.5001079 and https://arxiv.org/abs/1704.08963
///
///   Shurtz, Nicolai and Busquet 2000: https://doi.org/10.1063/1.1289512
///
class HeatFluxSNB {
public:
  /// Construct using the options in the "snb" section.
  HeatFluxSNB() : HeatFluxSNB(Options::root()["snb"]) {}

  /// Construct using options in given section.
  explicit HeatFluxSNB(Options &options) {
    invertpar = std::unique_ptr<InvertPar>{InvertPar::Create()};

    // Read options. Note that the defaults are initialised already
    r = options["r"].doc("Scaling of the electron-electron mean free path")
      .withDefault(r);
    beta_max = options["beta_max"].doc("Maximum energy group to consider (multiple of eT)")
      .withDefault(beta_max);
    ngroups = options["ngroups"].doc("Number of energy groups").withDefault(ngroups);
  }
  
  ~HeatFluxSNB() = default;

  HeatFluxSNB(HeatFluxSNB&&) = default;
  HeatFluxSNB& operator=(HeatFluxSNB&&) = default;

  // No copy constructors because invertpar is std::unique_ptr
  HeatFluxSNB(const HeatFluxSNB&) = delete;
  HeatFluxSNB& operator=(const HeatFluxSNB&) = delete;

  /// Calculate divergence of heat flux
  /// Te: Electron temperature in eV
  /// Ne: Electron density in m^-3
  ///
  /// Div_Q_SH_out : An optional output field to store the Spitzer-Harm heat flux
  ///
  /// Returns the divergence of heat flux in units of eV per cubic meter per second
  /// -> multiply by e=1.602e-19 to get Watts per cubic meter.
  Field3D divHeatFlux(const Field3D &Te, const Field3D &Ne, Field3D *Div_Q_SH_out = nullptr);
  
private:
  /// Parallel inversion of tridiagonal matrices
  std::unique_ptr<InvertPar> invertpar{nullptr};
  
  BoutReal Z{1}; ///< Average ion charge (1 = Hydrogen)
  BoutReal r{2}; ///< Electron-electron mean free path scaling factor
  BoutReal beta_max{10.0}; ///< Maximum energy group to consider (multiple of eT)
  int ngroups{40}; ///< Number of energy groups

  /// Indefinite integral of beta^4 * exp(-beta)
  /// with constant set to zero
  BoutReal int_beta4_exp(BoutReal beta) {
    return - exp(- beta) * (24 + beta * (24 + beta * (12 + beta * (4 + beta))));
  }
  
  /// (1/24) * Integral of beta^4 * exp(-beta) from beta_min to beta_max
  BoutReal groupWeight(BoutReal beta_min, BoutReal beta_max) {
    return (1./24) * (int_beta4_exp(beta_max) - int_beta4_exp(beta_min));
  }
};

} // namespace bout
