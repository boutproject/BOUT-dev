#include "../options.hxx"
#include "../invert_parderiv.hxx"

class HeatFluxSNB {
public:
  HeatFluxSNB(Options *options = nullptr) {
    invertpar = InvertPar::Create();

    if (!options) {
      options = &Options::root()["snb"];
    }

    // Read options. Note that the defaults are initialised already
    r = (*options)["r"].doc("Scaling of the electron-electron mean free path")
      .withDefault(r);
    beta_max = (*options)["beta_max"].doc("Maximum energy group to consider (multiple of eT)")
      .withDefault(beta_max);
    ngroups = (*options)["ngroups"].doc("Number of energy groups").withDefault(ngroups);
    
  }
  ~HeatFluxSNB() {
    delete invertpar;
  }

  /// Calculate divergence of heat flux
  /// Te: Electron temperature in eV
  /// Ne: Electron density in m^-3
  Field3D div_heatflux(const Field3D &Te, const Field3D &Ne);
  
private:
  InvertPar *invertpar;
  
  BoutReal Z{1}; // Average ion charge (1 = Hydrogen)
  BoutReal r{2}; // Electron-electron mean free path scaling factor
  BoutReal beta_max{5.0}; // Maximum energy group to consider (multiple of eT)
  int ngroups{20}; // Number of energy groups

  /// Indefinite integral of beta^4 * exp(-beta)
  /// with constant set to zero
  BoutReal int_beta4_exp(BoutReal beta) {
    return - exp(- beta) * (24 + beta * (24 + beta * (12 + beta * (4 + beta))));
  }
  
  /// Integral of beta^4 * exp(-beta) from beta_min to beta_max
  BoutReal groupWeight(BoutReal beta_min, BoutReal beta_max) {
    return int_beta4_exp(beta_max) - int_beta4_exp(beta_min);
  }
};
