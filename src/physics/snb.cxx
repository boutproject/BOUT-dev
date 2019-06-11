/// SNB model
///

#include "bout/snb.hxx"
#include "derivs.hxx"
#include "bout/fv_ops.hxx"
#include "bout/constants.hxx"

Field3D HeatFluxSNB::div_heatflux(const Field3D &Te, const Field3D &Ne) {
  Coordinates *coord = Te.getCoordinates();

  Field3D thermal_speed = sqrt(2.*SI::qe  * Te / SI::Me);

  BoutReal Y = SQ(SQ(SI::qe) / (SI::e0 * SI::Me)) / (4 * PI);
  Field3D coulomb_log = 6.6 - 0.5 * log(Ne * 1e-20) + 1.5 * log(Te);

  // Thermal electron-electron mean free path [m]
  Field3D lambda_ee_T = pow(thermal_speed, 3) / (Y * Ne  * coulomb_log);
  Field3D lambda_ei_T = lambda_ee_T / SQ(Z);
  

  // Thermal electron-ion collision time [s]
  Field3D tau_ei_T  = thermal_speed / lambda_ei_T;
  
  // Divergence of Spitzer-Harm heat flux
  Field3D Div_Q_SH = FV::Div_par_K_Grad_par((Ne * SI::qe * Te / SI::Me)
                                            * (0.25 * 3 * sqrt(PI) * tau_ei_T)
                                            * 13.58 * (Z + 0.24) / (Z + 4.2),
                                            Te);
  
  Field3D lambda_ee_Tprime = lambda_ee_T / r;
  Field3D lambda_ei_Tprime = lambda_ei_T * ((Z + 0.25) / (Z + 4.2));
  
  // Loop over energy groups.
  // beta = E / eT is the normalised group energy
  
  BoutReal beta_last = 0.0; // The last beta value calculated. Ths increases through the loop
  BoutReal dbeta = beta_max / ngroups; // Step in beta

  Field3D Div_Q = Div_Q_SH;  // Divergence of heat flux. Corrections added for each group

  for (int i = 0; i < ngroups; i++) {
    BoutReal beta = beta_last + dbeta;
    BoutReal weight = groupWeight(beta_last, beta);

    // Mean free paths for this group
    Field3D lambda_g_ee = SQ(beta) * lambda_ee_Tprime;
    Field3D lambda_g_ei = SQ(beta) * lambda_ei_Tprime;
    
    // Update coefficients in solver
    invertpar->setCoefA(1./lambda_g_ee);   // Constant term

    // The divergence term is implemented as a second derivative and first derivative correction
    Field3D coefB = (-1./3) * lambda_g_ei;
    invertpar->setCoefB(coefB); // Grad2_par2
    invertpar->setCoefE(DDY(coefB * coord->J / coord->g_22) / coord->J); // DDY

    // Solve to get H_g
    Field3D H_g = invertpar->solve((-weight) * Div_Q);

    // Add correction to divergence of heat flux
    // Note: The sum of weight over all groups approaches 24 as beta_max -> infinity
    Div_Q -= weight * Div_Q_SH - H_g / lambda_g_ee;
    
    // move to next group, updating lower limit 
    beta_last = beta;
  }

  return Div_Q;
}
