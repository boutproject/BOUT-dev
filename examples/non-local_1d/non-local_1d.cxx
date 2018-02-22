/*
 * 1D parallel transport problem
 * 
 */

#include <bout.hxx>
#include <bout/boutmain.hxx>

#include <bout/interpolation.hxx>
#include <cmath>
#include <sstream>

#include "non-local_parallel.hxx"

#include "sin-single-sources.cxx"

// #include <iomanip>

#define CUSTOMBCS
// #define CALCULATE_EFFECTIVE_ALPHAE
// #define LOCALHEATFLUX
// #define FLUXLIMITER
#define IONFLUXLIMITER
#define IONVISCOSITYLIMITER

#ifdef CALCULATE_EFFECTIVE_ALPHAE
  #ifdef LOCALHEATFLUX
    throw BoutException("Cannot calculate effective_alpha_e without non-local heat-flux");
  #endif
  FieldPerp effective_alpha_e;
#endif

#ifdef FLUXLIMITER
  #ifndef LOCALHEATFLUX
    throw BoutException("FLUXLIMITER can only be defined if LOCALHEATFLUX is also defined");
  #endif
  BoutReal electron_flux_limiter = 0.3;
  Field3D grad_par_T_electron;
  Field3D qSH_electron;
  Field3D qFS_electron;
  Field3D grad_par_electron_heat_flux;
#endif

#ifdef IONFLUXLIMITER
  BoutReal ion_flux_limiter = 0.1;
  Field3D grad_par_T_ion;
  Field3D qSH_ion;
  Field3D qFS_ion;
  Field3D grad_par_ion_heat_flux;
#endif

#ifdef IONVISCOSITYLIMITER
  BoutReal ion_viscosity_limiter = 0.5;
//   BoutReal ion_viscosity_limiter = 1.;
  Field3D grad_par_viscosity;
  Field3D viscosity;
  Field3D grad_par_V_centre;
  Field3D nTtau_ion_ylow;
  Field3D pi_Brag_ylow;
  Field3D pi_Brag_centre;
#endif

//   BoutReal gamma_factor = 1.;
//   BoutReal gamma_factor = 5./3.;
  BoutReal gamma_factor = 3.;

/********************************************************************************************************************************/

Field3D T_electron; // Electron temperature
// Ambipolarity: Assume that the electron density and velocity are equal to the ion density and velocity
Field3D n_ion; // Ion density
Field3D Vpar_ion; // Ion fluid velocity (parallel component since this is a 1d model)
Field3D T_ion; // Evolve the ion temperature separately (using Braginskii closure for ion heat flux)
Field3D j_parallel; // Parallel current electron_charge*n_ion*(Vpar_ion-Vpar_electron)

Field3D tau_ii; //Ion-ion collision time
Field3D tau_ei; //Electron-ion collision time
Field3D nTtau_ion; //used in calculating the viscosity tensor (parallel components) and ion heat flux
// Field3D VeminusVi;

BoutReal massunit;
BoutReal energyunit;
// Physical constants
BoutReal electron_charge;
BoutReal electron_mass;
BoutReal ion_charge;
BoutReal ion_mass;
BoutReal epsilon_0;
BoutReal logLambda;
int boundary_gradient_smoothing_length;
BoutReal boundary_condition_smoothing_range;

bindex position;

Sources sources;

NonLocalParallel nonlocal_parallel;
#ifndef LOCALHEATFLUX
  Field3D heat_flux_boundary_condition;
#endif

Field3D ratio;

/**********************************************************************************************************************************************/

int physics_init(bool restarting) {
  
  boundary_gradient_smoothing_length = mesh->GlobalNy/256; // Number of grid points to average over for the gradient used to extrapolate to the guard cells
  boundary_condition_smoothing_range = BoutReal(mesh->GlobalNy)/64.;
	
  // Physical constants
  electron_charge = -1.602176565e-19;	//C
    // Scale factors for units relative to SI (temperatures will be in energy units). All units apart from mass and energy are SI (i.e. m, s, C).
    massunit = abs(electron_charge); //the internal mass unit is 1.602176565e-19kg so that energy will be in eV
    energyunit = massunit; //internal energy unit is 1.602176565e-19J so that energies/temperatures are stored in eV
  electron_mass = 9.10938291e-31/massunit; //kg
  ion_charge = 1.602176565e-19; //C
  ion_mass =  3.34358348e-27/massunit;// kg
  epsilon_0 = 8.85418781762039e-12/pow(massunit,-1);	//C^2 s^2 kg^-1 m^-3
  logLambda = 16.0;	//Coulomb Logarithm --- logLambda=16.1 for (n = 10^18 m^-3) and (T = 100eV)
  //Recall: logLambda~25.3-1.15*log10(n/cm^3)+2.3*log10(T_e/eV) for T_e>50eV; logLambda~23.4-1.15*log10(n/cm^3)+3.45*log10(T_e/eV) for T_e<50eV -- Hazeltine+Meiss (2003)

  // Get the options for the model
  Options *options = Options::getRoot()->getSection("conduction");
  
  // Get the options for the sources
  options = Options::getRoot()->getSection("sources");
//   OPTION(options, sources.lower_source_yindex, 0);
//   OPTION(options, sources.upper_source_yindex, 0);
  OPTION(options, sources.source_length, 0);
  OPTION(options, sources.particle_amplitude ,0);
  OPTION(options, sources.electron_heat_amplitude ,0);
  OPTION(options, sources.ion_heat_amplitude ,0);
  OPTION(options, sources.on_time ,0);
  OPTION(options, sources.off_time ,0);
  OPTION(options, sources.particle_transient_amplitude ,0);
  OPTION(options, sources.electron_heat_transient_amplitude ,0);
  OPTION(options, sources.ion_heat_transient_amplitude ,0);
  
  sources.initialise();
  
  // Initialise nonlocal_parallel object
//   #ifndef LOCALHEATFLUX
  nonlocal_parallel.initialise(electron_charge, electron_mass, ion_mass, epsilon_0, logLambda, mesh->StaggerGrids, gamma_factor);
//   #endif
  #ifndef LOCALHEATFLUX
    heat_flux_boundary_condition.setLocation(CELL_YLOW);
    heat_flux_boundary_condition = 0.;
  #endif
  tau_ei = 0.;
  tau_ii = 0.;
  
  // Set non-default cell locations
  Vpar_ion.setLocation(CELL_YLOW); // Staggered relative to n_ion, etc.
  j_parallel.setLocation(CELL_YLOW);
  ratio.setLocation(CELL_YLOW);
  #ifdef IONVISCOSITYLIMITER
    grad_par_viscosity.setLocation(CELL_YLOW);
    pi_Brag_ylow.setLocation(CELL_YLOW);
  #endif
    
  // Tell BOUT++ which fields to evolve
  SOLVE_FOR(T_electron);
  SOLVE_FOR(n_ion);
  SOLVE_FOR(Vpar_ion);
  SOLVE_FOR(T_ion);
  SOLVE_FOR(j_parallel);
  
//   dump.add(mesh->dy,"dy",0);
//   dump.add(mesh->g_22,"g_22",0);
  #ifdef CALCULATE_HEATFLUX
    dump.add(nonlocal_parallel.electron_heat_flux,"heat_flux",1);
  #endif
  #ifdef CALCULATE_VISCOSITY
    dump.add(nonlocal_parallel.electron_viscosity,"viscosity",1);
  #endif
  #ifdef CALCULATE_FRICTION
    dump.add(nonlocal_parallel.electron_friction,"friction",1);
  #endif
  dump.add(ratio,"local_non-local_ratio",1);
  #ifdef CALCULATE_EFFECTIVE_ALPHAE
    effective_alpha_e = 0.;
    dump.add(effective_alpha_e[2][0],"effective_alpha_e",1);
  #endif
  
  output<<endl;
  
  #ifdef LOCALHEATFLUX
  output<<"Using local electron heat flux"<<endl;
  #endif
  
  #ifdef FLUXLIMITER
    output<<"Electron Flux Limiter is "<<electron_flux_limiter<<endl;
  #endif
    
  #ifdef IONFLUXLIMITER
    output<<"Ion Flux Limiter is "<<ion_flux_limiter<<endl;
  #endif
    
  #ifdef IONVISCOSITYLIMITER
    output<<"Ion Viscosity Limiter is "<<ion_viscosity_limiter<<endl;
  #endif
    
  output<<"gamma is "<<gamma_factor<<endl;
  
  output<<endl;
  
  return 0;
}

/***********************************************************************************************************************************************/

int physics_run(BoutReal t) {
  
  output<<"\r"<<t<<std::flush;
  
  mesh->communicate(T_electron, n_ion, Vpar_ion, T_ion);// Communicate guard cells
  
  // Done like this so we can deal with T_ion<0 at mesh->yend (which is a 'guard cell' for the staggered grid case) without throwing an exception
  for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
    for (int jz=0; jz<mesh->LocalNz; jz++)
      for (int jy=0; jy<mesh->LocalNy; jy++) {
	tau_ii(jx,jy,jz) = 3 * pow(PI,1.5) * pow(epsilon_0,2) * sqrt(ion_mass) * pow(2.,1.5) * pow(T_ion(jx,jy,jz),1.5) / n_ion(jx,jy,jz) / pow(ion_charge,4) / logLambda;
	tau_ei(jx,jy,jz) = 3 * pow(PI,1.5) * pow(epsilon_0,2) * sqrt(electron_mass) * pow(2.,1.5) * pow(T_electron(jx,jy,jz),1.5) / n_ion(jx,jy,jz) / pow(electron_charge,2) / pow(ion_charge,2) / logLambda;
      }

  #ifdef CUSTOMBCS
  // Fix the temperature to continue with the gradient given by the 4-point forward/backward difference estimate // and electron temperature gradient to give electron heat flux = 5.0 T_electron n_ion Vpar_ion at the boundaries, assuming that it were highly collisional (this is not actually used because the boundary condition is applied to the heat flux explicitly)
  
    for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	position.jx = rlow.ind;
	position.jy = mesh->ystart;
	position.jz = jz;
	calc_index(&position);
	BoutReal n_here = nonlocal_parallel.interp_to_point_YLOW(n_ion,position);
	BoutReal Te_here = nonlocal_parallel.interp_to_point_YLOW(T_electron,position);
	BoutReal Ti_here = nonlocal_parallel.interp_to_point_YLOW(T_ion,position);
	BoutReal tauei_here = nonlocal_parallel.interp_to_point_YLOW(tau_ei,position);
	#ifdef LOCALHEATFLUX
	  BoutReal sheath_potential = 0.5*Te_here*log(2*PI*electron_mass/ion_mass*(1+gamma_factor*Ti_here/Te_here));
	  #ifdef FLUXLIMITER
	    // Include the flux limiter in the boundary condition: i.e. find the gradient needed to give the right boundary heat-flux AFTER applying the flux-limiter
	    BoutReal qbc_over_n = -((2.0-2.5)*Te_here-sheath_potential)
				    *sqrt((Te_here+gamma_factor*Ti_here)/ion_mass); // -2.5*Te so that we subtract off the convective heat flux to leave just the conductive heat flux, not the total
	    BoutReal qBrag_over_n = 1./(1./qbc_over_n-1./electron_flux_limiter/( -sqrt(2./electron_mass)*pow(Te_here,1.5) ));
	    BoutReal gradient_T_electron = qBrag_over_n
					      /(-3.16*Te_here*tauei_here/electron_mass)
	      *coord->dy(rlow.ind,mesh->ystart)*sqrt(coord->g_22(rlow.ind,mesh->ystart));
#else
	    BoutReal gradient_T_electron = -((2.0-2.5)*Te_here-sheath_potential)
	      *sqrt((Te_here+gamma_factor*Ti_here)/ion_mass)
	      /(-3.16*Te_here*tauei_here/electron_mass)
	      *coord->dy(rlow.ind,mesh->ystart)*sqrt(coord->g_22(rlow.ind,mesh->ystart)); // -2.5*Te so that we subtract off the convective heat flux to leave just the conductive heat flux, not the total
	  #endif
	#else
  //       BoutReal gradient_T_electron = (-11.*T_electron[rlow.ind][mesh->ystart][jz] + 18.*T_electron[rlow.ind][mesh->ystart+1][jz] - 9.*T_electron[rlow.ind][mesh->ystart+2][jz] + 2.*T_electron[rlow.ind][mesh->ystart+3][jz]) / 6. / mesh->dy[rlow.ind][mesh->ystart] / sqrt((mesh->g_22[rlow.ind][mesh->ystart]+mesh->g_22[rlow.ind][mesh->ystart+1]+mesh->g_22[rlow.ind][mesh->ystart+2]+mesh->g_22[rlow.ind][mesh->ystart+3])/4.);
	    BoutReal gradient_T_electron = (-T_electron(rlow.ind,mesh->ystart,jz) + T_electron(rlow.ind,mesh->ystart+boundary_gradient_smoothing_length,jz))/BoutReal(boundary_gradient_smoothing_length);
	#endif
  //       BoutReal gradient_n = (-11.*n_ion[rlow.ind][mesh->ystart][jz] + 18.*n_ion[rlow.ind][mesh->ystart+1][jz] - 9.*n_ion[rlow.ind][mesh->ystart+2][jz] + 2.*n_ion[rlow.ind][mesh->ystart+3][jz]) / 6. / mesh->dy[rlow.ind][mesh->ystart] / sqrt((mesh->g_22[rlow.ind][mesh->ystart]+mesh->g_22[rlow.ind][mesh->ystart+1]+mesh->g_22[rlow.ind][mesh->ystart+2]+mesh->g_22[rlow.ind][mesh->ystart+3])/4.);
	    BoutReal gradient_n = (-n_ion(rlow.ind,mesh->ystart,jz) + n_ion(rlow.ind,mesh->ystart+boundary_gradient_smoothing_length,jz))/BoutReal(boundary_gradient_smoothing_length);
  //       BoutReal gradient_T_ion = (-11.*T_ion[rlow.ind][mesh->ystart][jz] + 18.*T_ion[rlow.ind][mesh->ystart+1][jz] - 9.*T_ion[rlow.ind][mesh->ystart+2][jz] + 2.*T_ion[rlow.ind][mesh->ystart+3][jz]) / 6. / mesh->dy[rlow.ind][mesh->ystart] / sqrt((mesh->g_22[rlow.ind][mesh->ystart]+mesh->g_22[rlow.ind][mesh->ystart+1]+mesh->g_22[rlow.ind][mesh->ystart+2]+mesh->g_22[rlow.ind][mesh->ystart+3])/4.);
	    BoutReal gradient_T_ion = (-T_ion(rlow.ind,mesh->ystart,jz) + T_ion(rlow.ind,mesh->ystart+boundary_gradient_smoothing_length,jz))/BoutReal(boundary_gradient_smoothing_length);
	for (int jy=mesh->ystart-1; jy>=0; jy--) {
	  n_ion(rlow.ind,jy,jz) = n_ion(rlow.ind,jy+1,jz) - gradient_n;
  // 	n_ion(rlow.ind,jy,jz) = n_ion(rlow.ind,jy+1,jz);
	  T_ion(rlow.ind,jy,jz) = T_ion(rlow.ind,jy+1,jz) - gradient_T_ion;
	  #ifndef LOCALHEATFLUX
	  T_electron(rlow.ind,jy,jz) = T_electron(rlow.ind,jy+1,jz) - gradient_T_electron;
	  #endif
	}
	#ifdef LOCALHEATFLUX
	// Set it up so that the fourth-order central finite difference derivative at CELL_YLOW of mesh->ystart is gradient_T_electron
	T_electron(rlow.ind,mesh->ystart-1,jz) = T_electron(rlow.ind,mesh->ystart,jz) - gradient_T_electron;
	T_electron(rlow.ind,mesh->ystart-2,jz) = T_electron(rlow.ind,mesh->ystart+1,jz) - 3.*gradient_T_electron;
	  for (int jy=mesh->ystart-3; jy>=0; jy--)
	    T_electron(rlow.ind,jy,jz) = T_electron(rlow.ind,jy+1,jz) - gradient_T_electron;
	#endif
      }

    if (mesh->StaggerGrids)
      for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
	for (int jz=0; jz<mesh->LocalNz; jz++) {
	  position.jx = rup.ind;
	  position.jy = mesh->yend;
	  position.jz = jz;
	  calc_index(&position);
	  BoutReal n_here = nonlocal_parallel.interp_to_point_YLOW(n_ion,position);
	  BoutReal Te_here = nonlocal_parallel.interp_to_point_YLOW(T_electron,position);
	  BoutReal Ti_here = nonlocal_parallel.interp_to_point_YLOW(T_ion,position);
	  BoutReal tauei_here = nonlocal_parallel.interp_to_point_YLOW(tau_ei,position);
	  #ifdef LOCALHEATFLUX
	    BoutReal sheath_potential = 0.5*Te_here*log(2*PI*electron_mass/ion_mass*(1+gamma_factor*Ti_here/Te_here));
	    #ifdef FLUXLIMITER
	      // Include the flux limiter in the boundary condition: i.e. find the gradient needed to give the right boundary heat-flux AFTER applying the flux-limiter
	      BoutReal qbc_over_n = ((2.0-2.5)*Te_here-sheath_potential)
				      *sqrt((Te_here+gamma_factor*Ti_here)/ion_mass); // -2.5*Te so that we subtract off the convective heat flux to leave just the conductive heat flux, not the total
	      BoutReal qBrag_over_n = 1./(1./qbc_over_n-1./electron_flux_limiter/( sqrt(2./electron_mass)*pow(Te_here,1.5) ));
	      BoutReal gradient_T_electron = qBrag_over_n
						/(-3.16*Te_here*tauei_here/electron_mass)
						*mesh->dy[rup.ind][mesh->yend-1]*sqrt(mesh->g_22[rup.ind][mesh->yend-1]);
	    #else
	      BoutReal gradient_T_electron = ((2.0-2.5)*Te_here-sheath_potential)
						*sqrt((Te_here+gamma_factor*Ti_here)/ion_mass)
						/(-3.16*Te_here*tauei_here/electron_mass)
						*mesh->dy[rup.ind][mesh->yend-1]*sqrt(mesh->g_22[rup.ind][mesh->yend-1]); // -2.5*Te so that we subtract off the convective heat flux to leave just the conductive heat flux, not the total
	    #endif
	  #else
  //         BoutReal gradient_T_electron = (11.*T_electron[rup.ind][mesh->yend-1][jz] - 18.*T_electron[rup.ind][mesh->yend-2][jz] + 9.*T_electron[rup.ind][mesh->yend-3][jz] - 2.*T_electron[rup.ind][mesh->yend-4][jz]) / 6. / mesh->dy[rup.ind][mesh->yend] / sqrt((mesh->g_22[rup.ind][mesh->yend-1]+mesh->g_22[rup.ind][mesh->yend-2]+mesh->g_22[rup.ind][mesh->yend-3]+mesh->g_22[rup.ind][mesh->yend-4])/4.);
	    BoutReal gradient_T_electron = (-T_electron[rup.ind][mesh->yend-1-boundary_gradient_smoothing_length][jz] + T_electron[rup.ind][mesh->yend-1][jz])/BoutReal(boundary_gradient_smoothing_length);
	  #endif
  // 	BoutReal gradient_n = (11.*n_ion[rup.ind][mesh->yend-1][jz] - 18.*n_ion[rup.ind][mesh->yend-2][jz] + 9.*n_ion[rup.ind][mesh->yend-3][jz] - 2.*n_ion[rup.ind][mesh->yend-4][jz]) / 6. / mesh->dy[rup.ind][mesh->yend] / sqrt((mesh->g_22[rup.ind][mesh->yend-1]+mesh->g_22[rup.ind][mesh->yend-2]+mesh->g_22[rup.ind][mesh->yend-3]+mesh->g_22[rup.ind][mesh->yend-4])/4.);
	  BoutReal gradient_n = (-n_ion[rup.ind][mesh->yend-1-boundary_gradient_smoothing_length][jz] + n_ion[rup.ind][mesh->yend-1][jz])/BoutReal(boundary_gradient_smoothing_length);
  // 	BoutReal gradient_T_ion = (11.*T_ion[rup.ind][mesh->yend-1][jz] - 18.*T_ion[rup.ind][mesh->yend-2][jz] + 9.*T_ion[rup.ind][mesh->yend-3][jz] - 2.*T_ion[rup.ind][mesh->yend-4][jz]) / 6. / mesh->dy[rup.ind][mesh->yend] / sqrt((mesh->g_22[rup.ind][mesh->yend-1]+mesh->g_22[rup.ind][mesh->yend-2]+mesh->g_22[rup.ind][mesh->yend-3]+mesh->g_22[rup.ind][mesh->yend-4])/4.);
	  BoutReal gradient_T_ion = (-T_ion[rup.ind][mesh->yend-1-boundary_gradient_smoothing_length][jz] + T_ion[rup.ind][mesh->yend-1][jz])/BoutReal(boundary_gradient_smoothing_length);
	  for (int jy=mesh->yend; jy<mesh->LocalNy; jy++) {
	    n_ion[rup.ind][jy][jz] = n_ion[rup.ind][jy-1][jz] + gradient_n;
	    T_ion[rup.ind][jy][jz] = T_ion[rup.ind][jy-1][jz] + gradient_T_ion;
	    #ifndef LOCALHEATFLUX
	      T_electron[rup.ind][jy][jz] = T_electron[rup.ind][jy-1][jz] + gradient_T_electron;
	    #endif
	  }
	  #ifdef LOCALHEATFLUX
	  // Set it up so that the fourth-order central finite difference derivative at CELL_YLOW of mesh->ystart is gradient_T_electron
	    T_electron[rup.ind][mesh->yend][jz] = T_electron[rup.ind][mesh->yend-1][jz] + gradient_T_electron;
	    T_electron[rup.ind][mesh->yend+1][jz] = T_electron[rup.ind][mesh->yend-2][jz] + 3.*gradient_T_electron;
	    for (int jy=mesh->yend+2; jy<mesh->LocalNy; jy++)
	      T_electron[rup.ind][jy][jz] = T_electron[rup.ind][jy-1][jz] + gradient_T_electron;
	  #endif
	}
    else
      for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
	for (int jz=0; jz<mesh->LocalNz; jz++) {
	  #ifdef LOCALHEATFLUX
	    BoutReal sheath_potential = 0.5*T_electron[rup.ind][mesh->yend][jz]*log(2*PI*electron_mass/ion_mass*(1+gamma_factor*T_ion[rup.ind][mesh->yend][jz]/T_electron[rup.ind][mesh->yend][jz]));
	    #ifdef FLUXLIMITER
	      // Include the flux limiter in the boundary condition: i.e. find the gradient needed to give the right boundary heat-flux AFTER applying the flux-limiter
	      BoutReal qbc_over_n = ((2.0-2.5)*T_electron[rup.ind][mesh->yend][jz]-sheath_potential)
				      *sqrt((T_electron[rup.ind][mesh->yend][jz]+gamma_factor*T_ion[rup.ind][mesh->yend][jz])/ion_mass); // -2.5*Te so that we subtract off the convective heat flux to leave just the conductive heat flux, not the total
	      BoutReal qBrag_over_n = 1./(1./qbc_over_n-1./electron_flux_limiter/( sqrt(2./electron_mass)*pow(T_electron[rup.ind][mesh->yend][jz],1.5) ));
	      BoutReal gradient_T_electron = qBrag_over_n
						/(-3.16*T_electron[rup.ind][mesh->yend][jz]*tau_ei[rup.ind][mesh->yend][jz]/electron_mass)
						*mesh->dy[rup.ind][mesh->yend-1]*sqrt(mesh->g_22[rup.ind][mesh->yend-1]);
	    #else
	      BoutReal gradient_T_electron = ((2.0-2.5)*T_electron[rup.ind][mesh->yend][jz]-sheath_potential)
						*sqrt((T_electron[rup.ind][mesh->yend][jz]+gamma_factor*T_ion[rup.ind][mesh->yend][jz])/ion_mass)
						/(-3.16*T_electron[rup.ind][mesh->yend][jz]*tau_ei[rup.ind][mesh->yend][jz]/electron_mass)
						*mesh->dy[rup.ind][mesh->yend]*sqrt(mesh->g_22[rup.ind][mesh->yend]); // -2.5*Te so that we subtract off the convective heat flux to leave just the conductive heat flux, not the total
	    #endif
	  #else
  //         BoutReal gradient_T_electron = (11.*T_electron[rup.ind][mesh->yend][jz] - 18.*T_electron[rup.ind][mesh->yend-1][jz] + 9.*T_electron[rup.ind][mesh->yend-2][jz] - 2.*T_electron[rup.ind][mesh->yend-3][jz]) / 6. / mesh->dy[rup.ind][mesh->yend] / sqrt((mesh->g_22[rup.ind][mesh->yend]+mesh->g_22[rup.ind][mesh->yend-1]+mesh->g_22[rup.ind][mesh->yend-2]+mesh->g_22[rup.ind][mesh->yend-3])/4.);
	    BoutReal gradient_T_electron = (-T_electron[rup.ind][mesh->yend-boundary_gradient_smoothing_length][jz] + T_electron[rup.ind][mesh->yend][jz])/BoutReal(boundary_gradient_smoothing_length);
	  #endif
  // 	BoutReal gradient_n = (11.*n_ion[rup.ind][mesh->yend][jz] - 18.*n_ion[rup.ind][mesh->yend-1][jz] + 9.*n_ion[rup.ind][mesh->yend-2][jz] - 2.*n_ion[rup.ind][mesh->yend-3][jz]) / 6. / mesh->dy[rup.ind][mesh->yend] / sqrt((mesh->g_22[rup.ind][mesh->yend]+mesh->g_22[rup.ind][mesh->yend-1]+mesh->g_22[rup.ind][mesh->yend-2]+mesh->g_22[rup.ind][mesh->yend-3])/4.);
	  BoutReal gradient_n = (-n_ion[rup.ind][mesh->yend-boundary_gradient_smoothing_length][jz] + n_ion[rup.ind][mesh->yend][jz])/BoutReal(boundary_gradient_smoothing_length);
  // 	BoutReal gradient_T_ion = (11.*T_ion[rup.ind][mesh->yend][jz] - 18.*T_ion[rup.ind][mesh->yend-1][jz] + 9.*T_ion[rup.ind][mesh->yend-2][jz] - 2.*T_ion[rup.ind][mesh->yend-3][jz]) / 6. / mesh->dy[rup.ind][mesh->yend] / sqrt((mesh->g_22[rup.ind][mesh->yend]+mesh->g_22[rup.ind][mesh->yend-1]+mesh->g_22[rup.ind][mesh->yend-2]+mesh->g_22[rup.ind][mesh->yend-3])/4.);
	  BoutReal gradient_T_ion = (-T_ion[rup.ind][mesh->yend-boundary_gradient_smoothing_length][jz] + T_ion[rup.ind][mesh->yend][jz])/BoutReal(boundary_gradient_smoothing_length);
	  for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++) {
	    n_ion[rup.ind][jy][jz] = n_ion[rup.ind][jy-1][jz] + gradient_n;
	    T_ion[rup.ind][jy][jz] = T_ion[rup.ind][jy-1][jz] + gradient_T_ion;
	    #ifndef LOCALHEATFLUX
	      T_electron[rup.ind][jy][jz] = T_electron[rup.ind][jy-1][jz] + gradient_T_electron;
	    #endif
	  }
	  #ifdef LOCALHEATFLUX
	  // Set it up so that the fourth-order central finite difference derivative at CELL_YLOW of mesh->ystart is gradient_T_electron
	    T_electron[rup.ind][mesh->yend+1][jz] = T_electron[rup.ind][mesh->yend][jz] + gradient_T_electron;
	    T_electron[rup.ind][mesh->yend+2][jz] = T_electron[rup.ind][mesh->yend-1][jz] + 3.*gradient_T_electron;
	    for (int jy=mesh->yend+3; jy<mesh->LocalNy; jy++)
	      T_electron[rup.ind][jy][jz] = T_electron[rup.ind][jy-1][jz] + gradient_T_electron;
	  #endif
	}
    // Enforce Vpar_ion=c_s at the boundaries, but only if it tends to increase the magnitude of the velocity (i.e. allow (hopefully temporary) supersonic velocities).
    // Apply the test at the point just inside the boundary so that the boundary point never needs to be evolved by the solver.
    for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	BoutReal boundarygradient;
	position.jx=rlow.ind;
	position.jy=mesh->ystart;
	position.jz=jz;
	calc_index(&position);
	BoutReal boundary_value_Vpar = -sqrt(nonlocal_parallel.interp_to_point_YLOW(T_electron,position)+gamma_factor*nonlocal_parallel.interp_to_point_YLOW(T_ion,position))/sqrt(ion_mass);
	Vpar_ion[rlow.ind][mesh->ystart][jz] = boundary_value_Vpar;
  //     BoutReal boundarygradient = (-11.*Vpar_ion[rlow.ind][mesh->ystart][jz] + 18.*Vpar_ion[rlow.ind][mesh->ystart+1][jz] - 9.*Vpar_ion[rlow.ind][mesh->ystart+2][jz] + 2.*Vpar_ion[rlow.ind][mesh->ystart+3][jz]) / 6. / mesh->dy[rlow.ind][mesh->ystart] / sqrt((mesh->g_22[rlow.ind][mesh->ystart]+mesh->g_22[rlow.ind][mesh->ystart+1]+mesh->g_22[rlow.ind][mesh->ystart+2]+mesh->g_22[rlow.ind][mesh->ystart+3])/4.);
	boundarygradient = (-Vpar_ion[rlow.ind][mesh->ystart][jz] + Vpar_ion[rlow.ind][mesh->ystart+boundary_gradient_smoothing_length][jz])/BoutReal(boundary_gradient_smoothing_length);
	for (int jy=mesh->ystart-1; jy>=0; jy--)
	  Vpar_ion[rlow.ind][jy][jz] = Vpar_ion[rlow.ind][jy+1][jz] - boundarygradient;
  // 	Vpar_ion[rlow.ind][jy][jz] = Vpar_ion[rlow.ind][jy+1][jz];
      }
    for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	BoutReal boundarygradient;
	position.jx=rup.ind;
	position.jy=mesh->yend;
	position.jz=jz;
	calc_index(&position);
	BoutReal boundary_value_Vpar = sqrt(nonlocal_parallel.interp_to_point_YLOW(T_electron,position)+gamma_factor*nonlocal_parallel.interp_to_point_YLOW(T_ion,position))/sqrt(ion_mass);
	Vpar_ion[rup.ind][mesh->yend][jz] = boundary_value_Vpar;
  //     BoutReal boundarygradient = (11.*Vpar_ion[rup.ind][mesh->yend][jz] - 18.*Vpar_ion[rup.ind][mesh->yend-1][jz] + 9.*Vpar_ion[rup.ind][mesh->yend-2][jz] - 2.*Vpar_ion[rup.ind][mesh->yend-3][jz]) / 6. / mesh->dy[rup.ind][mesh->yend] / sqrt((mesh->g_22[rup.ind][mesh->yend]+mesh->g_22[rup.ind][mesh->yend-1]+mesh->g_22[rup.ind][mesh->yend-2]+mesh->g_22[rup.ind][mesh->yend-3])/4.);
	boundarygradient = (-Vpar_ion[rup.ind][mesh->yend-boundary_gradient_smoothing_length][jz] + Vpar_ion[rup.ind][mesh->yend][jz])/BoutReal(boundary_gradient_smoothing_length);
	for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++)
	  Vpar_ion[rup.ind][jy][jz] = Vpar_ion[rup.ind][jy-1][jz] + boundarygradient;
  // 	Vpar_ion[rup.ind][jy][jz] = Vpar_ion[rup.ind][jy-1][jz];
      }
  #endif
  // The boundary condition on the temperature gradient sometimes makes the temperature in the guard cells negative.
  // If this happens set the temperature to zero there and also replace tau with an extrapolated version using forward/backward finite difference derivatives
  for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
    for (int jz=0; jz<mesh->LocalNz; jz++) 
      for (int jy=mesh->ystart-1; jy>=0; jy--) {
	if (T_ion[rlow.ind][jy][jz]<0.) {
	  T_ion[rlow.ind][jy][jz] = 0.;
	  tau_ii[rlow.ind][jy][jz] = tau_ii[rlow.ind][jy+1][jz]-(-11.*tau_ii[rlow.ind][jy+1][jz] + 18.*tau_ii[rlow.ind][jy+2][jz] - 9.*tau_ii[rlow.ind][jy+3][jz] + 2.*tau_ii[rlow.ind][jy+4][jz])/ 6.;
	}
	if (T_electron[rlow.ind][jy][jz]<0.) {
	  T_electron[rlow.ind][jy][jz] = 0.;
	  tau_ei[rlow.ind][jy][jz] = tau_ei[rlow.ind][jy+1][jz]-(-11.*tau_ei[rlow.ind][jy+1][jz] + 18.*tau_ei[rlow.ind][jy+2][jz] - 9.*tau_ei[rlow.ind][jy+3][jz] + 2.*tau_ei[rlow.ind][jy+4][jz])/ 6.;
	}
    }
  if (mesh->StaggerGrids)
    for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
      for (int jz=0; jz<mesh->LocalNz; jz++)
	for (int jy=mesh->yend; jy<mesh->LocalNy; jy++) {
	  if (T_ion[rup.ind][jy][jz]<0.) {
	    T_ion[rup.ind][jy][jz] = 0.;
	    tau_ii[rup.ind][jy][jz] = tau_ii[rup.ind][jy-1][jz] + (11.*tau_ii[rup.ind][jy-1][jz] - 18.*tau_ii[rup.ind][jy-2][jz] + 9.*tau_ii[rup.ind][jy-3][jz] - 2.*tau_ii[rup.ind][jy-4][jz])/ 6.;
	  }
	  if (T_electron[rup.ind][jy][jz]<0.) {
	    T_electron[rup.ind][jy][jz] = 0.;
	    tau_ei[rup.ind][jy][jz] = tau_ei[rup.ind][jy-1][jz] + (11.*tau_ei[rup.ind][jy-1][jz] - 18.*tau_ei[rup.ind][jy-2][jz] + 9.*tau_ei[rup.ind][jy-3][jz] - 2.*tau_ei[rup.ind][jy-4][jz])/ 6.;
	  }
      }
  else
    for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
      for (int jz=0; jz<mesh->LocalNz; jz++)
	for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++) {
	  if (T_ion(rup.ind,jy,jz)<0.) {
	    T_ion(rup.ind,jy,jz) = 0.;
	    tau_ii(rup.ind,jy,jz) = tau_ii(rup.ind,jy-1,jz) + (11.*tau_ii(rup.ind,jy-1,jz) - 18.*tau_ii(rup.ind,jy-2,jz) + 9.*tau_ii(rup.ind,jy-3,jz) - 2.*tau_ii(rup.ind,jy-4,jz))/ 6.;
	  }
	  if (T_electron(rup.ind,jy,jz)<0.) {
	    T_electron(rup.ind,jy,jz) = 0.;
	    tau_ei(rup.ind,jy,jz) = tau_ei(rup.ind,jy-1,jz) + (11.*tau_ei(rup.ind,jy-1,jz) - 18.*tau_ei(rup.ind,jy-2,jz) + 9.*tau_ei(rup.ind,jy-3,jz) - 2.*tau_ei(rup.ind,jy-4,jz))/ 6.;
	  }
      }
  
  nTtau_ion = n_ion * T_ion * sqrt(2) * tau_ii;
  
  #ifndef LOCALHEATFLUX
    for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	position.jx=rlow.ind;
	position.jy=mesh->ystart;
	position.jz=jz;
	calc_index(&position);
	BoutReal Te_here = nonlocal_parallel.interp_to_point_YLOW(T_electron,position);
	BoutReal V_here = Vpar_ion[rlow.ind][mesh->ystart][jz];
	BoutReal n_here = nonlocal_parallel.interp_to_point_YLOW(n_ion,position);
	BoutReal sheath_potential = 0.5*Te_here*log(2*PI*electron_mass*pow(V_here,2)/Te_here);
	heat_flux_boundary_condition[rlow.ind][mesh->ystart][jz] = ((2.0-2.5)*Te_here-sheath_potential)*n_here*V_here; // -2.5*Te so that we subtract off the convective heat flux to leave just the conductive heat flux, not the total
      }
  
    for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	position.jx=rup.ind;
	position.jy=mesh->yend;
	position.jz=jz;
	calc_index(&position);
	BoutReal Te_here = nonlocal_parallel.interp_to_point_YLOW(T_electron,position);
	BoutReal V_here = Vpar_ion[rup.ind][mesh->yend][jz];
	BoutReal n_here = nonlocal_parallel.interp_to_point_YLOW(n_ion,position);
	BoutReal sheath_potential = 0.5*Te_here*log(2*PI*electron_mass*pow(V_here,2)/Te_here);
	heat_flux_boundary_condition[rup.ind][mesh->yend][jz] = ((2.0-2.5)*Te_here-sheath_potential)*n_here*V_here; // -2.5*Te so that we subtract off the convective heat flux to leave just the conductive heat flux, not the total
      }
    
//     VeminusVi = -j_parallel/n_ion/electron_charge;
    Field3D zero = 0.;
    nonlocal_parallel.calculate_nonlocal_closures(n_ion, T_electron, Vpar_ion+j_parallel/electron_charge/interp_to(n_ion,CELL_YLOW), j_parallel, heat_flux_boundary_condition, zero);
//     nonlocal_parallel.calculate_nonlocal_closures(n_ion, T_electron, Vpar_ion+VeminusVi, VeminusVi, heat_flux_boundary_condition);
  
    #ifdef CUSTOMBCS
      nonlocal_parallel.set_boundary_gradients();
    #else
      nonlocal_parallel.set_neumann_boundary_conditions();
    #endif
  #endif

  #ifndef LOCALHEATFLUX
    try {
      ratio = (-3.16*interp_to(n_ion*T_electron*tau_ei/electron_mass,CELL_YLOW)*Grad_par(T_electron,CELL_YLOW))/nonlocal_parallel.electron_heat_flux;
    }
    catch (BoutException error) {
      // Don't really care if there are errors in the calculation of ratio as it is only a diagnostic variable: IGNORE
    }
  #endif
  
  #ifdef CALCULATE_EFFECTIVE_ALPHAE
    Field3D q_Brag = -3.16*interp_to(n_ion*T_electron*tau_ei,CELL_YLOW)/electron_mass*Grad_par(T_electron,CELL_YLOW);
    Field3D q_FS = interp_to(sqrt(2./electron_mass)*n_ion*(T_electron^1.5),CELL_YLOW);
    nonlocal_parallel.mean_over_y( 1./q_FS/abs( 1./nonlocal_parallel.electron_heat_flux-1./q_Brag) , effective_alpha_e);
  #endif
  
  #if defined FLUXLIMITER
    grad_par_T_electron = Grad_par(T_electron,CELL_CENTRE);
    qSH_electron = -3.16/electron_mass*n_ion*T_electron*tau_ei*grad_par_T_electron;
    qFS_electron = sqrt(2./electron_mass)*n_ion*(T_electron^1.5);
    grad_par_electron_heat_flux = electron_flux_limiter / ((abs(qSH_electron)/qFS_electron+electron_flux_limiter)^2)
			      * (qSH_electron*abs(qSH_electron)/(qFS_electron^2)*sqrt(2./electron_mass)*Grad_par(n_ion*(T_electron^1.5),CELL_CENTRE)
				  - electron_flux_limiter*3.16/electron_mass*(Grad_par(n_ion*T_electron*tau_ei,CELL_CENTRE)*grad_par_T_electron + n_ion*T_electron*tau_ei*Grad2_par2(T_electron,CELL_CENTRE)));
    ddt(T_electron) = -Vpar_Grad_par(Vpar_ion,T_electron,CELL_CENTRE) 
		      - 2./3.*T_electron*Grad_par(Vpar_ion,CELL_CENTRE) 
		      - 2./3./n_ion*grad_par_electron_heat_flux
		      + 2./3./n_ion * sources.electron_heat_source(t) - T_electron/n_ion*sources.particle_source(t);
  #elif defined LOCALHEATFLUX
    ddt(T_electron) = -Vpar_Grad_par(Vpar_ion,T_electron,CELL_CENTRE) 
			- 2./3.*T_electron*Grad_par(Vpar_ion,CELL_CENTRE) 
			- 2./3./n_ion*(-3.16/electron_mass*(Grad_par(n_ion*T_electron*tau_ei,CELL_CENTRE)*Grad_par(T_electron,CELL_CENTRE)+n_ion*T_electron*tau_ei*Grad2_par2(T_electron,CELL_CENTRE))) 
			+ 2./3./n_ion * sources.electron_heat_source(t) - T_electron/n_ion*sources.particle_source(t);
  #else
    ddt(T_electron) = -Vpar_Grad_par(Vpar_ion,T_electron,CELL_CENTRE) 
			- 2./3.*T_electron*Grad_par(Vpar_ion,CELL_CENTRE) 
			- 2./3./n_ion*Grad_par(nonlocal_parallel.electron_heat_flux,CELL_CENTRE) 
			+ 2./3./n_ion * sources.electron_heat_source(t) - T_electron/n_ion*sources.particle_source(t);
  #endif
  ddt(n_ion) = -Vpar_Grad_par(Vpar_ion,n_ion,CELL_CENTRE) 
		- n_ion*Grad_par(Vpar_ion,CELL_CENTRE)
		+ sources.particle_source(t);

  #ifdef IONVISCOSITYLIMITER
    nTtau_ion_ylow = interp_to(nTtau_ion,CELL_YLOW);
    pi_Brag_ylow = -4./3.*0.96*nTtau_ion_ylow*Grad_par(Vpar_ion,CELL_YLOW);
    grad_par_viscosity = ion_viscosity_limiter / (( abs(pi_Brag_ylow)/interp_to(n_ion*T_ion,CELL_YLOW) + ion_viscosity_limiter )^2)
			  * ( -ion_viscosity_limiter*4./3.*0.96*(Grad_par(nTtau_ion,CELL_YLOW)*Grad_par(Vpar_ion,CELL_YLOW) + nTtau_ion_ylow*Grad2_par2(Vpar_ion,CELL_YLOW))
			      +pi_Brag_ylow*abs(pi_Brag_ylow)/interp_to((n_ion*T_ion)^2,CELL_YLOW)*Grad_par(n_ion*T_ion,CELL_YLOW) );
  #endif
  ddt(Vpar_ion) = -Vpar_Grad_par(Vpar_ion,Vpar_ion,CELL_YLOW)
		    + interp_to(1./(ion_mass*n_ion),CELL_YLOW) * (-Grad_par(n_ion*(T_ion+T_electron),CELL_YLOW) 
    #ifdef IONVISCOSITYLIMITER
								  - grad_par_viscosity)
    #else
								  + 4./3.*0.96*Grad_par(nTtau_ion,CELL_YLOW)*Grad_par(Vpar_ion,CELL_YLOW))
		    + 1./ion_mass * 4./3.*0.96*interp_to(T_ion * sqrt(2) * tau_ii,CELL_YLOW) * Grad2_par2(Vpar_ion,CELL_YLOW)
    #endif
		    - Vpar_ion*interp_to(sources.particle_source(t)/n_ion,CELL_YLOW);		// The sources of heat and particles do not add any momentum, so V decreases in proportion to the increase in n.
  
  #ifdef IONFLUXLIMITER
    grad_par_T_ion = Grad_par(T_ion,CELL_CENTRE);
    qSH_ion = -3.9/ion_mass*nTtau_ion*grad_par_T_ion;
    qFS_ion = sqrt(2./ion_mass)*n_ion*(T_ion^1.5);
    grad_par_ion_heat_flux = ion_flux_limiter / ((abs(qSH_ion)/qFS_ion+ion_flux_limiter)^2)
			      * (qSH_ion*abs(qSH_ion)/(qFS_ion^2)*sqrt(2./ion_mass)*Grad_par(n_ion*(T_ion^1.5),CELL_CENTRE)
				  - ion_flux_limiter*3.9/ion_mass*(Grad_par(nTtau_ion,CELL_CENTRE)*grad_par_T_ion + nTtau_ion*Grad2_par2(T_ion,CELL_CENTRE)));
  #endif
  #ifdef IONVISCOSITYLIMITER
    grad_par_V_centre = Grad_par(Vpar_ion,CELL_CENTRE);
    pi_Brag_centre = -4./3.*0.96*nTtau_ion*grad_par_V_centre;
    viscosity = ion_viscosity_limiter*pi_Brag_centre / ( abs(pi_Brag_centre)/n_ion/T_ion + ion_viscosity_limiter );
  #endif
  ddt(T_ion) = -Vpar_Grad_par(Vpar_ion,T_ion,CELL_CENTRE) 
		- 2./3.*T_ion*Grad_par(Vpar_ion,CELL_CENTRE) 
    #ifdef IONFLUXLIMITER
		- 2./3./n_ion * grad_par_ion_heat_flux
    #else
		+ 2./3./n_ion * 3.9/ion_mass*(Grad_par(nTtau_ion,CELL_CENTRE)*Grad_par(T_ion,CELL_CENTRE) 
						+ nTtau_ion * Grad2_par2(T_ion,CELL_CENTRE))
    #endif
    #ifdef IONVISCOSITYLIMITER
		- 2./3./n_ion * viscosity*grad_par_V_centre
    #else
		+ 2./3./n_ion * 4./3.*0.96*nTtau_ion*(Grad_par(Vpar_ion,CELL_CENTRE)^2)
    #endif
		+ 2./3./n_ion * sources.ion_heat_source(t) - T_ion/n_ion*sources.particle_source(t);
  
  ddt(j_parallel) = 0.;
  
  #ifdef CUSTOMBCS
    for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	for (int jy=mesh->ystart-1; jy>=0; jy--) {
	  ddt(n_ion)(rlow.ind,jy,jz) = 0.;
	  ddt(T_ion)(rlow.ind,jy,jz) = 0.;
	  ddt(T_electron)(rlow.ind,jy,jz) = 0.;
	}
      }
    if (mesh->StaggerGrids)
      for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
	for (int jz=0; jz<mesh->LocalNz; jz++) {
	  for (int jy=mesh->yend; jy<mesh->LocalNy; jy++) {
	    ddt(n_ion)(rup.ind,jy,jz) = 0.;
	    ddt(T_ion)(rup.ind,jy,jz) = 0.;
	    ddt(T_electron)(rup.ind,jy,jz) = 0.;
	  }
	}
    else
      for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
	for (int jz=0; jz<mesh->LocalNz; jz++) {
	  for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++) {
	    ddt(n_ion)(rup.ind,jy,jz) = 0.;
	    ddt(T_ion)(rup.ind,jy,jz) = 0.;
	    ddt(T_electron)(rup.ind,jy,jz) = 0.;
	  }
	}

    for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	for (int jy=mesh->ystart; jy>=0; jy--)
	  ddt(Vpar_ion)(rlow.ind,jy,jz) = 0.;
      }
    for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	for (int jy=mesh->yend; jy<mesh->LocalNy; jy++)
	  ddt(Vpar_ion)(rup.ind,jy,jz) = 0.;
      }
  #endif
  
  return 0;
}
