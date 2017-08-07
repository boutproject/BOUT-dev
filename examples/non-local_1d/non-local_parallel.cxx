/*!
 * \file heat_flux_integration.cxx
 *
 * \brief Calculate non-local electron closures
 *
 *
 **************************************************************************
 * Copyright 2012 J.T.Omotani
 *
 * Contact: John Omotani, john.omotani@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
  * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "non-local_parallel.hxx"

// #include <iomanip>

// #define NOEDGETERMS

BoutReal NonLocalParallelgamma_factor = 3.;

/**********************************************************************************
 *                    HEAT FLUX INITIALISATION AND CREATION                       *
 **********************************************************************************/


NonLocalParallel::~NonLocalParallel() {
  delete [] eigenvalues;
  #ifdef CALCULATE_HEATFLUX
    #ifdef DRIVE_GRADT
      delete [] heatflux_gradT_coefficients;
    #endif
    #ifdef DRIVE_GRADV
      delete [] heatflux_gradV_coefficients;
    #endif
    #ifdef DRIVE_VEMINUSVI
      delete [] heatflux_VeminusVi_coefficients;
    #endif
    delete [] heatflux_lower_boundary_transients;
    delete [] heatflux_upper_boundary_transients;
  #endif
  #ifdef BC_HEATFLUX
    delete [] heatflux_transients_factors;
    delete [] W11_B_times_WinverseB_11;
    #ifdef BC_VISCOSITY
      delete [] W11_B_times_WinverseB_20;
    #endif
  #endif
  #ifdef CALCULATE_VISCOSITY
    #ifdef DRIVE_GRADT
      delete [] viscosity_gradT_coefficients;
    #endif
    #ifdef DRIVE_GRADV
      delete [] viscosity_gradV_coefficients;
    #endif
    #ifdef DRIVE_VEMINUSVI
      delete [] viscosity_VeminusVi_coefficients;
    #endif
    delete [] viscosity_lower_boundary_transients;
    delete [] viscosity_upper_boundary_transients;
  #endif
  #ifdef BC_VISCOSITY
    delete [] viscosity_transients_factors;
    delete [] W20_B_times_WinverseB_20;
    #ifdef BC_HEATFLUX
      delete [] W20_B_times_WinverseB_11;
    #endif
  #endif
  #ifdef CALCULATE_FRICTION
    #ifdef DRIVE_GRADT
      delete [] friction_gradT_coefficients;
    #endif
    #ifdef DRIVE_GRADV
      delete [] friction_gradV_coefficients;
    #endif
    #ifdef DRIVE_VEMINUSVI
      delete [] friction_VeminusVi_coefficients;
    #endif
    delete [] friction_lower_boundary_transients;
    delete [] friction_upper_boundary_transients;
    #ifdef BC_HEATFLUX
      delete [] C10_1k_dot_W1k_B_times_WinverseB_11;
    #endif
    #ifdef BC_VISCOSITY
      delete [] C10_1k_dot_W1k_B_times_WinverseB_20;
    #endif
  #endif
  if (is_lower_boundary) {
    delete [] exp_total_dimensionless_length_over_eigenvalue;
  }
  MPI_Comm_free(&comm_yprocs_minusone);
}

void NonLocalParallel::initialise(BoutReal pass_electron_charge, BoutReal pass_electron_mass, BoutReal pass_ion_mass, BoutReal pass_epsilon_0, BoutReal pass_logLambda, const bool pass_fluxes_location_is_ylow, BoutReal pass_gamma_factor) {
  fluxes_location_is_ylow = pass_fluxes_location_is_ylow;
  #if CHECK > 0
    calculated_before_setting_bcs=false;
  #endif
  is_lower_boundary=mesh->firstY();
  is_upper_boundary=mesh->lastY();
  boundary_gradient_smoothing_length = mesh->GlobalNy/256; // Number of grid points to average over for the gradient used to extrapolate to the guard cells
  boundary_condition_smoothing_range = BoutReal(mesh->GlobalNy)/64.; // Scale length of the exponential decay used to smoothly apply the boundary condition
  if (boundary_gradient_smoothing_length==0)
    boundary_gradient_smoothing_length = 1;
  
  if (fluxes_location_is_ylow && !mesh->StaggerGrids)
    throw BoutException("Trying to calculate the heat flux at CELL_YLOW while StaggerGrids=false is an error.");
  
  cubic_spline_inverse_lambdaC.initialise('y',true,fluxes_location_is_ylow);
  if (fluxes_location_is_ylow) {
    #ifdef DRIVE_GRADT
      cubic_spline_gradT_driveterm.initialise('y',true,true); // CELL_CENTRE quantity, so must be adjusted when the grids are staggered.
      cubic_spline_gradT.initialise('y',true,false); // will fix up the boundary guard cells with forward/backward derivatives, since at least one guard cell value will be needed
    #endif
    #ifdef DRIVE_GRADV
      cubic_spline_gradV_driveterm.initialise('y',false,true); // CELL_CENTRE quantity, so must be adjusted when the grids are staggered.
      cubic_spline_one.initialise('y',true,false);
      Field3D one = 1.;
      cubic_spline_one.calculate(one);
    #endif
    #ifdef DRIVE_VEMINUSVI
      cubic_spline_VeminusVi_driveterm.initialise('y',true,true);
      cubic_spline_jpar.initialise('y',true,false);
    #endif
  }
  else {
    #ifdef DRIVE_GRADT
      cubic_spline_gradT_driveterm.initialise('y',false,false); // false because gradT_driveterm includes a derivative, so we don't know it in the guard cells, hence calculate the interpolation excluding the guard cells at the target boundaries.
    #endif
    #ifdef DRIVE_GRADV
      cubic_spline_gradV_driveterm.initialise('y',false,false); // false because gradT_driveterm includes a derivative, so we don't know it in the guard cells, hence calculate the interpolation excluding the guard cells at the target boundaries.
    #endif
    #ifdef DRIVE_VEMINUSVI
      cubic_spline_VeminusVi_driveterm.initialise('y',true,false);
    #endif
  }
  
  // Get the options for the model
  Options *options = Options::getRoot()->getSection("non_local_parallel");
  OPTION(options, moments_number, 10);
  OPTION(options, NONLOCAL_PARALLEL_TAGBASE, 12381);
  
  position = new bindex;
  broadcast_request = MPI_REQUEST_NULL;

  electron_charge = pass_electron_charge;
  electron_mass = pass_electron_mass;
  ion_mass = pass_ion_mass;
  epsilon_0 = pass_epsilon_0;
  logLambda = pass_logLambda;
  NonLocalParallelgamma_factor = pass_gamma_factor;
  
  if (fluxes_location_is_ylow) {
    increasing_dimensionless_length.setLocation(CELL_YLOW);
    decreasing_dimensionless_length.setLocation(CELL_YLOW);
    #ifdef CALCULATE_HEATFLUX
      electron_heat_flux.setLocation(CELL_YLOW);
    #endif
    #ifdef CALCULATE_VISCOSITY
      electron_viscosity.setLocation(CELL_YLOW); // should really have this as CELL_CENTRE quantity, but then would have to calculate the integrals on both CELL_CENTRE and CELL_YLOW which is perfectly possible but has not been implemented yet
    #endif
    #ifdef CALCULATE_FRICTION
      electron_friction.setLocation(CELL_YLOW);
    #endif
    #ifdef DRIVE_GRADT
      gradT_electron.setLocation(CELL_YLOW);
      gradT_electron=0.;
    #endif
    #ifdef DRIVE_GRADV
      // No CELL_YLOW drive terms here
    #endif
    #ifdef DRIVE_VEMINUSVI
      // VeminusVi already CELL_YLOW
    #endif
  }
  
  #ifdef CALCULATE_HEATFLUX
    electron_heat_flux = 0.;
  #endif
  #ifdef CALCULATE_VISCOSITY
    electron_viscosity = 0.;
  #endif
  #ifdef CALCULATE_FRICTION
    electron_friction = 0.;
  #endif
  #ifdef DRIVE_GRADT
    gradT_driveterm = 0.;
  #endif
  #ifdef DRIVE_GRADV
    gradV_driveterm = 0.;
  #endif
  #ifdef DRIVE_VEMINUSVI
    VeminusVi_driveterm = 0.;
  #endif
  lambdaC_inverse = 0.;
  increasing_dimensionless_length = 0.;
  decreasing_dimensionless_length = 0.;
  dimensionless_length_deltas_above = 0.;
  if (fluxes_location_is_ylow)
    dimensionless_length_deltas_below=0.;
  total_dimensionless_length = 1./0.;
  
  integration.initialise(fluxes_location_is_ylow);
  
  // Get model's eigenvalues and coefficients from files
  BoutReal junk;
  #ifdef CALCULATE_HEATFLUX
    std::stringstream heatflux_infilename;
    heatflux_infilename<<"nonlocal_coefficients/heatfluxcoeffs"<<moments_number;
    std::ifstream heatflux_infile ( heatflux_infilename.str().c_str() );
    if (!heatflux_infile.is_open())
      throw BoutException("Could not open heatfluxcoeffs file");
    heatflux_infile>>number_of_negative_eigenvalues;
    #ifndef ALLOCATED_EIGENVALUES
      #define ALLOCATED_EIGENVALUES
      eigenvalues = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_GRADT
      heatflux_gradT_coefficients = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_GRADV
      heatflux_gradV_coefficients = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_VEMINUSVI
      heatflux_VeminusVi_coefficients = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_GRADT
      heatflux_infile>>heatflux_gradT_zerocoeff;
    #else
      heatflux_infile>>junk;
    #endif
    #ifdef DRIVE_GRADV
      heatflux_infile>>junk; // the zero eigenvector has only odd-l components, so the gradV drive (which is (2,0)) does not contribute
    #else
      heatflux_infile>>junk;
    #endif
    #ifdef DRIVE_VEMINUSVI
      heatflux_infile>>heatflux_VeminusVi_zerocoeff;
    #else
      heatflux_infile>>junk;
    #endif
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (heatflux_infile.eof())
	throw BoutException("reached end of heatfluxcoeffs file unexpectedly");
      heatflux_infile>>eigenvalues[i];
      #ifdef DRIVE_GRADT
	heatflux_infile>>heatflux_gradT_coefficients[i];
      #else
	heatflux_infile>>junk;
      #endif
      #ifdef DRIVE_GRADV
	heatflux_infile>>heatflux_gradV_coefficients[i];
      #else
	heatflux_infile>>junk;
      #endif
      #ifdef DRIVE_VEMINUSVI
	heatflux_infile>>heatflux_VeminusVi_coefficients[i];
      #else
	heatflux_infile>>junk;
      #endif
    }
    heatflux_infile.close();
    std::stringstream heatfluxbc_infilename;
    heatfluxbc_infilename<<"nonlocal_coefficients/heatfluxbc"<<moments_number;
    std::ifstream heatfluxbc_infile ( heatfluxbc_infilename.str().c_str() );
    if (!heatfluxbc_infile.is_open())
      throw BoutException("Could not open heatfluxbc file");
    #ifdef BC_HEATFLUX
      W11_B_times_WinverseB_11 = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef BC_VISCOSITY
      W11_B_times_WinverseB_20 = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef BC_HEATFLUX
      heatfluxbc_infile>>W11_dot_W11;
    #else
      heatfluxbc_infile>>junk;
    #endif
    #ifdef BC_VISCOSITY
      heatfluxbc_infile>>W11_dot_W20;
    #else
      heatfluxbc_infile>>junk;
    #endif
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (heatfluxbc_infile.eof()) {
	output<<"Error at i="<<i<<endl;
	throw BoutException("reached end of heatfluxbc file unexpectedly");
      }
      #ifdef BC_HEATFLUX
	heatfluxbc_infile>>W11_B_times_WinverseB_11[i];
      #else
	heatfluxbc_infile>>junk;
      #endif
      #ifdef BC_VISCOSITY
	heatfluxbc_infile>>W11_B_times_WinverseB_20[i];
      #else
	heatfluxbc_infile>>junk;
      #endif
    }
    heatfluxbc_infile.close();
  #endif
  #ifdef BC_HEATFLUX
    if (is_lower_boundary || is_upper_boundary) {
      pass_interim_upper_boundary_n11 = 0.;
      upper_boundary_condition_n11 = 0.;
    }
  #endif
  #ifdef CALCULATE_VISCOSITY
    std::stringstream viscosity_infilename;
    viscosity_infilename<<"nonlocal_coefficients/viscositycoeffs"<<moments_number;
    std::ifstream viscosity_infile ( viscosity_infilename.str().c_str() );
    if (!viscosity_infile.is_open())
      throw BoutException("Could not open viscositycoeffs file");
    viscosity_infile>>number_of_negative_eigenvalues;
    #ifndef ALLOCATED_EIGENVALUES
      #define ALLOCATED_EIGENVALUES
      eigenvalues = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_GRADT
      viscosity_gradT_coefficients = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_GRADV
      viscosity_gradV_coefficients = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_VEMINUSVI
      viscosity_VeminusVi_coefficients = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_GRADT
      viscosity_infile>>junk; // the zero eigenvector has only odd-l components, so does not contribute to the viscosity (2,0)
    #else
      viscosity_infile>>junk;
    #endif
    #ifdef DRIVE_GRADV
      viscosity_infile>>junk; // the zero eigenvector has only odd-l components, so does not contribute to the viscosity (2,0)
    #else
      viscosity_infile>>junk;
    #endif
    #ifdef DRIVE_GRADV
      viscosity_infile>>junk; // the zero eigenvector has only odd-l components, so does not contribute to the viscosity (2,0)
    #else
      viscosity_infile>>junk;
    #endif
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (viscosity_infile.eof())
	throw BoutException("reached end of viscositycoeffs file unexpectedly");
      viscosity_infile>>eigenvalues[i];
      #ifdef DRIVE_GRADT
	viscosity_infile>>viscosity_gradT_coefficients[i];
      #else
	viscosity_infile>>junk;
      #endif
      #ifdef DRIVE_GRADV
	viscosity_infile>>viscosity_gradV_coefficients[i];
      #else
	viscosity_infile>>junk;
      #endif
      #ifdef DRIVE_VEMINUSVI
	viscosity_infile>>viscosity_VeminusVi_coefficients[i];
      #else
	viscosity_infile>>junk;
      #endif
    }
    viscosity_infile.close();
    std::stringstream viscositybc_infilename;
    viscositybc_infilename<<"nonlocal_coefficients/viscositybc"<<moments_number;
    std::ifstream viscositybc_infile ( viscositybc_infilename.str().c_str() );
    if (!viscositybc_infile.is_open())
      throw BoutException("Could not open viscositybc file");
    #ifdef BC_VISCOSITY
      W20_B_times_WinverseB_20 = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef BC_HEATFLUX
      W20_B_times_WinverseB_11 = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef BC_VISCOSITY
      viscositybc_infile>>W20_dot_W20;
    #else
      viscositybc_infile>>junk;
    #endif
    #ifdef BC_HEATFLUX
      viscositybc_infile>>W20_dot_W11;
    #else
      viscositybc_infile>>junk;
    #endif
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (viscositybc_infile.eof())
	throw BoutException("reached end of viscositybc file unexpectedly");
      #ifdef BC_VISCOSITY
	viscositybc_infile>>W20_B_times_WinverseB_20[i];
      #else
	viscositybc_infile>>junk;
      #endif
      #ifdef BC_HEATFLUX
	viscositybc_infile>>W20_B_times_WinverseB_11[i];
      #else
	viscositybc_infile>>junk;
      #endif
    }
    viscositybc_infile.close();
  #endif
  #ifdef BC_VISCOSITY
    if (is_lower_boundary || is_upper_boundary) {
      pass_interim_upper_boundary_n20 = 0.;
      upper_boundary_condition_n20 = 0.;
    }
  #endif
  #ifdef CALCULATE_FRICTION
    std::stringstream friction_infilename;
    friction_infilename<<"nonlocal_coefficients/frictioncoeffs"<<moments_number;
    std::ifstream friction_infile ( friction_infilename.str().c_str() );
    if (!friction_infile.is_open())
      throw BoutException("Could not open frictioncoeffs file");
    friction_infile>>number_of_negative_eigenvalues;
    #ifndef ALLOCATED_EIGENVALUES
      #define ALLOCATED_EIGENVALUES
      eigenvalues = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_GRADT
      friction_gradT_coefficients = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_GRADV
      friction_gradV_coefficients = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_VEMINUSVI
      friction_VeminusVi_coefficients = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef DRIVE_GRADT
      friction_infile>>friction_gradT_zerocoeff;
    #else
      friction_infile>>junk;
    #endif
    #ifdef DRIVE_GRADV
      friction_infile>>junk; // the zero eigenvector has only odd-l components, so the gradV drive (which is (2,0)) does not contribute
    #else
      friction_infile>>junk;
    #endif
    #ifdef DRIVE_VEMINUSVI
      friction_infile>>friction_VeminusVi_zerocoeff;
    #else
      friction_infile>>junk;
    #endif
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (friction_infile.eof())
	throw BoutException("reached end of frictioncoeffs file unexpectedly");
      friction_infile>>eigenvalues[i];
      #ifdef DRIVE_GRADT
	friction_infile>>friction_gradT_coefficients[i];
      #else
	friction_infile>>junk;
      #endif
      #ifdef DRIVE_GRADV
	friction_infile>>friction_gradV_coefficients[i];
      #else
	friction_infile>>junk;
      #endif
      #ifdef DRIVE_VEMINUSVI
	friction_infile>>friction_VeminusVi_coefficients[i];
      #else
	friction_infile>>junk;
      #endif
    }
    friction_infile.close();
    std::stringstream frictionbc_infilename;
    frictionbc_infilename<<"nonlocal_coefficients/frictionbc"<<moments_number;
    std::ifstream frictionbc_infile ( frictionbc_infilename.str().c_str() );
    if (!frictionbc_infile.is_open())
      throw BoutException("Could not open frictionbc file");
    #ifdef BC_HEATFLUX
      C10_1k_dot_W1k_B_times_WinverseB_11 = new BoutReal[number_of_negative_eigenvalues];
    #endif
    #ifdef BC_VISCOSITY
      C10_1k_dot_W1k_B_times_WinverseB_20 = new BoutReal[number_of_negative_eigenvalues];
    #endif
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (frictionbc_infile.eof())
	throw BoutException("reached end of frictionbc file unexpectedly");
      #ifdef BC_HEATFLUX
      frictionbc_infile>>C10_1k_dot_W1k_B_times_WinverseB_11[i];
      #else
	frictionbc_infile>>junk;
      #endif
      #ifdef BC_VISCOSITY
	frictionbc_infile>>C10_1k_dot_W1k_B_times_WinverseB_20[i];
      #else
	frictionbc_infile>>junk;
      #endif
    }
    frictionbc_infile.close();
  #endif
  
  if (is_lower_boundary) {
    exp_total_dimensionless_length_over_eigenvalue = new BoutReal[number_of_negative_eigenvalues];
  }
      
  #ifdef CALCULATE_HEATFLUX
    heatflux_lower_boundary_transients = new FieldPerp[number_of_negative_eigenvalues];
    heatflux_upper_boundary_transients = new FieldPerp[number_of_negative_eigenvalues];
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      heatflux_lower_boundary_transients[i]= 0.;
      heatflux_upper_boundary_transients[i]= 0.;
    }
  #endif
  #ifdef BC_HEATFLUX
    heatflux_transients_factors = new BoutReal[2*(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)];
  #endif
  #ifdef CALCULATE_VISCOSITY
    viscosity_lower_boundary_transients = new FieldPerp[number_of_negative_eigenvalues];
    viscosity_upper_boundary_transients = new FieldPerp[number_of_negative_eigenvalues];
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      viscosity_lower_boundary_transients[i]= 0.;
      viscosity_upper_boundary_transients[i]= 0.;
    }
  #endif
  #ifdef BC_VISCOSITY
    viscosity_transients_factors = new BoutReal[2*(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)];
  #endif
  #ifdef CALCULATE_FRICTION
    friction_lower_boundary_transients = new FieldPerp[number_of_negative_eigenvalues];
    friction_upper_boundary_transients = new FieldPerp[number_of_negative_eigenvalues];
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      friction_lower_boundary_transients[i]= 0.;
      friction_upper_boundary_transients[i]= 0.;
    }
  #endif

  // Initialisation for stuff used in y_broadcast functions
  {
    MPI_Group group_yprocs;
    
    int n_yprocs = mesh->getNYPE()-1;
    int * indices_yprocs = new int[n_yprocs];
    for (int i=0; i<n_yprocs; i++)
      indices_yprocs[i] = (i+1) * mesh->getNXPE() + mesh->getXProcIndex();
    
    MPI_Group group_world;
    MPI_Comm_group(BoutComm::get(), &group_world); // Get the entire group
    MPI_Group_incl(group_world, n_yprocs, indices_yprocs, &group_yprocs);
    MPI_Group_free(&group_world);
    delete [] indices_yprocs;
    
    MPI_Comm_create(BoutComm::get(), group_yprocs, &comm_yprocs_minusone);
    MPI_Group_free(&group_yprocs);
  }
  
}

/**********************************************************************************
 *                    HEAT FLUX CALCULATION ROUTINES
 **********************************************************************************/

void NonLocalParallel::calculate_nonlocal_closures(const Field3D &n_electron, const Field3D &T_electron
						    #ifdef DRIVE_GRADV
						      , const Field3D &V_electron
						    #endif
						    #ifdef DRIVE_VEMINUSVI
						      , const Field3D &jpar
						    #endif
						    #ifdef BC_HEATFLUX
						      , const Field3D &heat_flux_boundary_condition
						    #endif
						    #ifdef BC_VISCOSITY
						      , const Field3D &viscosity_boundary_condition
						    #endif
						    ) {
  if (fluxes_location_is_ylow)
    calculate_nonlocal_closures_cell_ylow(n_electron, T_electron
					  #ifdef DRIVE_GRADV
					    , V_electron
					  #endif
					  #ifdef DRIVE_VEMINUSVI
					    , jpar
					  #endif
					  #ifdef BC_HEATFLUX
					    , heat_flux_boundary_condition
					  #endif
					  #ifdef BC_VISCOSITY
					    , viscosity_boundary_condition
					  #endif
					  );
  else
    calculate_nonlocal_closures_cell_centre(n_electron, T_electron
					    #ifdef DRIVE_GRADV
					      , V_electron
					    #endif
					    #ifdef DRIVE_VEMINUSVI
					      , jpar
					    #endif
					    #ifdef BC_HEATFLUX
					      , heat_flux_boundary_condition
					    #endif
					    #ifdef BC_VISCOSITY
					      , viscosity_boundary_condition
					    #endif
					  );
  #if CHECK > 0
    calculated_before_setting_bcs=true;
  #endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NonLocalParallel::calculate_nonlocal_closures_cell_centre(const Field3D &n_electron, const Field3D &T_electron
								#ifdef DRIVE_GRADV
								  , const Field3D &V_electron
								#endif
								#ifdef DRIVE_VEMINUSVI
								  , const Field3D &jpar
								#endif
								#ifdef BC_HEATFLUX
								  , const Field3D &heat_flux_boundary_condition
								#endif
								#ifdef BC_VISCOSITY
								  , const Field3D &viscosity_boundary_condition
								#endif
								) {
  
  lambdaC_inverse = n_electron * pow(electron_charge,4) * logLambda / 12 / pow(PI,1.5) / pow(epsilon_0,2) / (T_electron^2);
  
  #ifdef DRIVE_GRADT
    gradT_driveterm = 5./4. * n_electron / T_electron * Grad_par(T_electron) / lambdaC_inverse; //g^(1,1)
    mesh->communicate(gradT_driveterm);
  #endif
  #ifdef DRIVE_GRADV
    gradV_driveterm = -0.5 * n_electron / sqrt(2.*T_electron/electron_mass) * Grad_par(V_electron) / lambdaC_inverse;
    mesh->communicate(gradV_driveterm);
  #endif
  #ifdef DRIVE_VEMINUSVI
    VeminusVi_driveterm = -2./sqrt(PI) * (-jpar/electron_charge) / sqrt(2.*T_electron/electron_mass);
    mesh->communicate(VeminusVi_driveterm);
  #endif
  
  // Now calculate z and deltaz everywhere
  cubic_spline_inverse_lambdaC.calculate(lambdaC_inverse);
  
  start_index(position);
  
  if (mesh->UpXSplitIndex()!=0 || mesh->DownXSplitIndex()!=0) throw BoutException("This code cannot currently handle x-splitting of processors.");
    if (!is_lower_boundary) {
      FieldPerp pass_dimensionless_length;
      pass_dimensionless_length.allocate();
      {
	mesh->wait(mesh->irecvYInOutdest(*pass_dimensionless_length.getData(),mesh->LocalNx*(mesh->LocalNz),
				       NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz));
      }
      pass_dimensionless_length.setIndex(mesh->ystart);
      increasing_dimensionless_length = pass_dimensionless_length;
    }
    
  do {
    position->jy = mesh->ystart-1;
    calc_index(position);
    interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
    // d/dy(delta) = 1/lambdaC = a + b*t + c*t^2 + d*t^3; t=(ind-jy)=(y-y0)/(sqrt(g_22)*dy); ind is a notional continuous variable equal to jy at the gridpoints so at jy+1 t=1

    // Fetch coordinate system
    Coordinates *coord = mesh->coordinates();
    dimensionless_length_deltas_above[*position] /* = dy/dt*(a + 1/2*b + 1/3*c + 1/4*d) */
      = coord->dy(position->jx,position->jy)*sqrt(0.5*(coord->g_22(position->jx,position->jy)+coord->g_22(position->jx,position->jyp)))
      *(interp_coefficients[0] + interp_coefficients[1]/2. + interp_coefficients[2]/3. + interp_coefficients[3]/4.);
    next_index_y(position);
    
    do{
      interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
      // d/dy(delta) = 1/lambdaC = a + b*t + c*t^2 + d*t^3; t=(ind-jy)=(y-y0)/(sqrt(g_22)*dy); ind is a notional continuous variable equal to jy at the gridpoints so at jy+1 t=1
      dimensionless_length_deltas_above[*position] /* = dy/dt*(a + 1/2*b + 1/3*c + 1/4*d) */
	= coord->dy(position->jx,position->jy)*sqrt(0.5*(coord->g_22(position->jx,position->jy)+coord->g_22(position->jx,position->jyp)))
	*(interp_coefficients[0] + interp_coefficients[1]/2. + interp_coefficients[2]/3. + interp_coefficients[3]/4.);
      increasing_dimensionless_length(position->jx,position->jyp,position->jz) = increasing_dimensionless_length[*position] + dimensionless_length_deltas_above[*position];
    } while (next_index_y(position));
    
  } while (next_indexperp(position));
  
  {
    Timer timer("comms");
    mesh->sendYOutOutdest(*increasing_dimensionless_length.slice(mesh->yend+1).getData(),mesh->LocalNx*(mesh->LocalNz),
			  NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz);
  }
  
  // Send the total dimensionless_length at the upper boundary back to the other processors.
  if (is_upper_boundary)
    total_dimensionless_length = increasing_dimensionless_length.slice(mesh->yend);

  y_broadcast(*total_dimensionless_length.getData(), mesh->LocalNx*mesh->LocalNz, mesh->getNYPE()-1);
  
  decreasing_dimensionless_length = -increasing_dimensionless_length;
  for (int jy=mesh->ystart; jy<=mesh->yend; jy++) {
    total_dimensionless_length.setIndex(jy);
    decreasing_dimensionless_length += total_dimensionless_length;
  }
  
  #ifdef CALCULATE_HEATFLUX
    electron_heat_flux = 0.;
  #endif
  #ifdef CALCULATE_VISCOSITY
    electron_viscosity = 0.;
  #endif
  #ifdef CALCULATE_FRICTION
    electron_friction = 0.;
  #endif
  #ifdef DRIVE_GRADT
    #ifdef CALCULATE_HEATFLUX
      electron_heat_flux += -heatflux_gradT_zerocoeff * gradT_driveterm; //zero eigenvalue contribution to n^(1,1) //zero eigenvalue contribution to n^(1,1)/T^1.
    #endif
    // viscosity gets no zero eigenvalue contribution
    #ifdef CALCULATE_FRICTION
      electron_friction += -friction_gradT_zerocoeff * gradT_driveterm;
    #endif
    cubic_spline_gradT_driveterm.calculate(gradT_driveterm);
  #endif
  #ifdef DRIVE_GRADV
    // gradV drive has no zero eigenvalue contribution
    cubic_spline_gradV_driveterm.calculate(gradV_driveterm);
  #endif
  #ifdef DRIVE_VEMINUSVI
    #ifdef CALCULATE_HEATFLUX
      electron_heat_flux += -heatflux_VeminusVi_zerocoeff * VeminusVi_driveterm;
    #endif
    // viscosity gets no zero eigenvalue contribution
    #ifdef CALCULATE_FRICTION
      electron_friction += -friction_VeminusVi_zerocoeff * VeminusVi_driveterm;
    #endif
    cubic_spline_VeminusVi_driveterm.calculate(VeminusVi_driveterm);
  #endif

  for (int i=0; i<number_of_negative_eigenvalues; i++) {
    #ifdef DRIVE_GRADT
      integration.calculateIntegralBelow_cell_centre(eigenvalues[i], dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_gradT_driveterm, i+0);
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += heatflux_gradT_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += viscosity_gradT_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += friction_gradT_coefficients[i]*integration.integral_below;
      #endif
    #endif
    #ifdef DRIVE_GRADV
      integration.calculateIntegralBelow_cell_centre(eigenvalues[i], dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_gradV_driveterm, i+number_of_negative_eigenvalues);
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += heatflux_gradV_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += viscosity_gradV_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_FRICTION
      electron_friction += friction_gradV_coefficients[i]*integration.integral_below;
      #endif
    #endif
    #ifdef DRIVE_VEMINUSVI
      integration.calculateIntegralBelow_cell_centre(eigenvalues[i], dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_VeminusVi_driveterm, i+2*number_of_negative_eigenvalues);
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += heatflux_VeminusVi_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += viscosity_VeminusVi_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += friction_VeminusVi_coefficients[i]*integration.integral_below;
      #endif
    #endif
  }
    
  for (int i=0; i<number_of_negative_eigenvalues; i++) {
    #ifdef DRIVE_GRADT
      integration.calculateIntegralAbove_cell_centre(eigenvalues[i], dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_gradT_driveterm, i+0);
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += -heatflux_gradT_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += -viscosity_gradT_coefficients[i]*integration.integral_above;
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += -friction_gradT_coefficients[i]*integration.integral_above;
      #endif
    #endif
    #ifdef DRIVE_GRADV
      integration.calculateIntegralAbove_cell_centre(eigenvalues[i], dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_gradV_driveterm, i+number_of_negative_eigenvalues);
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += -heatflux_gradV_coefficients[i]*integration.integral_above; //gives \hat(n)^A
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += -viscosity_gradV_coefficients[i]*integration.integral_above;
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += -friction_gradV_coefficients[i]*integration.integral_above;
      #endif
    #endif
    #ifdef DRIVE_VEMINUSVI
      integration.calculateIntegralAbove_cell_centre(eigenvalues[i], dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_VeminusVi_driveterm, i+2*number_of_negative_eigenvalues);
      
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += -heatflux_VeminusVi_coefficients[i]*integration.integral_above; //gives \hat(n)^A
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += -viscosity_VeminusVi_coefficients[i]*integration.integral_above;
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += -friction_VeminusVi_coefficients[i]*integration.integral_above;
      #endif
    #endif
  }
  
//   throw BoutException("cell_centre version of boundary conditions code needs checking, at the moment it has just been cut&pasted from cell_ylow");
  for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
    for (int jz=0; jz<mesh->LocalNz; jz++) {
      position->jx=rup.ind;
      position->jy=mesh->yend;
      position->jz=jz;
      calc_index(position);
      BoutReal Te_here = interp_to_point_YLOW(T_electron,*position);
      #ifdef BC_HEATFLUX
	pass_interim_upper_boundary_n11[rup.ind][jz] = electron_heat_flux[rup.ind][mesh->yend][jz];
	upper_boundary_condition_n11[rup.ind][jz] = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition[rup.ind][mesh->yend][jz];
      #endif
      #ifdef BC_VISCOSITY
	pass_interim_upper_boundary_n20[rup.ind][jz] = electron_viscosity[rup.ind][mesh->yend][jz];
	upper_boundary_condition_n20[rup.ind][jz] = viscosity_boundary_condition[rup.ind][mesh->yend][jz]/Te_here;
      #endif
    }
  if (is_upper_boundary && mesh->getNYPE()>1) {
    #ifdef BC_HEATFLUX
      MPI_Request request1 = mesh->sendToProc(mesh->getXProcIndex(),0,
					      *pass_interim_upper_boundary_n11.getData(),
					      mesh->LocalNx*mesh->LocalNz,
					      NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex());
      MPI_Request request2 = mesh->sendToProc(mesh->getXProcIndex(),0,
					      *upper_boundary_condition_n11.getData(),
					      mesh->LocalNx*mesh->LocalNz,
					      NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 1);
      MPI_Waitall(1,&request1,MPI_STATUSES_IGNORE);
      MPI_Waitall(1,&request2,MPI_STATUSES_IGNORE);
    #endif
    #ifdef BC_VISCOSITY
      MPI_Request request3 = mesh->sendToProc(mesh->getXProcIndex(),0,
					      *pass_interim_upper_boundary_n20.getData(),
					      mesh->LocalNx*mesh->LocalNz,
					      NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 3);
      MPI_Request request4 = mesh->sendToProc(mesh->getXProcIndex(),0,
					      *upper_boundary_condition_n20.getData(),
					      mesh->LocalNx*mesh->LocalNz,
					      NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 4);
      MPI_Waitall(1,&request3,MPI_STATUSES_IGNORE);
      MPI_Waitall(1,&request4,MPI_STATUSES_IGNORE);
    #endif
  }
  if (is_lower_boundary && mesh->getNYPE()>1) {
    #ifdef BC_HEATFLUX
      mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(), mesh->getNYPE()-1,
					    *pass_interim_upper_boundary_n11.getData(),
					    mesh->LocalNx*mesh->LocalNz,
					    NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex()) );
      mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(), mesh->getNYPE()-1,
					    *upper_boundary_condition_n11.getData(),
					    mesh->LocalNx*mesh->LocalNz,
					    NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 1) );
    #endif
    #ifdef BC_VISCOSITY
      mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(), mesh->getNYPE()-1,
					    *pass_interim_upper_boundary_n20.getData(),
					    mesh->LocalNx*mesh->LocalNz,
					    NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 3) );
      mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(), mesh->getNYPE()-1,
					    *upper_boundary_condition_n20.getData(),
					    mesh->LocalNx*mesh->LocalNz,
					    NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 4) );
    #endif
  }
  
  if (is_lower_boundary) {
    for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	position->jx=jx;
	position->jy=mesh->ystart;
	position->jz=jz;
	calc_index(position);
	BoutReal Te_here = interp_to_point_YLOW(T_electron,*position);
	#if defined(BC_HEATFLUX) && !defined(BC_VISCOSITY)
	  BoutReal lower_boundary_n11 = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition[jx][mesh->ystart][jz];
	  BoutReal upper_boundary_n11 = upper_boundary_condition_n11[jx][jz];
	  BoutReal interim_lower_boundary_n11 = electron_heat_flux[jx][mesh->ystart][jz];
	  BoutReal interim_upper_boundary_n11 = pass_interim_upper_boundary_n11[jx][jz];
	  /*
	  electron_heat_flux is, at this point, the contribution to n11 from nhat_plus.
	  We want the actual heat flux at mesh->ystart to be boundary_heat_flux.
	  Thus the remainder must come from nhat_minus, which we will construct here just to give the right 1,1 component (could set number_of_negative_eigenvalues-1 more components if desired)
	  However, the transients from the other boundary do not necessarily decay to vanishing by the time they get to this boundary, so we must solve for both.
	  Fortunately this reduces to a single algebraic equation which makes it surprisingly easy to do (in the case of just a single condition being imposed at least).
	  */
	  
	  BoutReal sum_decayed_W11_W11_term = 0.;
	  for (int i=0; i<number_of_negative_eigenvalues; i++) {
	    exp_total_dimensionless_length_over_eigenvalue[i] = exp(total_dimensionless_length[jx][jz]/eigenvalues[i]);
	    sum_decayed_W11_W11_term += W11_B_times_WinverseB_11[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	  }
	  heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (lower_boundary_n11 - interim_lower_boundary_n11)*W11_dot_W11
								    - sum_decayed_W11_W11_term*(upper_boundary_n11 - interim_upper_boundary_n11) )
								    / ( pow(W11_dot_W11,2) - pow(sum_decayed_W11_W11_term,2) );
	  heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11 - interim_upper_boundary_n11)*W11_dot_W11
													    - sum_decayed_W11_W11_term*(lower_boundary_n11 - interim_lower_boundary_n11) )
													    / ( pow(W11_dot_W11,2) - pow(sum_decayed_W11_W11_term,2) );
	#elif defined(BC_VISCOSITY) && !defined(BC_HEATFLUX)
	  BoutReal lower_boundary_n20 = viscosity_boundary_condition[jx][mesh->ystart][jz]/Te_here;
	  BoutReal upper_boundary_n20 = upper_boundary_condition_n20[jx][jz];
	  BoutReal interim_lower_boundary_n20 = electron_viscosity[jx][mesh->ystart][jz];
	  BoutReal interim_upper_boundary_n20 = pass_interim_upper_boundary_n20[jx][jz];
	  /*
	  electron_viscosity is, at this point, the contribution to n20 from nhat_plus.
	  We want the actual viscosity at mesh->ystart to be boundary_viscosity.
	  Thus the remainder must come from nhat_minus, which we will construct here just to give the right 2,0 component (could set number_of_negative_eigenvalues-1 more components if desired)
	  However, the transients from the other boundary do not necessarily decay to vanishing by the time they get to this boundary, so we must solve for both.
	  Fortunately this reduces to a single algebraic equation which makes it surprisingly easy to do (in the case of just a single condition being imposed at least).
	  */
	  
	  BoutReal sum_decayed_W20_W20_term = 0.;
	  for (int i=0; i<number_of_negative_eigenvalues; i++) {
	    exp_total_dimensionless_length_over_eigenvalue[i] = exp(total_dimensionless_length[jx][jz]/eigenvalues[i]);
	    sum_decayed_W20_W20_term += W20_B_times_WinverseB_20[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	  }
	  viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (lower_boundary_n20 - interim_lower_boundary_n20)*W20_dot_W20
								    - sum_decayed_W20_W20_term*(upper_boundary_n20 - interim_upper_boundary_n20) )
								    / ( pow(W20_dot_W20,2) - pow(sum_decayed_W20_W20_term,2) );
	  viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n20 - interim_upper_boundary_n20)*W20_dot_W20
													    - sum_decayed_W20_W20_term*(lower_boundary_n20 - interim_lower_boundary_n20) )
													    / ( pow(W20_dot_W20,2) - pow(sum_decayed_W20_W20_term,2));
	#elif defined(BC_HEATFLUX) && defined(BC_VISCOSITY)
	  BoutReal lower_boundary_n11 = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition[jx][mesh->ystart][jz];
	  BoutReal upper_boundary_n11 = upper_boundary_condition_n11[jx][jz];
	  BoutReal interim_lower_boundary_n11 = electron_heat_flux[jx][mesh->ystart][jz];
	  BoutReal interim_upper_boundary_n11 = pass_interim_upper_boundary_n11[jx][jz];
	  BoutReal lower_boundary_n20 = viscosity_boundary_condition[jx][mesh->ystart][jz]/Te_here;
	  BoutReal upper_boundary_n20 = upper_boundary_condition_n20[jx][jz];
	  BoutReal interim_lower_boundary_n20 = electron_viscosity[jx][mesh->ystart][jz];
	  BoutReal interim_upper_boundary_n20 = pass_interim_upper_boundary_n20[jx][jz];
	  BoutReal sum_decayed_W11_W11_term = 0.;
	  BoutReal sum_decayed_W20_W20_term = 0.;
	  BoutReal sum_decayed_W11_W20_term = 0.;
	  BoutReal sum_decayed_W20_W11_term = 0.;
	  for (int i=0; i<number_of_negative_eigenvalues; i++) {
	    exp_total_dimensionless_length_over_eigenvalue[i] = exp(total_dimensionless_length[jx][jz]/eigenvalues[i]);
	    sum_decayed_W11_W11_term += W11_B_times_WinverseB_11[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	    sum_decayed_W20_W20_term += W20_B_times_WinverseB_20[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	    sum_decayed_W20_W11_term += W20_B_times_WinverseB_11[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	    sum_decayed_W11_W20_term += W11_B_times_WinverseB_20[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	  }
	  BoutReal det = pow(sum_decayed_W11_W20_term,2)*pow(sum_decayed_W20_W11_term,2)
			 - 2*sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
			 + pow(sum_decayed_W11_W11_term,2)*pow(sum_decayed_W20_W20_term,2)
			 - pow(sum_decayed_W20_W20_term,2)*pow(W11_dot_W11,2)
			 + 2*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
			 - pow(sum_decayed_W20_W11_term,2)*pow(W11_dot_W20,2)
			 - 2*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
			 + 2*sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
			 - pow(sum_decayed_W11_W20_term,2)*pow(W20_dot_W11,2)
			 + pow(W11_dot_W20,2)*pow(W20_dot_W11,2)
			 + 2*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
			 - 2*sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
			 + 2*sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
			 - 2*W11_dot_W11*W11_dot_W20*W20_dot_W11*W20_dot_W20
			 - pow(sum_decayed_W11_W11_term,2)*pow(W20_dot_W20,2)
			 + pow(W11_dot_W11,2)*pow(W20_dot_W20,2);
	  heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																+ sum_decayed_W11_W11_term*pow(sum_decayed_W20_W20_term,2)
																+ sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
																- sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
																+ sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
																- sum_decayed_W11_W11_term*pow(W20_dot_W20,2) )
										+ (upper_boundary_n20-interim_upper_boundary_n20)*( pow(sum_decayed_W11_W20_term,2)*sum_decayed_W20_W11_term
																    - sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term
																    + sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
																    - sum_decayed_W20_W11_term*pow(W11_dot_W20,2)
																    + sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W20
																    - sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W20 )
										+ (lower_boundary_n11-interim_lower_boundary_n11)*( -pow(sum_decayed_W20_W20_term,2)*W11_dot_W11
																    + sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																    - sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W20_dot_W11
																    + sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W20_dot_W20
																    - W11_dot_W20*W20_dot_W11*W20_dot_W20 + W11_dot_W11*pow(W20_dot_W20,2) )
										+ (lower_boundary_n20-interim_lower_boundary_n20)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11
																    + sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																    - pow(sum_decayed_W11_W20_term,2)*W20_dot_W11
																    + pow(W11_dot_W20,2)*W20_dot_W11
																    + sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W20
																    - W11_dot_W11*W11_dot_W20*W20_dot_W20 )
									    ) / det;
	  heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -pow(sum_decayed_W20_W20_term,2)*W11_dot_W11
																					  + sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																					  - sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W20_dot_W11
																					  + sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W20_dot_W20
																					  - W11_dot_W20*W20_dot_W11*W20_dot_W20
																					  + W11_dot_W11*pow(W20_dot_W20,2) )
															+ (upper_boundary_n20-interim_upper_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11
																					    - sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																					    + pow(sum_decayed_W11_W20_term,2)*W20_dot_W11
																					    - pow(W11_dot_W20,2)*W20_dot_W11
																					    - sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W20
																					    + W11_dot_W11*W11_dot_W20*W20_dot_W20 )
															+ (lower_boundary_n11-interim_lower_boundary_n11)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																					    + sum_decayed_W11_W11_term*pow(sum_decayed_W20_W20_term,2)
																					    + sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
																					    - sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
																					    + sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
																					    - sum_decayed_W11_W11_term*pow(W20_dot_W20,2) )
															+ (lower_boundary_n20-interim_lower_boundary_n20)*( -pow(sum_decayed_W11_W20_term,2)*sum_decayed_W20_W11_term
																					    + sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term
																					    - sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
																					    + sum_decayed_W20_W11_term*pow(W11_dot_W20,2)
																					    - sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W20
																					    + sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W20 )
														      ) / det;
	  viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( sum_decayed_W11_W20_term*pow(sum_decayed_W20_W11_term,2)
																 - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																 - sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
																 - sum_decayed_W11_W20_term*pow(W20_dot_W11,2)
																 + sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
																 + sum_decayed_W11_W11_term*W20_dot_W11*W20_dot_W20 )
										+ (upper_boundary_n20-interim_upper_boundary_n20)*( -sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term
																    + pow(sum_decayed_W11_W11_term,2)*sum_decayed_W20_W20_term
																    - sum_decayed_W20_W20_term*pow(W11_dot_W11,2)
																    + sum_decayed_W20_W11_term*W11_dot_W11*W11_dot_W20
																    - sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W11
																    + sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W11 )
										+ (lower_boundary_n11-interim_lower_boundary_n11)*( sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11
																    - pow(sum_decayed_W20_W11_term,2)*W11_dot_W20
																    + sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W20_dot_W11
																    + W11_dot_W20*pow(W20_dot_W11,2)
																    - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W20_dot_W20
																    - W11_dot_W11*W20_dot_W11*W20_dot_W20 )
										+ (lower_boundary_n20-interim_lower_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11
																    - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20
																    + sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11
																    - W11_dot_W11*W11_dot_W20*W20_dot_W11
																    - pow(sum_decayed_W11_W11_term,2)*W20_dot_W20
																    + pow(W11_dot_W11,2)*W20_dot_W20 )
									     ) / det;
	  viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11
																					   + pow(sum_decayed_W20_W11_term,2)*W11_dot_W20
																					   - sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W20_dot_W11
																					   - W11_dot_W20*pow(W20_dot_W11,2)
																					   + sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W20_dot_W20
																					   + W11_dot_W11*W20_dot_W11*W20_dot_W20 )
															+ (upper_boundary_n20-interim_upper_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11
																					    - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20
																					    + sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11
																					    - W11_dot_W11*W11_dot_W20*W20_dot_W11
																					    - pow(sum_decayed_W11_W11_term,2)*W20_dot_W20
																					    + pow(W11_dot_W11,2)*W20_dot_W20 )
															+ (lower_boundary_n11-interim_lower_boundary_n11)*( -sum_decayed_W11_W20_term*pow(sum_decayed_W20_W11_term,2)
																					    + sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																					    + sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
																					    + sum_decayed_W11_W20_term*pow(W20_dot_W11,2)
																					    - sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
																					    - sum_decayed_W11_W11_term*W20_dot_W11*W20_dot_W20 )
															+ (lower_boundary_n20-interim_lower_boundary_n20)*( -sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term
																					    + pow(sum_decayed_W11_W11_term,2)*sum_decayed_W20_W20_term
																					    - sum_decayed_W20_W20_term*pow(W11_dot_W11,2)
																					    + sum_decayed_W20_W11_term*W11_dot_W11*W11_dot_W20
																					    - sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W11
																					    + sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W11 )
														      ) / det;
	#endif
      }
  }

  #ifdef BC_HEATFLUX
    y_broadcast(heatflux_transients_factors, (mesh->xend-mesh->xstart+1)*(mesh->LocalNz)*2, 0);
  #endif
  #ifdef BC_VISCOSITY
    y_broadcast(viscosity_transients_factors, (mesh->xend-mesh->xstart+1)*(mesh->LocalNz)*2, 0);
  #endif
  
  for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
    for (int jz=0; jz<mesh->LocalNz; jz++) {
      #if defined(BC_HEATFLUX) && !defined(BC_VISCOSITY)
	for (int i=0; i<number_of_negative_eigenvalues; i++) {
	  #ifdef CALCULATE_HEATFLUX
	    heatflux_lower_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    heatflux_upper_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_VISCOSITY
	    viscosity_lower_boundary_transients[i][jx][jz] = W20_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    viscosity_upper_boundary_transients[i][jx][jz] = -W20_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_FRICTION
	    friction_lower_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    friction_upper_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	}
      #elif defined(BC_VISCOSITY) && !defined(BC_HEATFLUX)
	for (int i=0; i<number_of_negative_eigenvalues; i++) {
	  #ifdef CALCULATE_HEATFLUX
	    heatflux_lower_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    heatflux_upper_boundary_transients[i][jx][jz] = -W11_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_VISCOSITY
	    viscosity_lower_boundary_transients[i][jx][jz] = W20_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    viscosity_upper_boundary_transients[i][jx][jz] = W20_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_FRICTION
	    friction_lower_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    friction_upper_boundary_transients[i][jx][jz] = -C10_1k_dot_W1k_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	#endif
	}
      #elif defined(BC_HEATFLUX) && defined(BC_VISCOSITY)
	for (int i=0; i<number_of_negative_eigenvalues; i++) {
	  #ifdef CALCULATE_HEATFLUX
	    heatflux_lower_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							    + W11_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    heatflux_upper_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							    - W11_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_VISCOSITY
	    viscosity_lower_boundary_transients[i][jx][jz] = W20_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							      + W20_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    viscosity_upper_boundary_transients[i][jx][jz] = -W20_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							      + W20_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_FRICTION
	    friction_lower_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							    + C10_1k_dot_W1k_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    friction_upper_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							    - C10_1k_dot_W1k_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	}
      #else
	for (int i=0; i<number_of_negative_eigenvalues; i++) {
	  #ifdef CALCULATE_HEATFLUX
	    heatflux_lower_boundary_transients[i][jx][jz] = 0.;
	    heatflux_upper_boundary_transients[i][jx][jz] = 0.;
	  #endif
	  #ifdef CALCULATE_VISCOSITY
	    viscosity_lower_boundary_transients[i][jx][jz] = 0.;
	    viscosity_upper_boundary_transients[i][jx][jz] = 0.;
	  #endif
	  #ifdef CALCULATE_FRICTION
	    friction_lower_boundary_transients[i][jx][jz] = 0.;
	    friction_upper_boundary_transients[i][jx][jz] = 0.;
	  #endif
	}
      #endif
    }

  #ifndef NOEDGETERMS
  start_index(position);
  do {
    position->jy=mesh->ystart;
    calc_index(position);
    do {
      for (int i=0; i<number_of_negative_eigenvalues; i++) {
	BoutReal exp_increasing = exp(increasing_dimensionless_length[*position]/eigenvalues[i]);
	BoutReal exp_decreasing = exp(decreasing_dimensionless_length[*position]/eigenvalues[i]);
	#ifdef CALCULATE_HEATFLUX
	  electron_heat_flux[*position] += heatflux_lower_boundary_transients[i][position->jx][position->jz] * exp_increasing
						  + heatflux_upper_boundary_transients[i][position->jx][position->jz] * exp_decreasing;
	#endif
	#ifdef CALCULATE_VISCOSITY
	  electron_viscosity[*position] += viscosity_lower_boundary_transients[i][position->jx][position->jz] * exp_increasing
						  + viscosity_upper_boundary_transients[i][position->jx][position->jz] * exp_decreasing;
	#endif
	#ifdef CALCULATE_FRICTION
	  electron_friction[*position] += friction_lower_boundary_transients[i][position->jx][position->jz] * exp_increasing
						  + friction_upper_boundary_transients[i][position->jx][position->jz] * exp_decreasing;
	#endif
      }
      position->jy++;
      calc_index(position);
    } while (position->jy<mesh->yend+1);
  } while (next_indexperp(position));
  #endif
  
  #ifdef CALCULATE_HEATFLUX
    electron_heat_flux *= -5./4.*sqrt(2./electron_mass)*(T_electron^1.5); //now we have q=-5/4*v_Telectron*T_electron*n^(1,1)
    mesh->communicate(electron_heat_flux);
  #endif
  #ifdef CALCULATE_VISCOSITY
    electron_viscosity *= T_electron;
    mesh->communicate(electron_viscosity);
  #endif
  #ifdef CALCULATE_FRICTION
    electron_friction *= 2.*T_electron/3.*lambdaC_inverse;
    electron_friction += -2.*sqrt(2.*electron_mass*T_electron)*lambdaC_inverse*(-jpar/electron_charge); // Need to include also the friction due directly to the Maxwellian part of the distribution function
    mesh->communicate(electron_friction);
  #endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NonLocalParallel::calculate_nonlocal_closures_cell_ylow(const Field3D &n_electron, const Field3D &T_electron
								#ifdef DRIVE_GRADV
								  , const Field3D &V_electron
								#endif
								#ifdef DRIVE_VEMINUSVI
								  , const Field3D &jpar
								#endif
								#ifdef BC_HEATFLUX
								  , const Field3D &heat_flux_boundary_condition
								#endif
								#ifdef BC_VISCOSITY
								  , const Field3D &viscosity_boundary_condition
								#endif
							      ) {

  Coordinates *coord = mesh->coordinates();
  
  lambdaC_inverse = n_electron * pow(electron_charge,4) * logLambda / 12 / pow(PI,1.5) / pow(epsilon_0,2) / (T_electron^2);
  
  #ifdef DRIVE_GRADT
    gradT_driveterm = 5./4. * n_electron / T_electron / lambdaC_inverse; //g^(1,1)
    gradT_electron = Grad_par(T_electron,CELL_YLOW);
  #endif
  #ifdef DRIVE_GRADV
    gradV_driveterm = -0.5 * n_electron / sqrt(2.*T_electron/electron_mass) * Grad_par(V_electron,CELL_CENTRE) / lambdaC_inverse;
//     gradV_electron = Grad_par(V_electron,CELL_YLOW); // would be more consistent to put this on CELL_CENTRE, but that would make imposing a boundary condition a pain.
  #endif
  #ifdef DRIVE_VEMINUSVI
    VeminusVi_driveterm = -2./sqrt(PI) * (-1./electron_charge) / sqrt(2.*T_electron/electron_mass);
    // jpar is already CELL_YLOW.
  #endif
  // Calculate target boundary guard cell derivitives (at YLOW) for gradT_electron with 4th order forward/backward differences from T_electron (at CENTRE)
  // Also check for unphysical lambdaC_inverse
  for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
    for (int jz=0; jz<mesh->LocalNz; jz++) {
      for (int jy=mesh->ystart-1; jy>=0; jy--) {
	#ifdef DRIVE_GRADT
	gradT_electron(rlow.ind,jy,jz)=(-93.*T_electron(rlow.ind,jy,jz) + 229.*T_electron(rlow.ind,jy+1,jz) - 225.*T_electron(rlow.ind,jy+2,jz) + 111.*T_electron(rlow.ind,jy+3,jz) - 22.*T_electron(rlow.ind,jy+4,jz))/48./coord->dy(rlow.ind,jy)/sqrt((coord->g_22(rlow.ind,jy) + coord->g_22(rlow.ind,jy+1) + coord->g_22(rlow.ind,jy+2) + coord->g_22(rlow.ind,jy+3) + coord->g_22(rlow.ind,jy+4))/5.);
	if (abs(gradT_driveterm(rlow.ind,jy,jz))>1.e37 || gradT_driveterm(rlow.ind,jy,jz)!=gradT_driveterm(rlow.ind,jy,jz)) gradT_driveterm(rlow.ind,jy,jz) = 1.e37;
	#endif
	#ifdef DRIVE_GRADV
	  // Nothing to be done here: gradV_driveterm is CELL_CENTRE and the guard cell values are not used
	#endif
	#ifdef DRIVE_VEMINUSVI
	  if (abs(VeminusVi_driveterm[rlow.ind][jy][jz])>1.e37 || VeminusVi_driveterm[rlow.ind][jy][jz]!=VeminusVi_driveterm[rlow.ind][jy][jz]) VeminusVi_driveterm[rlow.ind][jy][jz] = 1.e37;
	#endif
	if (abs(lambdaC_inverse[rlow.ind][jy][jz])>1.e37 || lambdaC_inverse[rlow.ind][jy][jz]!=lambdaC_inverse[rlow.ind][jy][jz]) lambdaC_inverse[rlow.ind][jy][jz] = 1.e37;
      }
    }

  

  for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
    for (int jz=0; jz<mesh->LocalNz; jz++) {
      #ifdef DRIVE_GRADT
	gradT_electron[rup.ind][mesh->yend][jz] = (T_electron[rup.ind][mesh->yend-2][jz]-27.*T_electron[rup.ind][mesh->yend-1][jz]+27.*T_electron[rup.ind][mesh->yend][jz]-T_electron[rup.ind][mesh->yend+1][jz])/24./coord->dy[rup.ind][mesh->yend+1]/sqrt((coord->g_22[rup.ind][mesh->yend-1] + coord->g_22[rup.ind][mesh->yend] + coord->g_22[rup.ind][mesh->yend+1] + coord->g_22[rup.ind][mesh->yend+2])/4.);
      #endif
      for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++) {
	#ifdef DRIVE_GRADT
	gradT_electron[rup.ind][jy][jz]=(93.*T_electron(rup.ind,jy-1,jz) - 229.*T_electron(rup.ind,jy-2,jz) + 225.*T_electron(rup.ind,jy-3,jz) - 111.*T_electron(rup.ind,jy-4,jz) + 22.*T_electron(rup.ind,jy-5,jz))/48./coord->dy(rup.ind,jy-1)/sqrt((coord->g_22[rup.ind][jy-1] + coord->g_22[rup.ind][jy-2] + coord->g_22[rup.ind][jy-3] + coord->g_22[rup.ind][jy-4] + coord->g_22[rup.ind][jy-5])/5.);
	  if (abs(gradT_driveterm[rup.ind][jy][jz])>1.e37 || gradT_driveterm[rup.ind][jy][jz]!=gradT_driveterm[rup.ind][jy][jz]) gradT_driveterm[rup.ind][jy][jz] = 1.e37;
	#endif
	#ifdef DRIVE_GRADV
	  // Nothing to be done here: gradV_driveterm is CELL_CENTRE and the guard cell values are not used
	#endif
	#ifdef DRIVE_VEMINUSVI
	  if (abs(VeminusVi_driveterm[rup.ind][jy][jz])>1.e37 || VeminusVi_driveterm[rup.ind][jy][jz]!=VeminusVi_driveterm[rup.ind][jy][jz]) VeminusVi_driveterm[rup.ind][jy][jz] = 1.e37;
	#endif
	if (abs(lambdaC_inverse[rup.ind][jy][jz])>1.e37 || lambdaC_inverse[rup.ind][jy][jz]!=lambdaC_inverse[rup.ind][jy][jz]) lambdaC_inverse[rup.ind][jy][jz] = 1.e37;
      }
    }
        
  #ifdef DRIVE_GRADT
    mesh->communicate(gradT_electron);
  #endif
  #ifdef DRIVE_GRADV
    mesh->communicate(gradV_driveterm);
  #endif
  // No gradient term to communicate for VeminusVi
  
  // Now calculate z and deltaz everywhere
  cubic_spline_inverse_lambdaC.calculate(lambdaC_inverse);
  
  start_index(position);
  
  if (mesh->UpXSplitIndex()!=0 || mesh->DownXSplitIndex()!=0) throw BoutException("This code cannot currently handle x-splitting of processors.");
    if (!is_lower_boundary) {
      FieldPerp pass_dimensionless_length;
      pass_dimensionless_length.allocate();
      {
	mesh->wait(mesh->irecvYInOutdest(*pass_dimensionless_length.getData(),mesh->LocalNx*(mesh->LocalNz),
					NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz));
      }
      pass_dimensionless_length.setIndex(mesh->ystart);
      increasing_dimensionless_length = pass_dimensionless_length;
    }
    
  do {
    position->jy = mesh->ystart-1;
    calc_index(position);
    interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
    // dimensionless_length_deltas_above[jy] and dimensionless_length_deltas_below[jy] are the deltaz's for the half-step above and below, respectively, the CELL_CENTRE at jy
    // deltaz between position[jy](CELL_CENTRE), where t=0, and position[jyp](CELL_YLOW), where t=0.5
    
    Coordinates *coord = mesh->coordinates();
    
    dimensionless_length_deltas_above[position->jx][position->jy][position->jz] = coord->dy(position->jx,position->jy)*sqrt(coord->g_22(position->jx,position->jy))
										    *(interp_coefficients[0]/2. + interp_coefficients[1]/2./4. + interp_coefficients[2]/3./8. + interp_coefficients[3]/4./16.);
    // deltaz between position[jyp](CELL_YLOW), where t=0.5, and position[jyp](CELL_CENTRE), where t=1
    dimensionless_length_deltas_below[position->jx][position->jyp][position->jz] = coord->dy(position->jx,position->jy)*sqrt(coord->g_22(position->jx,position->jyp))
										    *(interp_coefficients[0]/2. + interp_coefficients[1]/2.*3./4. + interp_coefficients[2]/3.*7./8. + interp_coefficients[3]/4.*15./16.);
    next_index_y(position);
    
    do{
      interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
      // deltaz between position[jy](CELL_CENTRE), where t=0, and position[jyp](CELL_YLOW), where t=0.5
      
      dimensionless_length_deltas_above[position->jx][position->jy][position->jz] = coord->dy(position->jx,position->jy)*sqrt(coord->g_22(position->jx,position->jy))
										      *(interp_coefficients[0]/2. + interp_coefficients[1]/2./4. + interp_coefficients[2]/3./8. + interp_coefficients[3]/4./16.);
      // deltaz between position[jyp](CELL_YLOW), where t=0.5, and position[jyp](CELL_CENTRE), where t=1
      dimensionless_length_deltas_below[position->jx][position->jyp][position->jz] = coord->dy(position->jx,position->jyp)*sqrt(coord->g_22(position->jx,position->jyp))
										      *(interp_coefficients[0]/2. + interp_coefficients[1]/2.*3./4. + interp_coefficients[2]/3.*7./8. + interp_coefficients[3]/4.*15./16.);
      increasing_dimensionless_length[position->jx][position->jyp][position->jz] = increasing_dimensionless_length[*position] + dimensionless_length_deltas_below[*position] + dimensionless_length_deltas_above[*position];
      } while (next_index_y(position));
    
  } while (next_indexperp(position));
  
  {
    Timer timer("comms");
    mesh->sendYOutOutdest(*increasing_dimensionless_length.slice(mesh->yend+1).getData(),mesh->LocalNx*(mesh->LocalNz),
			    NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz);
  }
  
  // Send the total dimensionless_length at the upper boundary back to the other processors.
  if (is_upper_boundary) {
    total_dimensionless_length = increasing_dimensionless_length.slice(mesh->yend);
  }

  y_broadcast(*total_dimensionless_length.getData(), mesh->LocalNx*mesh->LocalNz, mesh->getNYPE()-1);
  
  decreasing_dimensionless_length = -increasing_dimensionless_length;
  for (int jy=mesh->ystart; jy<=mesh->yend; jy++) {
    total_dimensionless_length.setIndex(jy);
    decreasing_dimensionless_length += total_dimensionless_length;
  }
  
  #ifdef CALCULATE_HEATFLUX
    electron_heat_flux = 0.;
  #endif
  #ifdef CALCULATE_VISCOSITY
    electron_viscosity = 0.;
  #endif
  #ifdef CALCULATE_FRICTION
    electron_friction = 0.;
  #endif
  #ifdef DRIVE_GRADT
    #ifdef CALCULATE_HEATFLUX
      electron_heat_flux += -heatflux_gradT_zerocoeff * interp_to(gradT_driveterm,CELL_YLOW) * gradT_electron; //zero eigenvalue contribution to n^(1,1)/T^1.5
    #endif
    // viscosity gets no zero eigenvalue contribution
    #ifdef CALCULATE_FRICTION
      electron_friction += -friction_gradT_zerocoeff * interp_to(gradT_driveterm,CELL_YLOW) * gradT_electron; //zero eigenvalue contribution
    #endif
    cubic_spline_gradT_driveterm.calculate(gradT_driveterm);
    cubic_spline_gradT.calculate(gradT_electron);
  #endif
  #ifdef DRIVE_GRADV
    // gradV drive gives no zero eigenvalue contribution
    cubic_spline_gradV_driveterm.calculate(gradV_driveterm);
  #endif
  #ifdef DRIVE_VEMINUSVI
    #ifdef CALCULATE_HEATFLUX
      electron_heat_flux += -heatflux_VeminusVi_zerocoeff * interp_to(VeminusVi_driveterm,CELL_YLOW)*jpar;
    #endif
    // viscosity gets no zero eigenvalue contribution
    #ifdef CALCULATE_FRICTION
      electron_friction += -friction_VeminusVi_zerocoeff * interp_to(VeminusVi_driveterm,CELL_YLOW)*jpar;
    #endif
    cubic_spline_VeminusVi_driveterm.calculate(VeminusVi_driveterm);
    cubic_spline_jpar.calculate(jpar);
  #endif
  
  for (int i=0; i<number_of_negative_eigenvalues; i++) {
    #ifdef DRIVE_GRADT
      integration.calculateIntegralBelow_cell_ylow(eigenvalues[i], dimensionless_length_deltas_below, dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_gradT_driveterm, cubic_spline_gradT, i+0);
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += heatflux_gradT_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += viscosity_gradT_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += friction_gradT_coefficients[i]*integration.integral_below;
      #endif
    #endif
    #ifdef DRIVE_GRADV
      integration.calculateIntegralBelow_cell_ylow(eigenvalues[i], dimensionless_length_deltas_below, dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_gradV_driveterm, cubic_spline_one, i+number_of_negative_eigenvalues);
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += heatflux_gradV_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += viscosity_gradV_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += friction_gradV_coefficients[i]*integration.integral_below;
      #endif
    #endif
    #ifdef DRIVE_VEMINUSVI
      integration.calculateIntegralBelow_cell_ylow(eigenvalues[i], dimensionless_length_deltas_below,dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_VeminusVi_driveterm, cubic_spline_jpar, i+2*number_of_negative_eigenvalues); // The driveterm here is CELL_YLOW while there is no CELL_CENTRE quantity mulitplying it, so use a cubic spline of a constant field with value 1. everywhere instead (would be slightly more efficient to write a different calculateIntegralBelow_cell_ylow, but tedious).
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += heatflux_VeminusVi_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += viscosity_VeminusVi_coefficients[i]*integration.integral_below;
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += friction_VeminusVi_coefficients[i]*integration.integral_below;
      #endif
    #endif
  }
  for (int i=0; i<number_of_negative_eigenvalues; i++) {
    #ifdef DRIVE_GRADT
      integration.calculateIntegralAbove_cell_ylow(eigenvalues[i], dimensionless_length_deltas_below, dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_gradT_driveterm, cubic_spline_gradT, i+0);
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += -heatflux_gradT_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += -viscosity_gradT_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += -friction_gradT_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
    #endif
    #ifdef DRIVE_GRADV
      integration.calculateIntegralAbove_cell_ylow(eigenvalues[i], dimensionless_length_deltas_below, dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_gradV_driveterm, cubic_spline_one, i+number_of_negative_eigenvalues);
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += -heatflux_gradV_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += -viscosity_gradV_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += -friction_gradV_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
    #endif
    #ifdef DRIVE_VEMINUSVI
      integration.calculateIntegralAbove_cell_ylow(eigenvalues[i], dimensionless_length_deltas_below,dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, cubic_spline_VeminusVi_driveterm, cubic_spline_jpar, i+2*number_of_negative_eigenvalues); // The driveterm here is CELL_YLOW while there is no CELL_CENTRE quantity mulitplying it, so use a cubic spline of a constant field with value 1. everywhere instead (would be slightly more efficient to write a different calculateIntegralAbove_cell_ylow, but tedious).
      #ifdef CALCULATE_HEATFLUX
	electron_heat_flux += -heatflux_VeminusVi_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
      #ifdef CALCULATE_VISCOSITY
	electron_viscosity += -viscosity_VeminusVi_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
      #ifdef CALCULATE_FRICTION
	electron_friction += -friction_VeminusVi_coefficients[i]*integration.integral_above; //gives \hat(n)^A //gives \hat(n)^A/(T_electron^1.5)
      #endif
    #endif
  }
  
  for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
    for (int jz=0; jz<mesh->LocalNz; jz++) {
      position->jx=rup.ind;
      position->jy=mesh->yend;
      position->jz=jz;
      calc_index(position);
      BoutReal Te_here = interp_to_point_YLOW(T_electron,*position);
      #ifdef BC_HEATFLUX
	pass_interim_upper_boundary_n11[rup.ind][jz] = electron_heat_flux[rup.ind][mesh->yend][jz];
	upper_boundary_condition_n11[rup.ind][jz] = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition[rup.ind][mesh->yend][jz];
      #endif
      #ifdef BC_VISCOSITY
	pass_interim_upper_boundary_n20[rup.ind][jz] = electron_viscosity[rup.ind][mesh->yend][jz];
	upper_boundary_condition_n20[rup.ind][jz] = viscosity_boundary_condition[rup.ind][mesh->yend][jz]/Te_here;
      #endif
    }
  if (is_upper_boundary && mesh->getNYPE()>1) {
    #ifdef BC_HEATFLUX
      MPI_Request request1 = mesh->sendToProc(mesh->getXProcIndex(),0,
					      *pass_interim_upper_boundary_n11.getData(),
					      mesh->LocalNx*mesh->LocalNz,
					      NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex());
      MPI_Request request2 = mesh->sendToProc(mesh->getXProcIndex(),0,
					      *upper_boundary_condition_n11.getData(),
					      mesh->LocalNx*mesh->LocalNz,
					      NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 1);
      MPI_Waitall(1,&request1,MPI_STATUSES_IGNORE);
      MPI_Waitall(1,&request2,MPI_STATUSES_IGNORE);
    #endif
    #ifdef BC_VISCOSITY
      MPI_Request request3 = mesh->sendToProc(mesh->getXProcIndex(),0,
					      *pass_interim_upper_boundary_n20.getData(),
					      mesh->LocalNx*mesh->LocalNz,
					      NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 3);
      MPI_Request request4 = mesh->sendToProc(mesh->getXProcIndex(),0,
					      *upper_boundary_condition_n20.getData(),
					      mesh->LocalNx*mesh->LocalNz,
					      NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 4);
      MPI_Waitall(1,&request3,MPI_STATUSES_IGNORE);
      MPI_Waitall(1,&request4,MPI_STATUSES_IGNORE);
    #endif
  }
  if (is_lower_boundary && mesh->getNYPE()>1) {
    #ifdef BC_HEATFLUX
      mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(), mesh->getNYPE()-1,
					    *pass_interim_upper_boundary_n11.getData(),
					    mesh->LocalNx*mesh->LocalNz,
					    NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex()) );
      mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(), mesh->getNYPE()-1,
					    *upper_boundary_condition_n11.getData(),
					    mesh->LocalNx*mesh->LocalNz,
					    NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 1) );
    #endif
    #ifdef BC_VISCOSITY
      mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(), mesh->getNYPE()-1,
					    *pass_interim_upper_boundary_n20.getData(),
					    mesh->LocalNx*mesh->LocalNz,
					    NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 3) );
      mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(), mesh->getNYPE()-1,
					    *upper_boundary_condition_n20.getData(),
					    mesh->LocalNx*mesh->LocalNz,
					    NONLOCAL_PARALLEL_TAGBASE + mesh->getXProcIndex() + 4) );
    #endif
  }
  
  if (is_lower_boundary) {
    for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	position->jx=jx;
	position->jy=mesh->ystart;
	position->jz=jz;
	calc_index(position);
	BoutReal Te_here = interp_to_point_YLOW(T_electron,*position);
	#if defined(BC_HEATFLUX) && !defined(BC_VISCOSITY)
	  BoutReal lower_boundary_n11 = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition[jx][mesh->ystart][jz];
	  BoutReal upper_boundary_n11 = upper_boundary_condition_n11[jx][jz];
	  BoutReal interim_lower_boundary_n11 = electron_heat_flux[jx][mesh->ystart][jz];
	  BoutReal interim_upper_boundary_n11 = pass_interim_upper_boundary_n11[jx][jz];
	  /*
	  electron_heat_flux is, at this point, the contribution to n11 from nhat_plus.
	  We want the actual heat flux at mesh->ystart to be boundary_heat_flux.
	  Thus the remainder must come from nhat_minus, which we will construct here just to give the right 1,1 component (could set number_of_negative_eigenvalues-1 more components if desired)
	  However, the transients from the other boundary do not necessarily decay to vanishing by the time they get to this boundary, so we must solve for both.
	  Fortunately this reduces to a single algebraic equation which makes it surprisingly easy to do (in the case of just a single condition being imposed at least).
	  */
	  
	  BoutReal sum_decayed_W11_W11_term = 0.;
	  for (int i=0; i<number_of_negative_eigenvalues; i++) {
	    exp_total_dimensionless_length_over_eigenvalue[i] = exp(total_dimensionless_length[jx][jz]/eigenvalues[i]);
	    sum_decayed_W11_W11_term += W11_B_times_WinverseB_11[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	  }
	  heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (lower_boundary_n11 - interim_lower_boundary_n11)*W11_dot_W11
								    - sum_decayed_W11_W11_term*(upper_boundary_n11 - interim_upper_boundary_n11) )
								    / ( pow(W11_dot_W11,2) - pow(sum_decayed_W11_W11_term,2) );
	  heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11 - interim_upper_boundary_n11)*W11_dot_W11
													    - sum_decayed_W11_W11_term*(lower_boundary_n11 - interim_lower_boundary_n11) )
													    / ( pow(W11_dot_W11,2) - pow(sum_decayed_W11_W11_term,2) );
	#elif defined(BC_VISCOSITY) && !defined(BC_HEATFLUX)
	  BoutReal lower_boundary_n20 = viscosity_boundary_condition[jx][mesh->ystart][jz]/Te_here;
	  BoutReal upper_boundary_n20 = upper_boundary_condition_n20[jx][jz];
	  BoutReal interim_lower_boundary_n20 = electron_viscosity[jx][mesh->ystart][jz];
	  BoutReal interim_upper_boundary_n20 = pass_interim_upper_boundary_n20[jx][jz];
	  /*
	  electron_viscosity is, at this point, the contribution to n20 from nhat_plus.
	  We want the actual viscosity at mesh->ystart to be boundary_viscosity.
	  Thus the remainder must come from nhat_minus, which we will construct here just to give the right 2,0 component (could set number_of_negative_eigenvalues-1 more components if desired)
	  However, the transients from the other boundary do not necessarily decay to vanishing by the time they get to this boundary, so we must solve for both.
	  Fortunately this reduces to a single algebraic equation which makes it surprisingly easy to do (in the case of just a single condition being imposed at least).
	  */
	  
	  BoutReal sum_decayed_W20_W20_term = 0.;
	  for (int i=0; i<number_of_negative_eigenvalues; i++) {
	    exp_total_dimensionless_length_over_eigenvalue[i] = exp(total_dimensionless_length[jx][jz]/eigenvalues[i]);
	    sum_decayed_W20_W20_term += W20_B_times_WinverseB_20[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	  }
	  viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (lower_boundary_n20 - interim_lower_boundary_n20)*W20_dot_W20
								    - sum_decayed_W20_W20_term*(upper_boundary_n20 - interim_upper_boundary_n20) )
								    / ( pow(W20_dot_W20,2) - pow(sum_decayed_W20_W20_term,2) );
	  viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n20 - interim_upper_boundary_n20)*W20_dot_W20
													    - sum_decayed_W20_W20_term*(lower_boundary_n20 - interim_lower_boundary_n20) )
													    / ( pow(W20_dot_W20,2) - pow(sum_decayed_W20_W20_term,2));
	#elif defined(BC_HEATFLUX) && defined(BC_VISCOSITY)
	  BoutReal lower_boundary_n11 = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition[jx][mesh->ystart][jz];
	  BoutReal upper_boundary_n11 = upper_boundary_condition_n11[jx][jz];
	  BoutReal interim_lower_boundary_n11 = electron_heat_flux[jx][mesh->ystart][jz];
	  BoutReal interim_upper_boundary_n11 = pass_interim_upper_boundary_n11[jx][jz];
	  BoutReal lower_boundary_n20 = viscosity_boundary_condition[jx][mesh->ystart][jz]/Te_here;
	  BoutReal upper_boundary_n20 = upper_boundary_condition_n20[jx][jz];
	  BoutReal interim_lower_boundary_n20 = electron_viscosity[jx][mesh->ystart][jz];
	  BoutReal interim_upper_boundary_n20 = pass_interim_upper_boundary_n20[jx][jz];
	  BoutReal sum_decayed_W11_W11_term = 0.;
	  BoutReal sum_decayed_W20_W20_term = 0.;
	  BoutReal sum_decayed_W11_W20_term = 0.;
	  BoutReal sum_decayed_W20_W11_term = 0.;
	  for (int i=0; i<number_of_negative_eigenvalues; i++) {
	    exp_total_dimensionless_length_over_eigenvalue[i] = exp(total_dimensionless_length[jx][jz]/eigenvalues[i]);
	    sum_decayed_W11_W11_term += W11_B_times_WinverseB_11[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	    sum_decayed_W20_W20_term += W20_B_times_WinverseB_20[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	    sum_decayed_W20_W11_term += W20_B_times_WinverseB_11[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	    sum_decayed_W11_W20_term += W11_B_times_WinverseB_20[i]*exp_total_dimensionless_length_over_eigenvalue[i];
	  }
	  BoutReal det = pow(sum_decayed_W11_W20_term,2)*pow(sum_decayed_W20_W11_term,2)
			 - 2*sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
			 + pow(sum_decayed_W11_W11_term,2)*pow(sum_decayed_W20_W20_term,2)
			 - pow(sum_decayed_W20_W20_term,2)*pow(W11_dot_W11,2)
			 + 2*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
			 - pow(sum_decayed_W20_W11_term,2)*pow(W11_dot_W20,2)
			 - 2*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
			 + 2*sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
			 - pow(sum_decayed_W11_W20_term,2)*pow(W20_dot_W11,2)
			 + pow(W11_dot_W20,2)*pow(W20_dot_W11,2)
			 + 2*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
			 - 2*sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
			 + 2*sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
			 - 2*W11_dot_W11*W11_dot_W20*W20_dot_W11*W20_dot_W20
			 - pow(sum_decayed_W11_W11_term,2)*pow(W20_dot_W20,2)
			 + pow(W11_dot_W11,2)*pow(W20_dot_W20,2);
	  heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																+ sum_decayed_W11_W11_term*pow(sum_decayed_W20_W20_term,2)
																+ sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
																- sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
																+ sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
																- sum_decayed_W11_W11_term*pow(W20_dot_W20,2) )
										+ (upper_boundary_n20-interim_upper_boundary_n20)*( pow(sum_decayed_W11_W20_term,2)*sum_decayed_W20_W11_term
																    - sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term
																    + sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
																    - sum_decayed_W20_W11_term*pow(W11_dot_W20,2)
																    + sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W20
																    - sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W20 )
										+ (lower_boundary_n11-interim_lower_boundary_n11)*( -pow(sum_decayed_W20_W20_term,2)*W11_dot_W11
																    + sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																    - sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W20_dot_W11
																    + sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W20_dot_W20
																    - W11_dot_W20*W20_dot_W11*W20_dot_W20 + W11_dot_W11*pow(W20_dot_W20,2) )
										+ (lower_boundary_n20-interim_lower_boundary_n20)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11
																    + sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																    - pow(sum_decayed_W11_W20_term,2)*W20_dot_W11
																    + pow(W11_dot_W20,2)*W20_dot_W11
																    + sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W20
																    - W11_dot_W11*W11_dot_W20*W20_dot_W20 )
									    ) / det;
	  heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -pow(sum_decayed_W20_W20_term,2)*W11_dot_W11
																					  + sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																					  - sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W20_dot_W11
																					  + sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W20_dot_W20
																					  - W11_dot_W20*W20_dot_W11*W20_dot_W20
																					  + W11_dot_W11*pow(W20_dot_W20,2) )
															+ (upper_boundary_n20-interim_upper_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11
																					    - sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																					    + pow(sum_decayed_W11_W20_term,2)*W20_dot_W11
																					    - pow(W11_dot_W20,2)*W20_dot_W11
																					    - sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W20
																					    + W11_dot_W11*W11_dot_W20*W20_dot_W20 )
															+ (lower_boundary_n11-interim_lower_boundary_n11)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																					    + sum_decayed_W11_W11_term*pow(sum_decayed_W20_W20_term,2)
																					    + sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
																					    - sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
																					    + sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
																					    - sum_decayed_W11_W11_term*pow(W20_dot_W20,2) )
															+ (lower_boundary_n20-interim_lower_boundary_n20)*( -pow(sum_decayed_W11_W20_term,2)*sum_decayed_W20_W11_term
																					    + sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term
																					    - sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
																					    + sum_decayed_W20_W11_term*pow(W11_dot_W20,2)
																					    - sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W20
																					    + sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W20 )
														      ) / det;
	  viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( sum_decayed_W11_W20_term*pow(sum_decayed_W20_W11_term,2)
																 - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																 - sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
																 - sum_decayed_W11_W20_term*pow(W20_dot_W11,2)
																 + sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
																 + sum_decayed_W11_W11_term*W20_dot_W11*W20_dot_W20 )
										+ (upper_boundary_n20-interim_upper_boundary_n20)*( -sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term
																    + pow(sum_decayed_W11_W11_term,2)*sum_decayed_W20_W20_term
																    - sum_decayed_W20_W20_term*pow(W11_dot_W11,2)
																    + sum_decayed_W20_W11_term*W11_dot_W11*W11_dot_W20
																    - sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W11
																    + sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W11 )
										+ (lower_boundary_n11-interim_lower_boundary_n11)*( sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11
																    - pow(sum_decayed_W20_W11_term,2)*W11_dot_W20
																    + sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W20_dot_W11
																    + W11_dot_W20*pow(W20_dot_W11,2)
																    - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W20_dot_W20
																    - W11_dot_W11*W20_dot_W11*W20_dot_W20 )
										+ (lower_boundary_n20-interim_lower_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11
																    - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20
																    + sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11
																    - W11_dot_W11*W11_dot_W20*W20_dot_W11
																    - pow(sum_decayed_W11_W11_term,2)*W20_dot_W20
																    + pow(W11_dot_W11,2)*W20_dot_W20 )
									     ) / det;
	  viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11
																					   + pow(sum_decayed_W20_W11_term,2)*W11_dot_W20
																					   - sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W20_dot_W11
																					   - W11_dot_W20*pow(W20_dot_W11,2)
																					   + sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W20_dot_W20
																					   + W11_dot_W11*W20_dot_W11*W20_dot_W20 )
															+ (upper_boundary_n20-interim_upper_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11
																					    - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20
																					    + sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11
																					    - W11_dot_W11*W11_dot_W20*W20_dot_W11
																					    - pow(sum_decayed_W11_W11_term,2)*W20_dot_W20
																					    + pow(W11_dot_W11,2)*W20_dot_W20 )
															+ (lower_boundary_n11-interim_lower_boundary_n11)*( -sum_decayed_W11_W20_term*pow(sum_decayed_W20_W11_term,2)
																					    + sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																					    + sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
																					    + sum_decayed_W11_W20_term*pow(W20_dot_W11,2)
																					    - sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
																					    - sum_decayed_W11_W11_term*W20_dot_W11*W20_dot_W20 )
															+ (lower_boundary_n20-interim_lower_boundary_n20)*( -sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term
																					    + pow(sum_decayed_W11_W11_term,2)*sum_decayed_W20_W20_term
																					    - sum_decayed_W20_W20_term*pow(W11_dot_W11,2)
																					    + sum_decayed_W20_W11_term*W11_dot_W11*W11_dot_W20
																					    - sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W11
																					    + sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W11 )
														      ) / det;
	#endif
      }
  }

  #ifdef BC_HEATFLUX
    y_broadcast(heatflux_transients_factors, (mesh->xend-mesh->xstart+1)*(mesh->LocalNz)*2, 0);
  #endif
  #ifdef BC_VISCOSITY
    y_broadcast(viscosity_transients_factors, (mesh->xend-mesh->xstart+1)*(mesh->LocalNz)*2, 0);
  #endif
  
  for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
    for (int jz=0; jz<mesh->LocalNz; jz++) {
      #if defined(BC_HEATFLUX) && !defined(BC_VISCOSITY)
	for (int i=0; i<number_of_negative_eigenvalues; i++) {
	  #ifdef CALCULATE_HEATFLUX
	    heatflux_lower_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    heatflux_upper_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_VISCOSITY
	    viscosity_lower_boundary_transients[i][jx][jz] = W20_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    viscosity_upper_boundary_transients[i][jx][jz] = -W20_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_FRICTION
	    friction_lower_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    friction_upper_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	}
      #elif defined(BC_VISCOSITY) && !defined(BC_HEATFLUX)
	for (int i=0; i<number_of_negative_eigenvalues; i++) {
	  #ifdef CALCULATE_HEATFLUX
	    heatflux_lower_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    heatflux_upper_boundary_transients[i][jx][jz] = -W11_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_VISCOSITY
	    viscosity_lower_boundary_transients[i][jx][jz] = W20_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    viscosity_upper_boundary_transients[i][jx][jz] = W20_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_FRICTION
	    friction_lower_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    friction_upper_boundary_transients[i][jx][jz] = -C10_1k_dot_W1k_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	#endif
	}
      #elif defined(BC_HEATFLUX) && defined(BC_VISCOSITY)
	for (int i=0; i<number_of_negative_eigenvalues; i++) {
	  #ifdef CALCULATE_HEATFLUX
	    heatflux_lower_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							    + W11_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    heatflux_upper_boundary_transients[i][jx][jz] = W11_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							    - W11_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_VISCOSITY
	    viscosity_lower_boundary_transients[i][jx][jz] = W20_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							      + W20_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    viscosity_upper_boundary_transients[i][jx][jz] = -W20_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							      + W20_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	  #ifdef CALCULATE_FRICTION
	    friction_lower_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_11[i]*heatflux_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							    + C10_1k_dot_W1k_B_times_WinverseB_20[i]*viscosity_transients_factors[(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	    friction_upper_boundary_transients[i][jx][jz] = C10_1k_dot_W1k_B_times_WinverseB_11[i]*heatflux_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz]
							    - C10_1k_dot_W1k_B_times_WinverseB_20[i]*viscosity_transients_factors[(mesh->xend-mesh->xstart+1)*(mesh->LocalNz)+(jx-mesh->xstart)*(mesh->LocalNz)+jz];
	  #endif
	}
      #else
	for (int i=0; i<number_of_negative_eigenvalues; i++) {
	  #ifdef CALCULATE_HEATFLUX
	    heatflux_lower_boundary_transients[i][jx][jz] = 0.;
	    heatflux_upper_boundary_transients[i][jx][jz] = 0.;
	  #endif
	  #ifdef CALCULATE_VISCOSITY
	    viscosity_lower_boundary_transients[i][jx][jz] = 0.;
	    viscosity_upper_boundary_transients[i][jx][jz] = 0.;
	  #endif
	  #ifdef CALCULATE_FRICTION
	    friction_lower_boundary_transients[i][jx][jz] = 0.;
	    friction_upper_boundary_transients[i][jx][jz] = 0.;
	  #endif
	}
      #endif
    }

  #ifndef NOEDGETERMS
  start_index(position);
  do {
    position->jy=mesh->ystart;
    calc_index(position);
    do {
      for (int i=0; i<number_of_negative_eigenvalues; i++) {
	BoutReal exp_increasing = exp(increasing_dimensionless_length[*position]/eigenvalues[i]);
	BoutReal exp_decreasing = exp(decreasing_dimensionless_length[*position]/eigenvalues[i]);
	#ifdef CALCULATE_HEATFLUX
	  electron_heat_flux[*position] += heatflux_lower_boundary_transients[i][position->jx][position->jz] * exp_increasing
						  + heatflux_upper_boundary_transients[i][position->jx][position->jz] * exp_decreasing;
	#endif
	#ifdef CALCULATE_VISCOSITY
	  electron_viscosity[*position] += viscosity_lower_boundary_transients[i][position->jx][position->jz] * exp_increasing
						  + viscosity_upper_boundary_transients[i][position->jx][position->jz] * exp_decreasing;
	#endif
	#ifdef CALCULATE_FRICTION
	  electron_friction[*position] += friction_lower_boundary_transients[i][position->jx][position->jz] * exp_increasing
						  + friction_upper_boundary_transients[i][position->jx][position->jz] * exp_decreasing;
	#endif
      }
      position->jy++;
      calc_index(position);
    } while (position->jy<mesh->yend+1);
  } while (next_indexperp(position));
  #endif
  
  #ifdef CALCULATE_HEATFLUX
    electron_heat_flux *= -5./4.*sqrt(2./electron_mass)*interp_to(T_electron^1.5,CELL_YLOW); //now we have q=-5/4*v_Telectron*T_electron*n^(1,1)
    mesh->communicate(electron_heat_flux);
  #endif
  #ifdef CALCULATE_VISCOSITY
    electron_viscosity *= T_electron;
    mesh->communicate(electron_viscosity);
  #endif
  #ifdef CALCULATE_FRICTION
    electron_friction *= 2.*T_electron/3.*lambdaC_inverse;
    electron_friction += -2.*sqrt(2.*electron_mass*T_electron)*lambdaC_inverse*(-jpar/electron_charge);
  #endif
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NonLocalParallel::set_boundary_gradients() {
  for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
    for (int jz=0; jz<mesh->LocalNz; jz++) {
      #ifdef CALCULATE_HEATFLUX
//       BoutReal heat_flux_boundarygradient = (electron_heat_flux[rlow.ind][mesh->ystart][jz]-27.*electron_heat_flux[rlow.ind][mesh->ystart+1][jz]+27.*electron_heat_flux[rlow.ind][mesh->ystart+2][jz]-electron_heat_flux[rlow.ind][mesh->ystart+3][jz])/24.; // NB gradient in index space
//       BoutReal heat_flux_boundarygradient = (-11.*electron_heat_flux[rlow.ind][mesh->ystart][jz] + 18.*electron_heat_flux[rlow.ind][mesh->ystart+1][jz] - 9.*electron_heat_flux[rlow.ind][mesh->ystart+2][jz] + 2.*electron_heat_flux[rlow.ind][mesh->ystart+3][jz]) / 6. / coord->dy[rlow.ind][mesh->ystart] / sqrt((coord->g_22[rlow.ind][mesh->ystart]+coord->g_22[rlow.ind][mesh->ystart+1]+coord->g_22[rlow.ind][mesh->ystart+2]+coord->g_22[rlow.ind][mesh->ystart+3])/4.);
	BoutReal heat_flux_boundarygradient = (-electron_heat_flux[rlow.ind][mesh->ystart][jz] + electron_heat_flux[rlow.ind][mesh->ystart+boundary_gradient_smoothing_length][jz])/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      #endif
      #ifdef CALCULATE_VISCOSITY
//       BoutReal viscosity_boundarygradient = (electron_viscosity[rlow.ind][mesh->ystart][jz]-27.*electron_viscosity[rlow.ind][mesh->ystart+1][jz]+27.*electron_viscosity[rlow.ind][mesh->ystart+2][jz]-electron_viscosity[rlow.ind][mesh->ystart+3][jz])/24.; // NB gradient in index space
//       BoutReal viscosity_boundarygradient = (-11.*electron_viscosity[rlow.ind][mesh->ystart][jz] + 18.*electron_viscosity[rlow.ind][mesh->ystart+1][jz] - 9.*electron_viscosity[rlow.ind][mesh->ystart+2][jz] + 2.*electron_viscosity[rlow.ind][mesh->ystart+3][jz]) / 6. / coord->dy[rlow.ind][mesh->ystart] / sqrt((coord->g_22[rlow.ind][mesh->ystart]+coord->g_22[rlow.ind][mesh->ystart+1]+coord->g_22[rlow.ind][mesh->ystart+2]+coord->g_22[rlow.ind][mesh->ystart+3])/4.);
	BoutReal viscosity_boundarygradient = (-electron_viscosity[rlow.ind][mesh->ystart][jz] + electron_viscosity[rlow.ind][mesh->ystart+boundary_gradient_smoothing_length][jz])/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      #endif
      #ifdef CALCULATE_FRICTION
//       BoutReal friction_boundarygradient = (electron_friction[rlow.ind][mesh->ystart][jz]-27.*electron_friction[rlow.ind][mesh->ystart+1][jz]+27.*electron_friction[rlow.ind][mesh->ystart+2][jz]-electron_friction[rlow.ind][mesh->ystart+3][jz])/24.; // NB gradient in index space
//       BoutReal friction_boundarygradient = (-11.*electron_friction[rlow.ind][mesh->ystart][jz] + 18.*electron_friction[rlow.ind][mesh->ystart+1][jz] - 9.*electron_friction[rlow.ind][mesh->ystart+2][jz] + 2.*electron_friction[rlow.ind][mesh->ystart+3][jz]) / 6. / coord->dy[rlow.ind][mesh->ystart] / sqrt((coord->g_22[rlow.ind][mesh->ystart]+coord->g_22[rlow.ind][mesh->ystart+1]+coord->g_22[rlow.ind][mesh->ystart+2]+coord->g_22[rlow.ind][mesh->ystart+3])/4.);
	BoutReal friction_boundarygradient = (-electron_friction[rlow.ind][mesh->ystart][jz] + electron_friction[rlow.ind][mesh->ystart+boundary_gradient_smoothing_length][jz])/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      #endif
      for (int jy=mesh->ystart-1; jy>=0; jy--) {
	#ifdef CALCULATE_HEATFLUX
	  electron_heat_flux[rlow.ind][jy][jz] = electron_heat_flux[rlow.ind][jy+1][jz] - heat_flux_boundarygradient;
	#endif
	#ifdef CALCULATE_VISCOSITY
	  electron_viscosity[rlow.ind][jy][jz] = electron_viscosity[rlow.ind][jy+1][jz] - viscosity_boundarygradient;
	#endif
	#ifdef CALCULATE_FRICTION
	  electron_friction[rlow.ind][jy][jz] = electron_friction[rlow.ind][jy+1][jz] - friction_boundarygradient;
	#endif
      }
    }
  for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
    for (int jz=0; jz<mesh->LocalNz; jz++) {
      #ifdef CALCULATE_HEATFLUX
// 	BoutReal heat_flux_boundarygradient = (electron_heat_flux[rup.ind][mesh->yend-3][jz]-27.*electron_heat_flux[rup.ind][mesh->yend-2][jz]+27.*electron_heat_flux[rup.ind][mesh->yend-1][jz]-electron_heat_flux[rup.ind][mesh->yend][jz])/24.; // NB gradient in index space
//       BoutReal heat_flux_boundarygradient = (11.*electron_heat_flux[rup.ind][mesh->yend][jz] - 18.*electron_heat_flux[rup.ind][mesh->yend-1][jz] + 9.*electron_heat_flux[rup.ind][mesh->yend-2][jz] - 2.*electron_heat_flux[rup.ind][mesh->yend-3][jz]) / 6. / coord->dy[rup.ind][mesh->yend] / sqrt((coord->g_22[rup.ind][mesh->yend]+coord->g_22[rup.ind][mesh->yend-1]+coord->g_22[rup.ind][mesh->yend-2]+coord->g_22[rup.ind][mesh->yend-3])/4.);
	BoutReal heat_flux_boundarygradient = (-electron_heat_flux[rup.ind][mesh->yend-boundary_gradient_smoothing_length][jz] + electron_heat_flux[rup.ind][mesh->yend][jz])/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      #endif
      #ifdef CALCULATE_VISCOSITY
// 	BoutReal viscosity_boundarygradient = (electron_viscosity[rup.ind][mesh->yend-3][jz]-27.*electron_viscosity[rup.ind][mesh->yend-2][jz]+27.*electron_viscosity[rup.ind][mesh->yend-1][jz]-electron_viscosity[rup.ind][mesh->yend][jz])/24.; // NB gradient in index space
//       BoutReal viscosity_boundarygradient = (11.*electron_viscosity[rup.ind][mesh->yend][jz] - 18.*electron_viscosity[rup.ind][mesh->yend-1][jz] + 9.*electron_viscosity[rup.ind][mesh->yend-2][jz] - 2.*electron_viscosity[rup.ind][mesh->yend-3][jz]) / 6. / coord->dy[rup.ind][mesh->yend] / sqrt((coord->g_22[rup.ind][mesh->yend]+coord->g_22[rup.ind][mesh->yend-1]+coord->g_22[rup.ind][mesh->yend-2]+coord->g_22[rup.ind][mesh->yend-3])/4.);
	BoutReal viscosity_boundarygradient = (-electron_viscosity[rup.ind][mesh->yend-boundary_gradient_smoothing_length][jz] + electron_viscosity[rup.ind][mesh->yend][jz])/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      #endif
      #ifdef CALCULATE_FRICTION
// 	BoutReal friction_boundarygradient = (electron_friction[rup.ind][mesh->yend-3][jz]-27.*electron_friction[rup.ind][mesh->yend-2][jz]+27.*electron_friction[rup.ind][mesh->yend-1][jz]-electron_friction[rup.ind][mesh->yend][jz])/24.; // NB gradient in index space
//       BoutReal friction_boundarygradient = (11.*electron_friction[rup.ind][mesh->yend][jz] - 18.*electron_friction[rup.ind][mesh->yend-1][jz] + 9.*electron_friction[rup.ind][mesh->yend-2][jz] - 2.*electron_friction[rup.ind][mesh->yend-3][jz]) / 6. / coord->dy[rup.ind][mesh->yend] / sqrt((coord->g_22[rup.ind][mesh->yend]+coord->g_22[rup.ind][mesh->yend-1]+coord->g_22[rup.ind][mesh->yend-2]+coord->g_22[rup.ind][mesh->yend-3])/4.);
	BoutReal friction_boundarygradient = (-electron_friction[rup.ind][mesh->yend-boundary_gradient_smoothing_length][jz] + electron_friction[rup.ind][mesh->yend][jz])/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      #endif
      for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++) {
	#ifdef CALCULATE_HEATFLUX
	  electron_heat_flux[rup.ind][jy][jz] = electron_heat_flux[rup.ind][jy-1][jz] + heat_flux_boundarygradient;
	#endif
	#ifdef CALCULATE_VISCOSITY
	  electron_viscosity[rup.ind][jy][jz] = electron_viscosity[rup.ind][jy-1][jz] + viscosity_boundarygradient;
	#endif
	#ifdef CALCULATE_FRICTION
	  electron_friction[rup.ind][jy][jz] = electron_friction[rup.ind][jy-1][jz] + friction_boundarygradient;
	#endif
      }
    }
}

void NonLocalParallel::set_neumann_boundary_conditions() {
  for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
    for (int jy=0; jy<mesh->ystart; jy++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	#ifdef CALCULATE_HEATFLUX
	  electron_heat_flux[rlow.ind][jy][jz]= electron_heat_flux[rlow.ind][mesh->ystart][jz];
	#endif
	#ifdef CALCULATE_VISCOSITY
	  electron_viscosity[rlow.ind][jy][jz]= electron_viscosity[rlow.ind][mesh->ystart][jz];
	#endif
	#ifdef CALCULATE_FRICTION
	  electron_friction[rlow.ind][jy][jz]= electron_friction[rlow.ind][mesh->ystart][jz];
	#endif
      }
  for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
    for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++)
      for (int jz=0; jz<mesh->LocalNz; jz++) {
	#ifdef CALCULATE_HEATFLUX
	  electron_heat_flux[rup.ind][jy][jz]= electron_heat_flux[rup.ind][mesh->yend][jz];
	#endif
	#ifdef CALCULATE_VISCOSITY
	  electron_viscosity[rup.ind][jy][jz]= electron_viscosity[rup.ind][mesh->yend][jz];
	#endif
	#ifdef CALCULATE_FRICTION
	  electron_friction[rup.ind][jy][jz]= electron_friction[rup.ind][mesh->yend][jz];
	#endif
      }
}

void NonLocalParallel::y_broadcast(void* input_buffer, const int &size, const int &root_processor) {
  // NB Assumes that the mesh is BoutMesh
  Timer timer("comms");
  
  /// NOTE: This only works if there are no branch-cuts
  MPI_Comm comm_inner = mesh->getYcomm(0);
  
//  MPI_Bcast(input_buffer, size, PVEC_REAL_MPI_TYPE, root_processor, comm_yprocs);
  MPI_Bcast(input_buffer, size, MPI_DOUBLE, root_processor, comm_inner);
// Should use commented out version if method is transferred to boutmesh.cxx
}

void NonLocalParallel::rms_over_y(const Field3D &input_field, FieldPerp &output_field) {
  FieldPerp tempsum;
  tempsum = 0.;
  int ye = mesh->yend;
  if (mesh->StaggerGrids && input_field.getLocation()==CELL_CENTRE) ye--;
  for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
    for (int jz=0; jz<mesh->LocalNz; jz++)
      for (int jy=mesh->ystart; jy<=ye; jy++) {
	tempsum(jx,jz) += SQ(input_field(jx,jy,jz));
      }
  
  /// NOTE: This only works if there are no branch-cuts
  MPI_Comm comm_inner = mesh->getYcomm(0);
  
  MPI_Reduce(*tempsum.getData(),
	     *output_field.getData(),
	     mesh->LocalNx*mesh->LocalNz,
	     MPI_DOUBLE,
	     MPI_SUM,
	     mesh->getXProcIndex(), // Why?
	     comm_inner);

  // Don't really understand what this bit is supposed to do.
  if (mesh->getYProcIndex()==0) {
    int ny = mesh->GlobalNy;
    if (mesh->StaggerGrids && input_field.getLocation()==CELL_CENTRE) ny--;
    for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
      for (int jz=0; jz<mesh->LocalNz;jz++)
	output_field[jx][jz] = sqrt(output_field[jx][jz]/ny);
    mesh->sendToProc(mesh->getXProcIndex(),mesh->getNYPE()-1,*output_field.getData(),mesh->LocalNx*mesh->LocalNz,NONLOCAL_PARALLEL_TAGBASE);
  }
  else if (mesh->getYProcIndex()==mesh->getNYPE()-1) {
    mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(),0,*output_field.getData(),mesh->LocalNx*mesh->LocalNz,NONLOCAL_PARALLEL_TAGBASE));
  }
}

/*
  Calculates a mean over the Y (parallel) direction. 
  
  Should probably be replaced by a call to averageY (src/physics/smoothing.cxx)
 */
void NonLocalParallel::mean_over_y(const Field3D &input_field, FieldPerp &output_field, int exclude_edgecells) {
  FieldPerp tempsum;
  tempsum = 0.;
  int ys = mesh->ystart+exclude_edgecells;
  int ye = mesh->yend-exclude_edgecells;
  
  if (mesh->StaggerGrids && input_field.getLocation()==CELL_CENTRE && mesh->lastY()) ye--;
  
  for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
    for (int jz=0; jz<mesh->LocalNz; jz++)
      for (int jy=ys; jy<=ye; jy++) {
	tempsum(jx,jz)+=input_field(jx,jy,jz);
      }
  
  /// NOTE: This only works if there are no branch-cuts
  MPI_Comm comm_inner = mesh->getYcomm(0);

  MPI_Reduce(*tempsum.getData(),
	     *output_field.getData(),
	     mesh->LocalNx*mesh->LocalNz,
	     MPI_DOUBLE,
	     MPI_SUM,
	     mesh->getXProcIndex(),  // Why? 
	     comm_inner);
  if (mesh->getYProcIndex()==0) {
    int ny = mesh->GlobalNy;
    if (mesh->StaggerGrids && input_field.getLocation()==CELL_CENTRE) ny--;
    ny-=2*exclude_edgecells;
    for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
      for (int jz=0; jz<mesh->LocalNz;jz++)
	output_field[jx][jz] = output_field[jx][jz]/ny;
    mesh->sendToProc(mesh->getXProcIndex(),mesh->getNYPE()-1,*output_field.getData(),mesh->LocalNx*mesh->LocalNz,NONLOCAL_PARALLEL_TAGBASE);
  }
  else if (mesh->getYProcIndex()==mesh->getNYPE()-1) {
    mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(),0,*output_field.getData(),mesh->LocalNx*mesh->LocalNz,NONLOCAL_PARALLEL_TAGBASE));
  }
}

BoutReal NonLocalParallel::interp_to_point_YLOW(const Field3D &input, bindex &position) {
   if(mesh->StaggerGrids)
     return (9.*(input(position.jx,position.jym,position.jz)+input(position.jx,position.jy,position.jz))-(input(position.jx,position.jy2m,position.jz)+input(position.jx,position.jyp,position.jz)))/16.;
   else
     return input[position];
}
