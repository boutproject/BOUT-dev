/**************************************************************************
 * Calculate non-local electron closures
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
 **************************************************************************/

#ifndef __NONLOCALPARALLEL_H__
#define __NONLOCALPARALLEL_H__

#define CALCULATE_HEATFLUX
#define CALCULATE_VISCOSITY
#define CALCULATE_FRICTION

#define BC_HEATFLUX
#define BC_VISCOSITY

#define DRIVE_GRADT
#define DRIVE_GRADV
#define DRIVE_VEMINUSVI

#include <globals.hxx>

#include <sstream>
#include <cmath>

#include <bout_types.hxx>
#include <options.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <bout/sys/timer.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>
#include <fieldperp.hxx>
#include <field3d.hxx>
#include <field2d.hxx>
#include <stencils.hxx>
#include <interpolation.hxx>

#include "cubic_spline_local.hxx"
#include "non-local_parallel_integration.hxx"

class NonLocalParallel {
public:
//   NonLocalParallel();
  ~NonLocalParallel();
  int moments_number;
  void initialise(BoutReal pass_electron_charge, BoutReal pass_electron_mass, BoutReal pass_ion_mass, BoutReal pass_epsilon_0, BoutReal pass_logLambda, const bool pass_fluxes_location_is_ylow=false, BoutReal pass_gamma_factor=5./3.);
  #ifdef CALCULATE_HEATFLUX
    Field3D electron_heat_flux;
  #endif
  #ifdef CALCULATE_VISCOSITY
    Field3D electron_viscosity;
  #endif
  #ifdef CALCULATE_FRICTION
    Field3D electron_friction;
  #endif
  void calculate_nonlocal_closures(const Field3D &n_electron, const Field3D &T_electron
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
				   );
  void set_boundary_gradients();
  void set_neumann_boundary_conditions();
  
  // Bit of a hack: should put the method somewhere more sensible (like in boutmesh.cxx)
  void y_broadcast(void* input_buffer, const int &size, const int &root_processor); // NB Assumes that the mesh is BoutMesh
  
  void rms_over_y(const Field3D &input_field, FieldPerp &output_field);
  void mean_over_y(const Field3D &input_field, FieldPerp &output_field, int exclude_edgecells=0);
  BoutReal interp_to_point_YLOW(const Field3D &input, bindex &position);
  
private:
  Field3D lambdaC_inverse;	// inverse collision length, i.e. 1/lambdaC
  Field3D increasing_dimensionless_length;	// d/dl(increasing_dimensionless_length)=1/lambdaC
  Field3D decreasing_dimensionless_length;
  Field3D dimensionless_length_deltas_above; // change in dimensionless_length between jy and jy+1, used in calculating integrals. (More accurate than first adding up and then taking differences)
  Field3D dimensionless_length_deltas_below;
  FieldPerp total_dimensionless_length;
  NonLocalParallelIntegration integration;
  int number_of_negative_eigenvalues;
  BoutReal * eigenvalues;
  BoutReal * exp_total_dimensionless_length_over_eigenvalue;
  BoutReal electron_charge;
  BoutReal electron_mass;
  BoutReal ion_mass;
  BoutReal epsilon_0;
  BoutReal logLambda;
  int boundary_gradient_smoothing_length;
  BoutReal boundary_condition_smoothing_range;
  CubicSpline cubic_spline_inverse_lambdaC;
  BoutReal * interp_coefficients;
  #ifdef DRIVE_GRADT
    Field3D gradT_driveterm;	// drive from Maxwellian moments for heat flux calculation. g^(1,1)/T^1.5 as a function of l.
					// NB slightly bad notation: drive_term will INCLUDE grad(T) for cell-centred version but EXCLUDE grad(T) (i.e. be only the cell-centred prefactor) for the ylow staggered version, so that grad(T) can be interpolated separately
    Field3D gradT_electron;	// only used for staggered-grids case
    CubicSpline cubic_spline_gradT_driveterm;
    CubicSpline cubic_spline_gradT;
  #endif
  #ifdef DRIVE_GRADV
    Field3D gradV_driveterm;
    Field3D gradV_electron;
    CubicSpline cubic_spline_gradV_driveterm;
    CubicSpline cubic_spline_one;
  #endif
  #ifdef DRIVE_VEMINUSVI
    Field3D VeminusVi_driveterm;
    CubicSpline cubic_spline_VeminusVi_driveterm;
    CubicSpline cubic_spline_jpar;
  #endif
  #ifdef CALCULATE_HEATFLUX
    #ifdef DRIVE_GRADT
      BoutReal heatflux_gradT_zerocoeff;
      BoutReal * heatflux_gradT_coefficients;
    #endif
    #ifdef DRIVE_GRADV
      BoutReal * heatflux_gradV_coefficients;
    #endif
    #ifdef DRIVE_VEMINUSVI
      BoutReal heatflux_VeminusVi_zerocoeff;
      BoutReal * heatflux_VeminusVi_coefficients;
    #endif
    FieldPerp * heatflux_lower_boundary_transients;
    FieldPerp * heatflux_upper_boundary_transients;
    #ifdef BC_VISCOSITY
      BoutReal * W11_B_times_WinverseB_20;
      BoutReal W11_dot_W20;
    #endif
  #endif
  #ifdef BC_HEATFLUX
    BoutReal * W11_B_times_WinverseB_11;
    BoutReal W11_dot_W11;
    BoutReal * heatflux_transients_factors;
    FieldPerp pass_interim_upper_boundary_n11;
    FieldPerp upper_boundary_condition_n11;
  #endif
  #ifdef CALCULATE_VISCOSITY
    #ifdef DRIVE_GRADT
      BoutReal * viscosity_gradT_coefficients;
    #endif
    #ifdef DRIVE_GRADV
      BoutReal * viscosity_gradV_coefficients;
    #endif
    #ifdef DRIVE_VEMINUSVI
      BoutReal * viscosity_VeminusVi_coefficients;
    #endif
    FieldPerp * viscosity_lower_boundary_transients;
    FieldPerp * viscosity_upper_boundary_transients;
    #ifdef BC_HEATFLUX
      BoutReal * W20_B_times_WinverseB_11;
      BoutReal W20_dot_W11;
    #endif
  #endif
  #ifdef BC_VISCOSITY
    BoutReal * W20_B_times_WinverseB_20;
    BoutReal W20_dot_W20;
    BoutReal * viscosity_transients_factors;
    FieldPerp pass_interim_upper_boundary_n20;
    FieldPerp upper_boundary_condition_n20;
  #endif
  #ifdef CALCULATE_FRICTION
    #ifdef DRIVE_GRADT
      BoutReal friction_gradT_zerocoeff;
      BoutReal * friction_gradT_coefficients;
    #endif
    #ifdef DRIVE_GRADV
      BoutReal * friction_gradV_coefficients;
    #endif
    #ifdef DRIVE_VEMINUSVI
      BoutReal friction_VeminusVi_zerocoeff;
      BoutReal * friction_VeminusVi_coefficients;
    #endif
    #ifdef BC_HEATFLUX
      BoutReal * C10_1k_dot_W1k_B_times_WinverseB_11;
    #endif
    #ifdef BC_VISCOSITY
      BoutReal * C10_1k_dot_W1k_B_times_WinverseB_20;
    #endif
    FieldPerp * friction_lower_boundary_transients;
    FieldPerp * friction_upper_boundary_transients;
  #endif
  int NONLOCAL_PARALLEL_TAGBASE;
  MPI_Request broadcast_request;
  MPI_Comm comm_yprocs_minusone;
  bindex* position;
  bool fluxes_location_is_ylow;
  bool is_lower_boundary;
  bool is_upper_boundary;
  void calculate_nonlocal_closures_cell_centre(const Field3D &n_electron, const Field3D &T_electron
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
  );
  void calculate_nonlocal_closures_cell_ylow(const Field3D &n_electron, const Field3D &T_electron
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
  );
  #if CHECK > 0
    bool calculated_before_setting_bcs;
  #endif
};

#endif // __NONLOCALPARALLEL_H__
