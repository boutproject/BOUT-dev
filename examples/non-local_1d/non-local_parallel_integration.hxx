/**************************************************************************
 * Perform integral needed to calculate non-local electron closures
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

#ifndef __NONLOCALPARALLELINTEGRATION_H__
#define __NONLOCALPARALLELINTEGRATION_H__

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

// #include "cubic_spline.hxx"
#include "cubic_spline_local.hxx"

/// Class for non-local heat flux integration
class NonLocalParallelIntegration {
public:
  NonLocalParallelIntegration();
  ~NonLocalParallelIntegration();
  void initialise(const bool pass_electron_heat_flux_location_is_ylow);
  
  Field3D integral_below;
  Field3D integral_above;
  void	calculateIntegralBelow_cell_centre(const BoutReal &eigenvalue, const Field3D &dimensionless_length_deltas_above, CubicSpline &cubic_spline_inverse_lambdaC, CubicSpline &cubic_spline_drive_term, const int &counter);
  void	calculateIntegralAbove_cell_centre(const BoutReal &eigenvalue, const Field3D &dimensionless_length_deltas_above, CubicSpline &cubic_spline_inverse_lambdaC, CubicSpline &cubic_spline_drive_term, const int &counter);
  void	calculateIntegralBelow_cell_ylow(const BoutReal &eigenvalue, const Field3D &dimensionless_length_deltas_below, const Field3D &dimensionless_length_deltas_above, CubicSpline &cubic_spline_inverse_lambdaC, CubicSpline &cubic_spline_drive_term, CubicSpline &cubic_spline_gradT, const int &counter);
  void	calculateIntegralAbove_cell_ylow(const BoutReal &eigenvalue, const Field3D &dimensionless_length_deltas_below, const Field3D &dimensionless_length_deltas_above, CubicSpline &cubic_spline_inverse_lambdaC, CubicSpline &cubic_spline_drive_term, CubicSpline &cubic_spline_gradT, const int &counter);
 
private:
  bindex * position;
  BoutReal * deltal;
  BoutReal * interp_coeffs_lambdaC_inverse;
  BoutReal * interp_coeffs_drive_term;
  BoutReal * interp_coeffs_gradT;
  BoutReal * integral_coeffs;
  BoutReal * integral_parts;
  int NONLOCAL_PARALLEL_INTEGRATION_TAGBASE;
  bool electron_heat_flux_location_is_ylow;
};

#endif //__NONLOCALPARALLELINTEGRATION_H__
