/*!
 * \file heat_flux_integration.cxx
 *
 * \brief Perform the integral needed to calculate the non-local heat flux
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

#include "heat_flux_integration.hxx"

/**********************************************************************************
 *                    INTEGRATION INITIALISATION AND CREATION
 **********************************************************************************/

HeatFluxIntegration::HeatFluxIntegration() {
  position = new bindex;
  deltal = new BoutReal;
  integral_coeffs = new BoutReal[4];
  integral_parts = new BoutReal[4];
  }

HeatFluxIntegration::~HeatFluxIntegration() {
  delete position;
  delete deltal;
  delete [] integral_coeffs;
  delete [] integral_parts;
}

void HeatFluxIntegration::initialise(const bool pass_electron_heat_flux_location_is_ylow=false) {
  electron_heat_flux_location_is_ylow = pass_electron_heat_flux_location_is_ylow;
  
  if (electron_heat_flux_location_is_ylow) {
    integral_below.setLocation(CELL_YLOW);
    integral_above.setLocation(CELL_YLOW);
  }
  
  *deltal = 0.;
  for (int i=0; i<4; i++) {
    integral_coeffs[i]=0.;
    integral_parts[i]=0.;
  }
  integral_below = 0.0;
  integral_above = 0.0;
  
  // Get the options for the model
  Options *options = Options::getRoot()->getSection("electron_heat_flux");
  OPTION(options, HEATFLUX_INTEGRATION_TAGBASE, 16381);
}


/**********************************************************************************
 *                            INTEGRATION ROUTINES
 **********************************************************************************/

void HeatFluxIntegration::calculateIntegralBelow_cell_centre(BoutReal eigenvalue, const Field3D &dimensionless_length_deltas_above, CubicSpline &cubic_spline_inverse_lambdaC, CubicSpline &cubic_spline_drive_term, const int &counter) {
  TRACE("HeatFluxIntegration::calculateIntegralBelow()");

  start_index(position);
  do {
    position->jy=mesh->ystart-1;
    calc_index(position);
    if (mesh->firstY())
      next_index_y(position);
    else {
      // Set the value at ystart-1 equal to the value at yend on the previous processor.
      if (position->jx<mesh->DownXSplitIndex()) {
	mesh->wait(mesh->irecvYInIndest(&integral_below[*position],1,HEATFLUX_INTEGRATION_TAGBASE + mesh->DownXSplitIndex()*mesh->LocalNz*counter + mesh->LocalNz*position->jx + position->jz));
      }
      else {
	mesh->wait(mesh->irecvYInOutdest(&integral_below[*position],1,HEATFLUX_INTEGRATION_TAGBASE + (mesh->LocalNx-mesh->DownXSplitIndex())*mesh->LocalNz*counter + mesh->LocalNz*(position->jx-mesh->DownXSplitIndex()) + position->jz));
      }
    }
    do {
      *deltal = mesh->dy[position->jx][position->jyp]*sqrt((mesh->g_22[position->jx][position->jy]+mesh->g_22[position->jx][position->jyp])/2.);
      
      interp_coeffs_drive_term = cubic_spline_drive_term.coefficients(position);
      interp_coeffs_lambdaC_inverse = cubic_spline_inverse_lambdaC.coefficients(position);
      
      //calculate the coefficients of drive_term expanded in z rather than l (from the expansion of drive_term in l and l in z [inverted from z in l to third order])
      integral_coeffs[0] = interp_coeffs_drive_term[0];
      integral_coeffs[1] = interp_coeffs_drive_term[1] / *deltal/interp_coeffs_lambdaC_inverse[0];
      integral_coeffs[2] = interp_coeffs_drive_term[2] / pow(*deltal*interp_coeffs_lambdaC_inverse[0],2)
			    - interp_coeffs_drive_term[1] * interp_coeffs_lambdaC_inverse[1]/2./pow(*deltal,2)/pow(interp_coeffs_lambdaC_inverse[0],3);
      integral_coeffs[3] = interp_coeffs_drive_term[3] / pow(*deltal*interp_coeffs_lambdaC_inverse[0],3)
			    - interp_coeffs_drive_term[2] * interp_coeffs_lambdaC_inverse[1]/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4)
			    + interp_coeffs_drive_term[1] * (0.5*pow(interp_coeffs_lambdaC_inverse[1],2)/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],5) - interp_coeffs_lambdaC_inverse[2]/3./pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4));

      // Calculate analytically the integral from jy to jyp of z^n*exp( (zupper-z)/zeta ), firstly allowing for zeta to be large and secondly allowing for it to be small.
      BoutReal exp_delta_over_eigenvalue = exp(dimensionless_length_deltas_above[*position]/eigenvalue);
      if (dimensionless_length_deltas_above[*position]>abs(eigenvalue)) {
	integral_parts[0] = eigenvalue * (exp_delta_over_eigenvalue-1);
	integral_parts[1] = pow(eigenvalue,2) * (exp_delta_over_eigenvalue-1)
			      - eigenvalue * dimensionless_length_deltas_above[*position];
	integral_parts[2] = 2*pow(eigenvalue,3) * (exp_delta_over_eigenvalue-1)
			      - 2*pow(eigenvalue,2) * dimensionless_length_deltas_above[*position] - eigenvalue * pow(dimensionless_length_deltas_above[*position],2);
	integral_parts[3] = 6*pow(eigenvalue,4) * (exp_delta_over_eigenvalue-1)
			      - 6*pow(eigenvalue,3) * dimensionless_length_deltas_above[*position] - 3*pow(eigenvalue,2) * pow(dimensionless_length_deltas_above[*position],2) - eigenvalue * pow(dimensionless_length_deltas_above[*position],3);
      }
      else {
	BoutReal exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue = exp(dimensionless_length_deltas_above[*position]/eigenvalue/2.) * 2.*sinh(dimensionless_length_deltas_above[*position]/eigenvalue/2.);
	integral_parts[0] = eigenvalue * exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue;
	integral_parts[1] = pow(eigenvalue,2) * exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      - eigenvalue* dimensionless_length_deltas_above[*position];
	integral_parts[2] = 2*pow(eigenvalue,3) * exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      - 2*pow(eigenvalue,2) * dimensionless_length_deltas_above[*position] - eigenvalue * pow(dimensionless_length_deltas_above[*position],2);
	integral_parts[3] = 6*pow(eigenvalue,4) * exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      - 6*pow(eigenvalue,3) * dimensionless_length_deltas_above[*position] - 3*pow(eigenvalue,2) * pow(dimensionless_length_deltas_above[*position],2) - eigenvalue * pow(dimensionless_length_deltas_above[*position],3);
      }

      //add up the contributions to the integral at jy first from the integral at jy-1 and then from the expansion of the integral between jy-1 and jy
      integral_below[position->jx][position->jyp][position->jz] = integral_below[*position] * exp_delta_over_eigenvalue;
      for (int i=0; i<4; i++) integral_below[position->jx][position->jyp][position->jz] += integral_coeffs[i]*integral_parts[i];
      
      next_index_y(position);
    } while (position->jy < mesh->yend);
    
    // Send the value at yend to the next processor.
    if (position->jx < mesh->UpXSplitIndex()) {
      Timer timer("comms");
      mesh->sendYOutIndest(&integral_below[*position],1,HEATFLUX_INTEGRATION_TAGBASE + mesh->UpXSplitIndex()*mesh->LocalNz*counter + mesh->LocalNz*position->jx + position->jz);
    }
    else {
      Timer timer("comms");
      mesh->sendYOutOutdest(&integral_below[*position],1,HEATFLUX_INTEGRATION_TAGBASE + (mesh->LocalNx-mesh->UpXSplitIndex())*mesh->LocalNz*counter + mesh->LocalNz*(position->jx-mesh->UpXSplitIndex()) + position->jz);
    }
    
  } while (next_indexperp(position));
}

void HeatFluxIntegration::calculateIntegralAbove_cell_centre(BoutReal eigenvalue, const Field3D &dimensionless_length_deltas_above, CubicSpline &cubic_spline_inverse_lambdaC, CubicSpline &cubic_spline_drive_term, const int &counter) {
  TRACE("HeatFluxIntegration::calculateIntegralAbove()");

  start_index_lasty(position);
  do {
    position->jy=mesh->yend;
    calc_index(position);
    if (mesh->lastY())
      previous_index_y(position);
    else {
      // Set the value at yend+1 to the value at ystart on the previous processor.
      if (position->jx<mesh->UpXSplitIndex()) {
	mesh->wait(mesh->irecvYOutIndest(&integral_above[position->jx][position->jyp][position->jz],1,HEATFLUX_INTEGRATION_TAGBASE + mesh->UpXSplitIndex()*mesh->LocalNz*counter + mesh->LocalNz*position->jx + position->jz));
      }
      else {
	mesh->wait(mesh->irecvYOutOutdest(&integral_above[position->jx][position->jyp][position->jz],1,HEATFLUX_INTEGRATION_TAGBASE + (mesh->LocalNx - mesh->UpXSplitIndex())*mesh->LocalNz*counter + mesh->LocalNz*(position->jx-mesh->UpXSplitIndex()) + position->jz));
      }
    }
    do {
      *deltal = mesh->dy[position->jx][position->jy]*sqrt((mesh->g_22[position->jx][position->jy]+mesh->g_22[position->jx][position->jyp])/2.);
      
      interp_coeffs_drive_term = cubic_spline_drive_term.coefficients(position);
      interp_coeffs_lambdaC_inverse = cubic_spline_inverse_lambdaC.coefficients(position);
      
      //calculate the coefficients of drive_term expanded in z rather than l (from the expansion of drive_term in l and l in z [inverted from z in l to third order])
      integral_coeffs[0] = interp_coeffs_drive_term[0];
      integral_coeffs[1] = interp_coeffs_drive_term[1] / *deltal/interp_coeffs_lambdaC_inverse[0];
      integral_coeffs[2] = interp_coeffs_drive_term[2] / pow(*deltal*interp_coeffs_lambdaC_inverse[0],2)
			    - interp_coeffs_drive_term[1] * interp_coeffs_lambdaC_inverse[1]/2./pow(*deltal,2)/pow(interp_coeffs_lambdaC_inverse[0],3);
      integral_coeffs[3] = interp_coeffs_drive_term[3] / pow(*deltal*interp_coeffs_lambdaC_inverse[0],3)
			    - interp_coeffs_drive_term[2] * interp_coeffs_lambdaC_inverse[1]/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4)
			    + interp_coeffs_drive_term[1] * (0.5*pow(interp_coeffs_lambdaC_inverse[1],2)/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],5) - interp_coeffs_lambdaC_inverse[2]/3./pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4));

      // Calculate analytically the integral from jyp to jy of z^n*exp( (zupper-z)/zeta ), firstly allowing for zeta to be large and secondly allowing for it to be small.
      BoutReal exp_delta_over_eigenvalue = exp(dimensionless_length_deltas_above[*position]/eigenvalue);
      if (dimensionless_length_deltas_above[*position]>abs(eigenvalue)) {
	integral_parts[0] = -eigenvalue * (exp_delta_over_eigenvalue-1);
	integral_parts[1] = pow(eigenvalue,2) * (exp_delta_over_eigenvalue-1)
			      - eigenvalue* dimensionless_length_deltas_above[*position] * exp_delta_over_eigenvalue;
	integral_parts[2] = -2*pow(eigenvalue,3) * (exp_delta_over_eigenvalue-1)
			      + 2*pow(eigenvalue,2) * dimensionless_length_deltas_above[*position] * exp_delta_over_eigenvalue
			      - eigenvalue * pow(dimensionless_length_deltas_above[*position],2) * exp_delta_over_eigenvalue;
	integral_parts[3] = 6*pow(eigenvalue,4) * (exp_delta_over_eigenvalue-1)
			      - 6*pow(eigenvalue,3) * dimensionless_length_deltas_above[*position] * exp_delta_over_eigenvalue
			      + 3*pow(eigenvalue,2) * pow(dimensionless_length_deltas_above[*position],2) * exp_delta_over_eigenvalue
			      - eigenvalue * pow(dimensionless_length_deltas_above[*position],3) * exp_delta_over_eigenvalue;
      }
      else {
	BoutReal exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue = exp(dimensionless_length_deltas_above[*position]/eigenvalue/2.) * 2.*sinh(dimensionless_length_deltas_above[*position]/eigenvalue/2.);
	integral_parts[0] = -eigenvalue * exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue;
	integral_parts[1] = pow(eigenvalue,2) * exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      - eigenvalue* dimensionless_length_deltas_above[*position] * exp_delta_over_eigenvalue;
	integral_parts[2] = -2*pow(eigenvalue,3) * exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      + 2*pow(eigenvalue,2) * dimensionless_length_deltas_above[*position] * exp_delta_over_eigenvalue
			      - eigenvalue * pow(dimensionless_length_deltas_above[*position],2) * exp_delta_over_eigenvalue;
	integral_parts[3] = 6*pow(eigenvalue,4) * exp_half_delta_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      - 6*pow(eigenvalue,3) * dimensionless_length_deltas_above[*position] * exp_delta_over_eigenvalue
			      + 3*pow(eigenvalue,2) * pow(dimensionless_length_deltas_above[*position],2) * exp_delta_over_eigenvalue
			      - eigenvalue * pow(dimensionless_length_deltas_above[*position],3) * exp_delta_over_eigenvalue;
      }
 
      //add up the contributions to the integral at jy first from the integral at jy-1 and then from the expansion of the integral between jy-1 and jy
      integral_above[*position] = integral_above[position->jx][position->jyp][position->jz] * exp_delta_over_eigenvalue;
      for (int i=0; i<4; i++) integral_above[*position] += integral_coeffs[i]*integral_parts[i];

    } while (previous_index_y(position));
    
    // Send the value at ystart to the next processor.
    if (position->jx < mesh->DownXSplitIndex()) {
      Timer timer("comms");
      mesh->sendYInIndest(&integral_above[*position],1,HEATFLUX_INTEGRATION_TAGBASE + mesh->DownXSplitIndex()*mesh->LocalNz*counter + mesh->LocalNz*position->jx + position->jz);
    }
    else {
      Timer timer("comms");
      mesh->sendYInOutdest(&integral_above[*position],1,HEATFLUX_INTEGRATION_TAGBASE + (mesh->LocalNx - mesh->DownXSplitIndex())*mesh->LocalNz*counter + mesh->LocalNz*(position->jx-mesh->DownXSplitIndex()) + position->jz);
    }
    
  } while (next_indexperp(position));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HeatFluxIntegration::calculateIntegralBelow_cell_ylow(BoutReal eigenvalue, const Field3D &dimensionless_length_deltas_below, const Field3D &dimensionless_length_deltas_above, CubicSpline &cubic_spline_inverse_lambdaC, CubicSpline &cubic_spline_drive_term, CubicSpline &cubic_spline_gradT, const int &counter) {
  TRACE("HeatFluxIntegration::calculateIntegralBelow()");

  start_index(position);
  do {
    position->jy=mesh->ystart;
    calc_index(position);
    if (!mesh->firstY()) {
      // Set the value at ystart to the value at yend+1 of the previous processor.
      if (position->jx<mesh->DownXSplitIndex()) {
	mesh->wait(mesh->irecvYInIndest(&integral_below[*position],1,HEATFLUX_INTEGRATION_TAGBASE + mesh->DownXSplitIndex()*mesh->LocalNz*counter + mesh->LocalNz*position->jx + position->jz));
      }
      else {
	mesh->wait(mesh->irecvYInOutdest(&integral_below[*position],1,HEATFLUX_INTEGRATION_TAGBASE + (mesh->LocalNx-mesh->DownXSplitIndex())*mesh->LocalNz*counter + mesh->LocalNz*(position->jx-mesh->DownXSplitIndex()) + position->jz));
      }
    }
    do {
      *deltal = mesh->dy[position->jx][position->jy]*sqrt(mesh->g_22[position->jx][position->jy]);
      
      BoutReal deltazbelow = dimensionless_length_deltas_below[*position];
      BoutReal deltazabove = dimensionless_length_deltas_above[*position];
      BoutReal deltaz = dimensionless_length_deltas_below[*position] + dimensionless_length_deltas_above[*position];
      
      interp_coeffs_gradT = cubic_spline_gradT.coefficients(position);
      
      // Contribution to the integral at jy from points up to jy-1
      integral_below[position->jx][position->jyp][position->jz] = integral_below[*position] * exp(deltaz/eigenvalue);
      
      // Calculate the contribution from the lower part of the interval (below CELL_CENTRE)
      // The integration variable ('z-prime') goes from (dimensionless_length_deltas_above[jy-1]) to (dimensionless_length_deltas_above[jy-1]+dimensionless_length_deltas_below[jy]) with the final z-value (that goes into the exponential) being (dimensionless_length_deltas_above[jy-1]+dimensionless_length_deltas_below[jy]+dimensionless_length_deltas_above[jy])
      BoutReal deltazabovefromjym = dimensionless_length_deltas_above[position->jx][position->jym][position->jz];
      interp_coeffs_drive_term = cubic_spline_drive_term.coefficients(position->jx,position->jym,position->jz);
      interp_coeffs_lambdaC_inverse = cubic_spline_inverse_lambdaC.coefficients(position->jx,position->jym,position->jz);
      
      BoutReal integrand_coefficient0 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[0] - 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] - 0.125*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient1 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[1] - interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[0] - 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] - 0.125*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient2 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[2] - 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[1] - interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[0] - 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] - 0.125*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient3 = interp_coeffs_drive_term[0]*interp_coeffs_gradT[3]
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[2] - 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[1] - interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[3]*(interp_coeffs_gradT[0] - 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] - 0.125*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient4 = interp_coeffs_drive_term[1]*interp_coeffs_gradT[3]
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[2] - 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[3]*(interp_coeffs_gradT[1] - interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient5 = interp_coeffs_drive_term[2]*interp_coeffs_gradT[3]
					+ interp_coeffs_drive_term[3]*(interp_coeffs_gradT[2] - 1.5*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient6 = interp_coeffs_drive_term[3]*interp_coeffs_gradT[3];
      
      // Approximate the sixth order polynomial by the least-squares-fit cubic between t=1/2 and t=1
      BoutReal cubic_integrand_coefficient0 = integrand_coefficient0 - 321./1120.*integrand_coefficient4 - 97./112.*integrand_coefficient5 - 2225./1344.*integrand_coefficient6;
      BoutReal cubic_integrand_coefficient1 = integrand_coefficient1 + 45./28.*integrand_coefficient4 + 3065./672.*integrand_coefficient5 + 939./112.*integrand_coefficient6;
      BoutReal cubic_integrand_coefficient2 = integrand_coefficient2 - 93./28.*integrand_coefficient4 - 235./28.*integrand_coefficient5 - 3245./224.*integrand_coefficient6;
      BoutReal cubic_integrand_coefficient3 = integrand_coefficient3 + 3.*integrand_coefficient4 + 205./36.*integrand_coefficient5 + 35./4.*integrand_coefficient6;
      
      //calculate the coefficients of the integrand expanded in z rather than l (from the expansion of drive_term in l and l in z [inverted from z in l to third order])
      integral_coeffs[0] = cubic_integrand_coefficient0;
      integral_coeffs[1] = cubic_integrand_coefficient1 / *deltal/interp_coeffs_lambdaC_inverse[0];
      integral_coeffs[2] = cubic_integrand_coefficient2 / pow(*deltal*interp_coeffs_lambdaC_inverse[0],2)
			    - cubic_integrand_coefficient1 * interp_coeffs_lambdaC_inverse[1]/2./pow(*deltal,2)/pow(interp_coeffs_lambdaC_inverse[0],3);
      integral_coeffs[3] = cubic_integrand_coefficient3 / pow(*deltal*interp_coeffs_lambdaC_inverse[0],3)
			    - cubic_integrand_coefficient2 * interp_coeffs_lambdaC_inverse[1]/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4)
			    + cubic_integrand_coefficient1 * (0.5*pow(interp_coeffs_lambdaC_inverse[1],2)/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],5) - interp_coeffs_lambdaC_inverse[2]/3./pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4));
      
      // Calculate analytically the integral from jy-1/2 to jy of z^n*exp( (deltaz+dimensionless_length_deltas_above[jy-1]-z)/zeta ), firstly allowing for zeta to be large and secondly allowing for it to be small.
      BoutReal exp_deltaabove_over_eigenvalue = exp(deltazabove/eigenvalue);
      if (deltazbelow>abs(eigenvalue)) {
	BoutReal exp_deltabelow_over_eigenvalue = exp(deltazbelow/eigenvalue);
	integral_parts[0] = exp_deltaabove_over_eigenvalue * eigenvalue * (exp_deltabelow_over_eigenvalue-1);
	integral_parts[1] = exp_deltaabove_over_eigenvalue * eigenvalue * ((eigenvalue + deltazabovefromjym) * (exp_deltabelow_over_eigenvalue-1)
									- deltazbelow);
	integral_parts[2] = exp_deltaabove_over_eigenvalue * eigenvalue * ((pow(deltazabovefromjym,2) + 2.*eigenvalue*(deltazabovefromjym+eigenvalue)) * (exp_deltabelow_over_eigenvalue-1)
									- pow(deltazbelow,2)-2.*deltazabovefromjym*deltazbelow - 2.*eigenvalue*deltazbelow);
	integral_parts[3] = exp_deltaabove_over_eigenvalue * eigenvalue * ((pow(deltazabovefromjym,3) + 3.*eigenvalue*pow(deltazabovefromjym,2) + 6.*pow(eigenvalue,2)*deltazabovefromjym + 6.*pow(eigenvalue,3)) * (exp_deltabelow_over_eigenvalue-1)
									- 3.*pow(deltazabovefromjym,2)*deltazbelow - 3.*deltazabovefromjym*pow(deltazbelow,2) - pow(deltazbelow,3) - 3.*eigenvalue*(2.*deltazabovefromjym*deltazbelow + pow(deltazbelow,2)) - 6.*pow(eigenvalue,2)*deltazbelow);
      }
      else {
	BoutReal exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue = exp(deltazbelow/2./eigenvalue) * 2.*sinh(deltazbelow/eigenvalue/2.);
	integral_parts[0] = exp_deltaabove_over_eigenvalue * eigenvalue * exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue;
	integral_parts[1] = exp_deltaabove_over_eigenvalue * eigenvalue * ((eigenvalue + deltazabovefromjym) * exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
									- deltazbelow);
	integral_parts[2] = exp_deltaabove_over_eigenvalue * eigenvalue * ((pow(deltazabovefromjym,2) + 2.*eigenvalue*(deltazabovefromjym+eigenvalue)) * exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
									- pow(deltazbelow,2)-2.*deltazabovefromjym*deltazbelow - 2.*eigenvalue*deltazbelow);
	integral_parts[3] = exp_deltaabove_over_eigenvalue * eigenvalue * ((pow(deltazabovefromjym,3) + 3.*eigenvalue*pow(deltazabovefromjym,2) + 6.*pow(eigenvalue,2)*deltazabovefromjym + 6.*pow(eigenvalue,3)) * exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
									- 3.*pow(deltazabovefromjym,2)*deltazbelow - 3.*deltazabovefromjym*pow(deltazbelow,2) - pow(deltazbelow,3) - 3.*eigenvalue*(2.*deltazabovefromjym*deltazbelow + pow(deltazbelow,2)) - 6.*pow(eigenvalue,2)*deltazbelow);
      }
      
      //add on the contributions to the integral at jy from between jy-1 and jy-1/2
      for (int i=0; i<4; i++) integral_below[position->jx][position->jyp][position->jz] += integral_coeffs[i]*integral_parts[i];
       
      // Calculate the contribution from the upper part of the interval (above CELL_CENTRE)
      // The integration variable ('z-prime') goes from 0 to (dimensionless_length_deltas_above[jy]) with the final z-value (that goes into the exponential) being (dimensionless_length_deltas_above[jy])
      interp_coeffs_drive_term = cubic_spline_drive_term.coefficients(position);
      interp_coeffs_lambdaC_inverse = cubic_spline_inverse_lambdaC.coefficients(position);
      
      integrand_coefficient0 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[0] + 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] + 0.125*interp_coeffs_gradT[3]);
      integrand_coefficient1 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[1] + interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[0] + 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] + 0.125*interp_coeffs_gradT[3]);
      integrand_coefficient2 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[2] + 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[1] + interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[0] + 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] + 0.125*interp_coeffs_gradT[3]);
      integrand_coefficient3 = interp_coeffs_drive_term[0]*interp_coeffs_gradT[3]
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[2] + 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[1] + interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[3]*(interp_coeffs_gradT[0] + 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] + 0.125*interp_coeffs_gradT[3]);
      integrand_coefficient4 = interp_coeffs_drive_term[1]*interp_coeffs_gradT[3]
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[2] + 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[3]*(interp_coeffs_gradT[1] + interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3]);
      integrand_coefficient5 = interp_coeffs_drive_term[2]*interp_coeffs_gradT[3]
					+ interp_coeffs_drive_term[3]*(interp_coeffs_gradT[2] + 1.5*interp_coeffs_gradT[3]);
      integrand_coefficient6 = interp_coeffs_drive_term[3]*interp_coeffs_gradT[3];
      
      // Approximate the sixth order polynomial by the least-squares-fit cubic between t=0 and t=1/2
      cubic_integrand_coefficient0 = integrand_coefficient0 - 1./1120.*integrand_coefficient4 - 1./1008.*integrand_coefficient5 - 1./1344.*integrand_coefficient6;
      cubic_integrand_coefficient1 = integrand_coefficient1 + 1./28.*integrand_coefficient4 + 25./672.*integrand_coefficient5 + 3./112.*integrand_coefficient6;
      cubic_integrand_coefficient2 = integrand_coefficient2 - 9./28.*integrand_coefficient4 - 25./84.*integrand_coefficient5 - 45./224.*integrand_coefficient6;
      cubic_integrand_coefficient3 = integrand_coefficient3 + integrand_coefficient4 + 25./36.*integrand_coefficient5 + 5./12.*integrand_coefficient6;
      
      //calculate the coefficients of the integrand expanded in z rather than l (from the expansion of drive_term in l and l in z [inverted from z in l to third order])
      integral_coeffs[0] = cubic_integrand_coefficient0;
      integral_coeffs[1] = cubic_integrand_coefficient1 / *deltal/interp_coeffs_lambdaC_inverse[0];
      integral_coeffs[2] = cubic_integrand_coefficient2 / pow(*deltal*interp_coeffs_lambdaC_inverse[0],2)
			    - cubic_integrand_coefficient1 * interp_coeffs_lambdaC_inverse[1]/2./pow(*deltal,2)/pow(interp_coeffs_lambdaC_inverse[0],3);
      integral_coeffs[3] = cubic_integrand_coefficient3 / pow(*deltal*interp_coeffs_lambdaC_inverse[0],3)
			    - cubic_integrand_coefficient2 * interp_coeffs_lambdaC_inverse[1]/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4)
			    + cubic_integrand_coefficient1 * (0.5*pow(interp_coeffs_lambdaC_inverse[1],2)/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],5) - interp_coeffs_lambdaC_inverse[2]/3./pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4));

      // Calculate analytically the integral from jy to jy+1/2 of z^n*exp( (deltazabove-z)/zeta ), firstly allowing for zeta to be large and secondly allowing for it to be small.
      if (deltazabove>abs(eigenvalue)) {
	BoutReal exp_deltaabove_over_eigenvalue = exp(deltazabove/eigenvalue);
	integral_parts[0] = eigenvalue * (exp_deltaabove_over_eigenvalue-1);
	integral_parts[1] = eigenvalue * ( eigenvalue*(exp_deltaabove_over_eigenvalue-1)
			      - deltazabove);
	integral_parts[2] = eigenvalue * (2.*pow(eigenvalue,2) * (exp_deltaabove_over_eigenvalue-1)
			      - pow(deltazabove,2) - 2.*deltazabove*eigenvalue);
	integral_parts[3] = eigenvalue*(6.*pow(eigenvalue,3) * (exp_deltaabove_over_eigenvalue-1)
			      - pow(deltazabove,3) - 3.*pow(deltazabove,2)*eigenvalue - 6.*deltazabove*pow(eigenvalue,2));
      }
      else {
	BoutReal exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue = exp(deltazabove/2./eigenvalue) * 2.*sinh(deltazabove/eigenvalue/2.);
	integral_parts[0] = eigenvalue * exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue;
	integral_parts[1] = eigenvalue * ( eigenvalue*exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      - deltazabove);
	integral_parts[2] = eigenvalue * (2.*pow(eigenvalue,2) * exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      - pow(deltazabove,2) - 2.*deltazabove*eigenvalue);
	integral_parts[3] = eigenvalue*(6.*pow(eigenvalue,3) * exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      - pow(deltazabove,3) - 3.*pow(deltazabove,2)*eigenvalue - 6.*deltazabove*pow(eigenvalue,2));
      }
 
      //add on the contributions to the integral at jy from the expansion of the integral between jy-1/2 and jy
      for (int i=0; i<4; i++) integral_below[position->jx][position->jyp][position->jz] += integral_coeffs[i]*integral_parts[i];
       
      position->jy++;
      calc_index(position);
    } while (position->jy < mesh->yend+1);
    
    // Send the value at yend+1 to the next processor.
    if (position->jx < mesh->UpXSplitIndex()) {
      Timer timer("comms");
      mesh->sendYOutIndest(&integral_below[*position],1,HEATFLUX_INTEGRATION_TAGBASE + mesh->UpXSplitIndex()*mesh->LocalNz*counter + mesh->LocalNz*position->jx + position->jz);
    }
    else {
      Timer timer("comms");
      mesh->sendYOutOutdest(&integral_below[*position],1,HEATFLUX_INTEGRATION_TAGBASE + (mesh->LocalNx-mesh->UpXSplitIndex())*mesh->LocalNz*counter + mesh->LocalNz*(position->jx-mesh->UpXSplitIndex()) + position->jz);
    }
    
  } while (next_indexperp(position));
}

void HeatFluxIntegration::calculateIntegralAbove_cell_ylow(BoutReal eigenvalue, const Field3D &dimensionless_length_deltas_below, const Field3D &dimensionless_length_deltas_above, CubicSpline &cubic_spline_inverse_lambdaC, CubicSpline &cubic_spline_drive_term, CubicSpline &cubic_spline_gradT, const int &counter) {
  TRACE("HeatFluxIntegration::calculateIntegralAbove()");

  start_index_lasty(position);
  do {
    position->jy=mesh->yend;
    calc_index(position);
    if (!mesh->lastY()) {
      // Set the value at yend+1 equal to the value at ystart of the previous processor.
      if (position->jx<mesh->UpXSplitIndex()) {
	mesh->wait(mesh->irecvYOutIndest(&integral_above[position->jx][position->jyp][position->jz],1,HEATFLUX_INTEGRATION_TAGBASE + mesh->UpXSplitIndex()*mesh->LocalNz*counter + mesh->LocalNz*position->jx + position->jz));
      }
      else {
	mesh->wait(mesh->irecvYOutOutdest(&integral_above[position->jx][position->jyp][position->jz],1,HEATFLUX_INTEGRATION_TAGBASE + (mesh->LocalNx - mesh->UpXSplitIndex())*mesh->LocalNz*counter + mesh->LocalNz*(position->jx-mesh->UpXSplitIndex()) + position->jz));
      }
    }
    else {
      // Last grid point for CELL_YLOW quantities is mesh->yend on the last y-processor, so start the calculation from mesh->yend-1 if mesh->lastY()
      position->jy=mesh->yend-1;
      calc_index(position);
    }
    do {
      *deltal = mesh->dy[position->jx][position->jy]*sqrt(mesh->g_22[position->jx][position->jy]);
      
      BoutReal deltazbelow = dimensionless_length_deltas_below[*position];
      BoutReal deltazabove = dimensionless_length_deltas_above[*position];
      BoutReal deltaz = dimensionless_length_deltas_below[*position] + dimensionless_length_deltas_above[*position];
      
      interp_coeffs_gradT = cubic_spline_gradT.coefficients(position);
      
      // Contribution to the integral at jy from points up to jy+1
      integral_above[*position] = integral_above[position->jx][position->jyp][position->jz] * exp(deltaz/eigenvalue);
      
      // Calculate the contribution from the lower part of the interval (below CELL_CENTRE)
      // The integration variable ('z-prime') goes from (dimensionless_length_deltas_above[jy-1]+dimensionless_length_deltas_below[jy]) to (dimensionless_length_deltas_above[jy-1]) with the final z-value (that goes into the exponential) being (dimensionless_length_deltas_above[jy-1])
      BoutReal deltazabovefromjym = dimensionless_length_deltas_above[position->jx][position->jym][position->jz];
      interp_coeffs_drive_term = cubic_spline_drive_term.coefficients(position->jx,position->jym,position->jz);
      interp_coeffs_lambdaC_inverse = cubic_spline_inverse_lambdaC.coefficients(position->jx,position->jym,position->jz);
      
      BoutReal integrand_coefficient0 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[0] - 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] - 0.125*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient1 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[1] - interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[0] - 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] - 0.125*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient2 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[2] - 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[1] - interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[0] - 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] - 0.125*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient3 = interp_coeffs_drive_term[0]*interp_coeffs_gradT[3]
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[2] - 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[1] - interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[3]*(interp_coeffs_gradT[0] - 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] - 0.125*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient4 = interp_coeffs_drive_term[1]*interp_coeffs_gradT[3]
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[2] - 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[3]*(interp_coeffs_gradT[1] - interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient5 = interp_coeffs_drive_term[2]*interp_coeffs_gradT[3]
					+ interp_coeffs_drive_term[3]*(interp_coeffs_gradT[2] - 1.5*interp_coeffs_gradT[3]);
      BoutReal integrand_coefficient6 = interp_coeffs_drive_term[3]*interp_coeffs_gradT[3];
      
      // Approximate the sixth order polynomial by the least-squares-fit cubic between t=1/2 and t=1
      BoutReal cubic_integrand_coefficient0 = integrand_coefficient0 - 321./1120.*integrand_coefficient4 - 97./112.*integrand_coefficient5 - 2225./1344.*integrand_coefficient6;
      BoutReal cubic_integrand_coefficient1 = integrand_coefficient1 + 45./28.*integrand_coefficient4 + 3065./672.*integrand_coefficient5 + 939./112.*integrand_coefficient6;
      BoutReal cubic_integrand_coefficient2 = integrand_coefficient2 - 93./28.*integrand_coefficient4 - 235./28.*integrand_coefficient5 - 3245./224.*integrand_coefficient6;
      BoutReal cubic_integrand_coefficient3 = integrand_coefficient3 + 3.*integrand_coefficient4 + 205./36.*integrand_coefficient5 + 35./4.*integrand_coefficient6;
      
      //calculate the coefficients of the integrand expanded in z rather than l (from the expansion of drive_term in l and l in z [inverted from z in l to third order])
      integral_coeffs[0] = cubic_integrand_coefficient0;
      integral_coeffs[1] = cubic_integrand_coefficient1 / *deltal/interp_coeffs_lambdaC_inverse[0];
      integral_coeffs[2] = cubic_integrand_coefficient2 / pow(*deltal*interp_coeffs_lambdaC_inverse[0],2)
			    - cubic_integrand_coefficient1 * interp_coeffs_lambdaC_inverse[1]/2./pow(*deltal,2)/pow(interp_coeffs_lambdaC_inverse[0],3);
      integral_coeffs[3] = cubic_integrand_coefficient3 / pow(*deltal*interp_coeffs_lambdaC_inverse[0],3)
			    - cubic_integrand_coefficient2 * interp_coeffs_lambdaC_inverse[1]/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4)
			    + cubic_integrand_coefficient1 * (0.5*pow(interp_coeffs_lambdaC_inverse[1],2)/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],5) - interp_coeffs_lambdaC_inverse[2]/3./pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4));
      
      // Calculate analytically the integral from jy+1/2 to jy of z^n*exp( -(deltazabovefromjym-z)/zeta ), firstly allowing for zeta to be large and secondly allowing for it to be small.
      if (deltazbelow>abs(eigenvalue)) {
	BoutReal exp_deltabelow_over_eigenvalue = exp(deltazbelow/eigenvalue);
	integral_parts[0] = -eigenvalue * (exp_deltabelow_over_eigenvalue-1);
	integral_parts[1] = -eigenvalue * ((-eigenvalue +  deltazabovefromjym) * (exp_deltabelow_over_eigenvalue-1)
					    + deltazbelow*exp_deltabelow_over_eigenvalue);
	integral_parts[2] = -eigenvalue * ((pow(deltazabovefromjym,2) - 2.*eigenvalue*(deltazabovefromjym-eigenvalue)) * (exp_deltabelow_over_eigenvalue-1)
					    + (pow(deltazbelow,2)+2.*deltazabovefromjym*deltazbelow - 2.*eigenvalue*deltazbelow)*exp_deltabelow_over_eigenvalue);
	integral_parts[3] = -eigenvalue * ((pow(deltazabovefromjym,3) - 3.*eigenvalue*pow(deltazabovefromjym,2) + 6.*pow(eigenvalue,2)*deltazabovefromjym - 6.*pow(eigenvalue,3)) * (exp_deltabelow_over_eigenvalue-1)
					    + (3.*pow(deltazabovefromjym,2)*deltazbelow + 3.*deltazabovefromjym*pow(deltazbelow,2) + pow(deltazbelow,3) - 3.*eigenvalue*(2.*deltazabovefromjym*deltazbelow + pow(deltazbelow,2)) + 6.*pow(eigenvalue,2)*deltazbelow)*exp_deltabelow_over_eigenvalue);
      }
      else {
	BoutReal exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue = exp(deltazbelow/2./eigenvalue) * 2.*sinh(deltazbelow/eigenvalue/2.);
	integral_parts[0] = -eigenvalue * exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue;
	integral_parts[1] = -eigenvalue * ((-eigenvalue + deltazbelow + deltazabovefromjym) * exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
					    + deltazbelow);
	integral_parts[2] = -eigenvalue * ((pow(deltazabovefromjym+deltazbelow,2) - 2.*eigenvalue*(deltazabovefromjym+deltazbelow-eigenvalue)) * exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
					    + pow(deltazbelow,2)+2.*deltazabovefromjym*deltazbelow - 2.*eigenvalue*deltazbelow);
	integral_parts[3] = -eigenvalue * ((pow(deltazabovefromjym + deltazbelow,3) - 3.*eigenvalue*pow(deltazabovefromjym + deltazbelow,2) + 6.*pow(eigenvalue,2)*(deltazabovefromjym + deltazbelow) - 6.*pow(eigenvalue,3)) * exp_half_deltabelow_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
					    + 3.*pow(deltazabovefromjym,2)*deltazbelow + 3.*deltazabovefromjym*pow(deltazbelow,2) + pow(deltazbelow,3) - 3.*eigenvalue*(2.*deltazabovefromjym*deltazbelow + pow(deltazbelow,2)) + 6.*pow(eigenvalue,2)*deltazbelow);
      }
      
      //add on the contributions to the integral at jy from between jy and jy+1/2
      for (int i=0; i<4; i++) integral_above[*position] += integral_coeffs[i]*integral_parts[i];
      
      
      // Calculate the contribution from the upper part of the interval (above CELL_CENTRE)
      // The integration variable ('z-prime') goes from (dimensionless_length_deltas_above[jy]) to 0 with the final z-value (that goes into the exponential) being (-dimensionless_length_deltas_below[jy])
      interp_coeffs_drive_term = cubic_spline_drive_term.coefficients(position);
      interp_coeffs_lambdaC_inverse = cubic_spline_inverse_lambdaC.coefficients(position);
      
      integrand_coefficient0 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[0] + 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] + 0.125*interp_coeffs_gradT[3]);
      integrand_coefficient1 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[1] + interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[0] + 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] + 0.125*interp_coeffs_gradT[3]);
      integrand_coefficient2 = interp_coeffs_drive_term[0]*(interp_coeffs_gradT[2] + 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[1] + interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[0] + 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] + 0.125*interp_coeffs_gradT[3]);
      integrand_coefficient3 = interp_coeffs_drive_term[0]*interp_coeffs_gradT[3]
					  + interp_coeffs_drive_term[1]*(interp_coeffs_gradT[2] + 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[1] + interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[3]*(interp_coeffs_gradT[0] + 0.5*interp_coeffs_gradT[1] + 0.25*interp_coeffs_gradT[2] + 0.125*interp_coeffs_gradT[3]);
      integrand_coefficient4 = interp_coeffs_drive_term[1]*interp_coeffs_gradT[3]
					  + interp_coeffs_drive_term[2]*(interp_coeffs_gradT[2] + 1.5*interp_coeffs_gradT[3])
					  + interp_coeffs_drive_term[3]*(interp_coeffs_gradT[1] + interp_coeffs_gradT[2] + 0.75*interp_coeffs_gradT[3]);
      integrand_coefficient5 = interp_coeffs_drive_term[2]*interp_coeffs_gradT[3]
					+ interp_coeffs_drive_term[3]*(interp_coeffs_gradT[2] + 1.5*interp_coeffs_gradT[3]);
      integrand_coefficient6 = interp_coeffs_drive_term[3]*interp_coeffs_gradT[3];
      
      // Approximate the sixth order polynomial by the least-squares-fit cubic between t=0 and t=1/2
      cubic_integrand_coefficient0 = integrand_coefficient0 - 1./1120.*integrand_coefficient4 - 1./1008.*integrand_coefficient5 - 1./1344.*integrand_coefficient6;
      cubic_integrand_coefficient1 = integrand_coefficient1 + 1./28.*integrand_coefficient4 + 25./672.*integrand_coefficient5 + 3./112.*integrand_coefficient6;
      cubic_integrand_coefficient2 = integrand_coefficient2 - 9./28.*integrand_coefficient4 - 25./84.*integrand_coefficient5 - 45./224.*integrand_coefficient6;
      cubic_integrand_coefficient3 = integrand_coefficient3 + integrand_coefficient4 + 25./36.*integrand_coefficient5 + 5./12.*integrand_coefficient6;
      
      //calculate the coefficients of the integrand expanded in z rather than l (from the expansion of drive_term in l and l in z [inverted from z in l to third order])
      integral_coeffs[0] = cubic_integrand_coefficient0;
      integral_coeffs[1] = cubic_integrand_coefficient1 / *deltal/interp_coeffs_lambdaC_inverse[0];
      integral_coeffs[2] = cubic_integrand_coefficient2 / pow(*deltal*interp_coeffs_lambdaC_inverse[0],2)
			    - cubic_integrand_coefficient1 * interp_coeffs_lambdaC_inverse[1]/2./pow(*deltal,2)/pow(interp_coeffs_lambdaC_inverse[0],3);
      integral_coeffs[3] = cubic_integrand_coefficient3 / pow(*deltal*interp_coeffs_lambdaC_inverse[0],3)
			    - cubic_integrand_coefficient2 * interp_coeffs_lambdaC_inverse[1]/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4)
			    + cubic_integrand_coefficient1 * (0.5*pow(interp_coeffs_lambdaC_inverse[1],2)/pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],5) - interp_coeffs_lambdaC_inverse[2]/3./pow(*deltal,3)/pow(interp_coeffs_lambdaC_inverse[0],4));

      // Calculate analytically the integral from jy+1 to jy+1/2 of z^n*exp( -(-deltazbelow-z)/zeta ), firstly allowing for zeta to be large and secondly allowing for it to be small.
      BoutReal exp_deltabelow_over_eigenvalue = exp(deltazbelow/eigenvalue);
      if (deltazabove>abs(eigenvalue)) {
	BoutReal exp_deltaabove_over_eigenvalue = exp(deltazabove/eigenvalue);
	integral_parts[0] = -exp_deltabelow_over_eigenvalue * eigenvalue * (exp_deltaabove_over_eigenvalue-1);
	integral_parts[1] = -exp_deltabelow_over_eigenvalue * eigenvalue * (-eigenvalue * (exp_deltaabove_over_eigenvalue-1)
			      + deltazabove*exp_deltaabove_over_eigenvalue);
	integral_parts[2] = -exp_deltabelow_over_eigenvalue * eigenvalue * (2.*pow(eigenvalue,2) * (exp_deltaabove_over_eigenvalue-1)
			      + (pow(deltazabove,2) - 2.*eigenvalue*deltazabove)*exp_deltaabove_over_eigenvalue);
	integral_parts[3] = -exp_deltabelow_over_eigenvalue * eigenvalue * (-6.*pow(eigenvalue,3) * (exp_deltaabove_over_eigenvalue-1)
			      + (pow(deltazabove,3) - 3.*eigenvalue*pow(deltazabove,2) + 6.*pow(eigenvalue,2)*deltazabove)*exp_deltaabove_over_eigenvalue);
      }
      else {
	BoutReal exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue = exp(deltazabove/2./eigenvalue) * 2.*sinh(deltazabove/eigenvalue/2.);
	integral_parts[0] = -exp_deltabelow_over_eigenvalue * eigenvalue * exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue;
	integral_parts[1] = -exp_deltabelow_over_eigenvalue * eigenvalue * ((-eigenvalue + deltazabove) * exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      + deltazabove);
	integral_parts[2] = -exp_deltabelow_over_eigenvalue * eigenvalue * ((pow(deltazabove,2) - 2.*eigenvalue*deltazabove + 2.*pow(eigenvalue,2)) * exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      + pow(deltazabove,2) - 2.*eigenvalue*deltazabove);
	integral_parts[3] = -exp_deltabelow_over_eigenvalue * eigenvalue * ((pow(deltazabove,3) - 3.*eigenvalue*pow(deltazabove,2) + 6.*pow(eigenvalue,2)*deltazabove - 6.*pow(eigenvalue,3)) * exp_half_deltaabove_over_eigenvalue_times_two_sinh_half_delta_over_eigenvalue
			      + pow(deltazabove,3) - 3.*eigenvalue*pow(deltazabove,2) + 6.*pow(eigenvalue,2)*deltazabove);
      }
      
      //add on the contributions to the integral at jy from the expansion of the integral between jy+1 and jy+1/2
      for (int i=0; i<4; i++) integral_above[*position] += integral_coeffs[i]*integral_parts[i];

    } while (previous_index_y(position));
    
    // Send the value at ystart to the next processor
    if (position->jx < mesh->DownXSplitIndex()) {
      Timer timer("comms");
      mesh->sendYInIndest(&integral_above[*position],1,HEATFLUX_INTEGRATION_TAGBASE + mesh->DownXSplitIndex()*mesh->LocalNz*counter + mesh->LocalNz*position->jx + position->jz);
    }
    else {
      Timer timer("comms");
      mesh->sendYInOutdest(&integral_above[*position],1,HEATFLUX_INTEGRATION_TAGBASE + (mesh->LocalNx - mesh->DownXSplitIndex())*mesh->LocalNz*counter + mesh->LocalNz*(position->jx-mesh->DownXSplitIndex()) + position->jz);
    }
    
  } while (next_indexperp(position));
}
