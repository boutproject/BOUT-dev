/*!
 * \file conservative_flux_div.hxx
 *
 * \brief Conservative flux divergence operator for finite-volume methods
 *
 * \author BOUT++ Team
 *
 **************************************************************************
 * Copyright 2024 BOUT++ Contributors
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
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

#pragma once
#ifndef BOUT_CONSERVATIVE_FLUX_DIV_H
#define BOUT_CONSERVATIVE_FLUX_DIV_H

#include "bout/facefield3d.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/coordinates.hxx"
#include "bout/region.hxx"

namespace FV {

/*!
 * \brief Compute conservative divergence of face-centered fluxes
 * 
 * This function implements a finite-volume divergence operator that
 * ensures exact conservation. Given fluxes on cell faces, it computes
 * the divergence using:
 * 
 * div(F) = 1/V * sum(F·n·dA)
 * 
 * where the sum is over all cell faces, V is the cell volume,
 * n is the outward normal, and dA is the face area.
 * 
 * \param F Face-centered flux field with components on cell faces
 * \param bndry_flux If true, use fluxes at domain boundaries (default: false)
 * \return Cell-centered divergence field
 * 
 * \note This operator ensures exact conservation: the sum of div(F)*V
 *       over all cells equals the net flux through domain boundaries
 */
Field3D ConservativeFluxDiv(const FaceField3D& F, bool bndry_flux = false);

/*!
 * \brief Compute conservative divergence with specified boundary flux
 * 
 * This variant allows specifying the flux values at domain boundaries
 * rather than using zero flux (default) or extrapolation.
 * 
 * \param F Face-centered flux field
 * \param bndry_flux_x Flux at X boundaries (optional)
 * \param bndry_flux_y Flux at Y boundaries (optional) 
 * \param bndry_flux_z Flux at Z boundaries (optional)
 * \return Cell-centered divergence field
 */
Field3D ConservativeFluxDiv(const FaceField3D& F,
                           const Field3D* bndry_flux_x,
                           const Field3D* bndry_flux_y,
                           const Field3D* bndry_flux_z);

} // namespace FV

#endif // BOUT_CONSERVATIVE_FLUX_DIV_H