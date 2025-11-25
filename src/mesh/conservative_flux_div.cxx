/*!
 * \file conservative_flux_div.cxx
 *
 * \brief Implementation of conservative flux divergence operator
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

#include "bout/conservative_flux_div.hxx"
#include "bout/assert.hxx"
#include "bout/msg_stack.hxx"
#include "bout/utils.hxx"

namespace FV {

Field3D ConservativeFluxDiv(const FaceField3D& F, bool bndry_flux) {
  TRACE("FV::ConservativeFluxDiv(FaceField3D)");
  
  return ConservativeFluxDiv(F, 
                           bndry_flux ? &F.x() : nullptr,
                           bndry_flux ? &F.y() : nullptr,
                           bndry_flux ? &F.z() : nullptr);
}

Field3D ConservativeFluxDiv(const FaceField3D& F,
                           const Field3D* bndry_flux_x,
                           const Field3D* bndry_flux_y,
                           const Field3D* bndry_flux_z) {
  TRACE("FV::ConservativeFluxDiv(FaceField3D, boundaries)");
  
  Mesh* mesh = F.getMesh();
  ASSERT1(mesh != nullptr);
  
  // Get the result field
  Field3D result{zeroFrom(F.x())};
  result.setLocation(CELL_CENTRE);
  
  // Get coordinates
  Coordinates* coord = mesh->getCoordinates();
  ASSERT1(coord != nullptr);
  
  // Get mesh spacing
  // Note: These may be Field2D or BoutReal depending on grid
  const auto& dx = coord->dx;
  const auto& dy = coord->dy;
  const auto& dz = coord->dz;
  
  // Jacobian
  const auto& J = coord->J;
  
  // Ensure staggered grids are enabled for face fields
  if (!mesh->StaggerGrids) {
    // Enable staggered grids for face field operations
    // This is safe as it only affects how locations are handled
    output_warn.write("ConservativeFluxDiv: Enabling StaggerGrids for face field operations\n");
    mesh->StaggerGrids = true;
  }
  
  // Loop over interior points
  // For now, use simple loops to avoid potential issues with BOUT_FOR
  for (int ix = mesh->xstart; ix <= mesh->xend; ix++) {
    for (int iy = mesh->ystart; iy <= mesh->yend; iy++) {
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
    
    // X-direction flux difference
    // F.x() is located at CELL_XLOW, so F.x(i,j,k) is flux at left face of cell (i,j,k)
    // F.x(i+1,j,k) is flux at right face of cell (i,j,k)
    BoutReal flux_diff_x = F.x()(ix+1, iy, iz) - F.x()(ix, iy, iz);
    
    // Y-direction flux difference
    BoutReal flux_diff_y = F.y()(ix, iy+1, iz) - F.y()(ix, iy, iz);
    
    // Z-direction flux difference
    BoutReal flux_diff_z = F.z()(ix, iy, iz+1) - F.z()(ix, iy, iz);
    
    // Get metric factors at cell center
    BoutReal J_c = J(ix, iy, iz);
    BoutReal dx_c = dx(ix, iy, iz);
    BoutReal dy_c = dy(ix, iy, iz);
    BoutReal dz_c = dz(ix, iy, iz);
    
    // Compute divergence
    // div(F) = 1/J * (dF_x/dx + dF_y/dy + dF_z/dz)
    // Note: Face fluxes should already include appropriate metric factors
    result(ix, iy, iz) = (flux_diff_x / dx_c + 
                         flux_diff_y / dy_c + 
                         flux_diff_z / dz_c) / J_c;
      }
    }
  }
  
  // Handle boundary conditions
  if (mesh->firstX()) {
    // X lower boundary
    for (int ix = 0; ix < mesh->xstart; ix++) {
      for (int iy = mesh->ystart; iy <= mesh->yend; iy++) {
        for (int iz = 0; iz < mesh->LocalNz; iz++) {
          if (bndry_flux_x) {
            // Use specified boundary flux
            BoutReal flux_diff_x = F.x()(ix+1, iy, iz) - (*bndry_flux_x)(ix, iy, iz);
            BoutReal flux_diff_y = F.y()(ix, iy+1, iz) - F.y()(ix, iy, iz);
            BoutReal flux_diff_z = F.z()(ix, iy, iz+1) - F.z()(ix, iy, iz);
            
            result(ix, iy, iz) = (flux_diff_x / dx(ix, iy, iz) + 
                                 flux_diff_y / dy(ix, iy, iz) + 
                                 flux_diff_z / dz(ix, iy, iz)) / J(ix, iy, iz);
          } else {
            // Zero flux boundary condition
            result(ix, iy, iz) = 0.0;
          }
        }
      }
    }
  }
  
  if (mesh->lastX()) {
    // X upper boundary
    for (int ix = mesh->xend+1; ix < mesh->LocalNx; ix++) {
      for (int iy = mesh->ystart; iy <= mesh->yend; iy++) {
        for (int iz = 0; iz < mesh->LocalNz; iz++) {
          if (bndry_flux_x) {
            // Use specified boundary flux
            BoutReal flux_diff_x = (*bndry_flux_x)(ix, iy, iz) - F.x()(ix, iy, iz);
            BoutReal flux_diff_y = F.y()(ix, iy+1, iz) - F.y()(ix, iy, iz);
            BoutReal flux_diff_z = F.z()(ix, iy, iz+1) - F.z()(ix, iy, iz);
            
            result(ix, iy, iz) = (flux_diff_x / dx(ix, iy, iz) + 
                                 flux_diff_y / dy(ix, iy, iz) + 
                                 flux_diff_z / dz(ix, iy, iz)) / J(ix, iy, iz);
          } else {
            // Zero flux boundary condition
            result(ix, iy, iz) = 0.0;
          }
        }
      }
    }
  }
  
  // Similar handling for Y boundaries if needed
  // Note: Y boundaries are typically handled differently in BOUT++ due to field alignment
  
  return result;
}

} // namespace FV