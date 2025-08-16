/*!
 * \file mesh_facefield_comm.cxx
 *
 * \brief Implementation notes for face-centered field communication
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

#include "bout/mesh_facefield_comm.hxx"
#include "bout/globals.hxx"

/*!
 * \brief Implementation notes for face field communication
 * 
 * The communication of face-centered fields relies on BOUT++'s existing
 * support for staggered grids. When a Field3D has a location of CELL_XLOW,
 * CELL_YLOW, or CELL_ZLOW, the communication routines automatically handle
 * the half-cell offset.
 * 
 * Key points:
 * 
 * 1. StaggerGrids Support:
 *    - If Mesh::StaggerGrids is true, BOUT++ already handles staggered locations
 *    - The communicate() method checks field locations and adjusts accordingly
 *    - No additional modification to core communication is needed
 * 
 * 2. Face Field Specifics:
 *    - X-face fields (CELL_XLOW) have nx+1 points in X direction
 *    - Y-face fields (CELL_YLOW) have ny+1 points in Y direction  
 *    - Z-face fields (CELL_ZLOW) have nz+1 points in Z direction
 *    - Ghost cells are exchanged to maintain flux continuity
 * 
 * 3. Conservative Properties:
 *    - Face values at processor boundaries are shared between domains
 *    - This ensures conservation when computing divergence
 *    - The "owner" of a boundary face is determined by BOUT++ conventions
 * 
 * 4. Enabling Face Communication:
 *    - Set mesh:staggergrids = true in input file, OR
 *    - The conservative divergence operator will enable it internally
 *    - Face fields work even without global StaggerGrids enabled
 * 
 * The implementation in mesh_facefield_comm.hxx provides convenience
 * methods that simply call the standard communicate() on each component.
 * The actual staggered communication logic is in the core BOUT++ mesh.
 */

// Note: No additional implementation needed here as all methods are inline
// in the header file. This file exists primarily for documentation and
// to maintain consistency with BOUT++ structure.