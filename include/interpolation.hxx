/**************************************************************************
 * Functions to interpolate between cell locations (e.g. lower Y and centred)
 * 
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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
 **************************************************************************/

#ifndef __INTERP_H__
#define __INTERP_H__

#include "field3d.hxx"
#include "bout_types.hxx"

/// Interpolate to a give cell location
const Field3D interp_to(const Field3D &var, CELL_LOC loc);
const Field2D interp_to(const Field2D &var, CELL_LOC loc);

/// Print out the cell location (for debugging)
void printLocation(const Field3D &var);

const char* strLocation(CELL_LOC loc);


/// Interpolate a field onto a perturbed set of points
const Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z);

const Field3D interpolate(const Field2D &f, const Field3D &delta_x, const Field3D &delta_z);
const Field3D interpolate(const Field2D &f, const Field3D &delta_x);

#endif // __INTERP_H__
