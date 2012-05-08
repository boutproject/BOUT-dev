/**************************************************************
 * Smoothing operators
 *
 **************************************************************
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
 **************************************************************/

#ifndef __SMOOTHING_H__
#define __SMOOTHING_H__

#include "field3d.hxx"

/// Smooth in X using simple 1-2-1 filter
const Field3D smooth_x(const Field3D &f, bool BoutRealspace = true);

/// Smooth in Y using 1-2-1 filter
const Field3D smooth_y(const Field3D &f);

/// Smooth using a stencil in X and Y
const Field3D smoothXY(const Field3D &f);

/// Average over Y
const Field2D averageY(const Field2D &f);
const Field3D averageY(const Field3D &f);

/// Non-linear filter to remove grid-scale oscillations
const Field3D nl_filter_x(const Field3D &f, BoutReal w=1.0);
const Field3D nl_filter_y(const Field3D &f, BoutReal w=1.0);
const Field3D nl_filter_z(const Field3D &f, BoutReal w=1.0);
const Field3D nl_filter(const Field3D &f, BoutReal w=1.0);

#endif // __SMOOTHING_H__
