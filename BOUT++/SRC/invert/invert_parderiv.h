/************************************************************************
 * Inversion of parallel derivatives
 * Intended for use in preconditioner for reduced MHD
 * 
 * Inverts a matrix of the form 
 *
 * A + B * Grad2_par2
 * 
 * Stages:
 * - Problem trivially parallel in X, so gather all data for fixed X onto
 *   a single processor. Split MXSUB locations between NYPE processors
 * - Use nearest neighbour for twist-shift location.
 *   This splits the problem into one or more (cyclic) tridiagonal problems
 *   (number depends on the q value)
 * - Solve each of these tridiagonal systems O(Nz*Ny)
 * - Scatter data back 
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
 ************************************************************************/

#ifndef __INV_PAR_H__
#define __INV_PAR_H__

#include "field3d.h"
#include "field2d.h"

namespace invpar {
  const Field3D invert_parderiv(const Field2D &A, const Field2D &B, const Field3D &r);
  const Field3D invert_parderiv(real val, const Field2D &B, const Field3D &r);
  const Field3D invert_parderiv(const Field2D &A, real val, const Field3D &r);
  const Field3D invert_parderiv(real val, real val2, const Field3D &r);
}

using invpar::invert_parderiv;


#endif // __INV_PAR_H__

