/************************************************************************
 * Inversion of parallel derivatives
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
 * Author: Ben Dudson, University of York, June 2009
 * 
 * Known issues:
 * ------------
 *
 * - Assumes all points are in the core: no handling of SOL / PF regions
 * - This algorithm can only use MXSUB processors, so if NYPE > MXSUB
 *   (i.e. Np > Nx) then efficiency will fall.
 * - When higher-order interpolation is used at the twist-shift location,
 *   errors can become large
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

#ifndef __INV_PAR_SIMPLE_H__
#define __INV_PAR_SIMPLE_H__

#include <invert_parderiv.hxx>

class InvertParSimple : public InvertPar {
public:
  const Field3D solve(const Field3D &f);
  
  void setCoefA(const Field2D &f) {A = f;}
  void setCoefB(const Field2D &f) {B = f;}
  
private:
  Field2D A, B;
  
  void cyclicSolve(int ysize, int xpos, BoutReal *data, 
                   BoutReal *result, bool coef3d = false);
  
};

#endif // __INV_PAR_SIMPLE_H__
