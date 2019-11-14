/************************************************************************
 * Inversion of parallel derivatives
 * 
 * Inverts a matrix of the form 
 *
 * (A + B * Grad2_par2) x = r
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

#include <invert_parderiv.hxx>
#include "impls/cyclic/cyclic.hxx"

const Field2D InvertPar::solve(const Field2D &f) {
  Field3D var(f);
  
  var = solve(var);
  return DC(var);
}
