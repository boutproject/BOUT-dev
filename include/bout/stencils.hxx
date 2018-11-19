/*!************************************************************************
 * \file stencils.hxx
 * 
 * Sets stencils for differencing
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

#ifndef __STENCILS_H__
#define __STENCILS_H__

#include "bout/bout_types.hxx"

/// Defines a set of values in 1D in the neighbourhood of an index
/// Used for calculating derivatives
struct stencil {
  /// stencil 2 each side of the centre -- in effect means M?G > 2 is not supported
  BoutReal mm = BoutNaN, m = BoutNaN, c = BoutNaN, p = BoutNaN, pp = BoutNaN; 
};

#endif /* __STENCILS_H__ */
