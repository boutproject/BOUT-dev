/*******************************************************************************
 * Inversion routine using GMRES
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
 *******************************************************************************/

#ifndef __INVERT_GMRES_H__
#define __INVERT_GMRES_H__

#include "fieldperp.hxx"

typedef int (*opfunc2D) (FieldPerp &b, FieldPerp &x, void *data);

int iter_solve(FieldPerp &b, FieldPerp &x, opfunc2D A, void *extra);
int iter_solve_bndry(FieldPerp &b, FieldPerp &x, opfunc2D A, int flags, void *extra);

int iter_solve(Field3D &b, Field3D &x, opfunc2D A, void *extra);
int iter_solve_bndry(Field3D &b, Field3D &x, opfunc2D A, int flags, void *extra);

#endif // __INVERT_GMRES_H__

