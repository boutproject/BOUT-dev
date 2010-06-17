/* 
 * This switches between different solver header files based 
 * on command-line defines. 
 *
 * Need a better solution later...
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
 */

#ifndef __SOLVER_H__
#define __SOLVER_H__

#ifdef IDA

#include "ida_solver.h"

#else
// Using a CVODE solver
#ifdef CVODE

#include "sundials_solver.h"

#else

#ifdef PETSC

#include "petsc_solver.h"

#else

#include "cvode_solver.h"

#endif // CVODE

#endif // PETSC

#endif // IDA

#endif // __SOLVER_H__
