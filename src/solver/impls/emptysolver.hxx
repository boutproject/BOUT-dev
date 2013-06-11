/**************************************************************************
 * Empty solver class for throwing errors
 * 
 * NOTE: Only one solver can currently be compiled in
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

class EmptySolver;

#ifndef __EMPTY_SOLVER_H__
#define __EMPTY_SOLVER_H__

#include <bout/solver.hxx>
#include <boutexception.hxx>

class EmptySolver : public Solver {
public:
  EmptySolver(Options *opt = NULL) {throw BoutException("Solver not enabled!");}
  
  int run() {return 0;}
  BoutReal run(BoutReal tout, int &ncalls, BoutReal &rhstime) {return 0;}

};

#endif // __EMPTY_SOLVER_H__

