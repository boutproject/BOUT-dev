/**************************************************************************
 * Class for non-linear inversion problems
 *
 * Uses GMRES to solve problems of form F(x)=b, either on a single processor
 * or in parallel
 *
 * Changelog: 
 *
 * 2007-10 Ben Dudson <bd512@york.ac.uk>
 *    * Initial version. Not working yet.
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

class Inverter;

#ifndef __INVERTER_H__
#define __INVERTER_H__

#include "fieldperp.h"

#include <vector>

class Inverter {
 public:
  Inverter();
  ~Inverter();
  
  /// Main function to start the inversion
  int solve(const FieldPerp &b, FieldPerp &x,  
	    int flags=0, 
	    int restart=10, int itmax=100,
	    real tol=1.e-7);
  
  /// Solve for a 3D field, one slice at a time
  int solve(const Field3D &b, Field3D &x, 
	    int flags=0, 
	    int restart=10, int itmax=100,
	    real tol=1.e-7);
  
  /// User must implement this function
  virtual const FieldPerp function(const FieldPerp &x) = 0;
  
 private:
  
  void A(real *b, real *x); ///< Calculates b = Ax
  
  // GMRES solver code
  
  real norm_vector(real *b, int n);
  real dot_product(real *a, real *b, int n);
  void Update(real *x, int it, real **h, real *s, real *y, real **v, int n);
  void GeneratePlaneRotation(real dx, real dy, real *cs, real *sn);
  void ApplyPlaneRotation(real *dx, real *dy, real cs, real sn);
  int gmres_solve(real *b, real *x, int n, int m, int itmax, real tol, int &iterations, real &residual);
};

// Variable Flags


// Solve Flags

const int INVERT_BOUNDARY = 1; // Set boundary condition for each variable

#endif // __INVERTER_H__

