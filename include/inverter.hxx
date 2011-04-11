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

#include "fieldperp.hxx"
#include "dcomplex.hxx"

#include <vector>

class Inverter {
 public:
  Inverter();
  ~Inverter();
  
  /// Main function to start the inversion
  int solve(const FieldPerp &b, FieldPerp &x,  
	    int flags=0, 
	    int restart=10, int itmax=100,
	    BoutReal tol=1.e-7);
  
  /// Solve for a 3D field, one slice at a time
  int solve(const Field3D &b, Field3D &x, 
	    int flags=0, 
	    int restart=10, int itmax=100,
	    BoutReal tol=1.e-7);
  
  /// User must implement this function
  virtual const FieldPerp function(const FieldPerp &x) = 0;
  
  void setBoundaryFlags(int flags) {bndry_flags = flags; }
  void setXCoupling(bool xcouple) {parallel = xcouple & nxgt1; }
  
 protected:
  /// Apply a boundary condition to a FieldPerp variable
  void applyBoundary(FieldPerp &f, int flags);

 private:
  
  void calcBoundary(dcomplex **cdata, int n, BoutReal *h, int flags);

  bool nxgt1; // True if NXPE > 1
  bool parallel; // True if need to communicate
  int bndry_flags; // Boundary condition flags
  
  void A(BoutReal *b, BoutReal *x); ///< Calculates b = Ax
  
  // GMRES solver code
  
  BoutReal norm_vector(BoutReal *b, int n);
  BoutReal dot_product(BoutReal *a, BoutReal *b, int n);
  void Update(BoutReal *x, int it, BoutReal **h, BoutReal *s, BoutReal *y, BoutReal **v, int n);
  void GeneratePlaneRotation(BoutReal dx, BoutReal dy, BoutReal *cs, BoutReal *sn);
  void ApplyPlaneRotation(BoutReal *dx, BoutReal *dy, BoutReal cs, BoutReal sn);
  int gmres_solve(BoutReal *b, BoutReal *x, int n, int m, int itmax, BoutReal tol, int &iterations, BoutReal &residual);
};

// Variable Flags


// Solve Flags

const int INVERT_BOUNDARY = 1; // Set boundary condition for each variable

#endif // __INVERTER_H__

