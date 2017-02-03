/**************************************************************************
 * Class for non-linear Laplacian inversion
 *
 * Equation solved is: \nabla^2_\perp x + \nabla_perp c\cdot\nabla_\perp x + a x = b
 *
 * i.e. the same as invert_laplace.hxx
 * 
 * The difference is that here a and c can be 3D variables, not 2D
 * i.e. this solver does not need to make the Boussinesq approximation for
 * vorticity equation inversion.
 *
 * Optionally uses invert_laplace method as a preconditioner
 * 
 * Changelog: 
 *
 * 2010-05-04 Ben Dudson <bd512@york.ac.uk>
 *    * Initial version
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

#ifndef __INVERT_LAP_GMRES_H__
#define __INVERT_LAP_GMRES_H__

#include "inverter.hxx"

class LaplaceGMRES : public Inverter {
 public:
  /// Main solver function. Pass NULL to omit terms
  const Field3D invert(const Field3D &b, const Field3D &start, int inv_flags, bool precon=true, Field3D *a=NULL, Field3D *c=NULL);
  
  /// Implement the function to be inverted
  const FieldPerp function(const FieldPerp &x);
 private:
  int flags;
  
  bool enable_a, enable_c; // Terms enabled
  Field3D a3d, c3d;

  Field2D a2d, c2d;  // DC components (for preconditioner)
  Field2D *aptr, *cptr; // Pointers to the 2D variables (for passing to preconditioner)
};

#endif // __INVERT_LAP_GMRES_H__


