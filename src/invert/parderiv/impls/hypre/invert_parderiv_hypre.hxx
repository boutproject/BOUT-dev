/************************************************************************
 * Inversion of parallel derivatives
 * 
 * Inverts a matrix of the form 
 *
 * A + B * Grad2_par2
 * 
 * - Each flux surface can be solved independently
 * - By taking FFT in Z direction, toroidal modes can be
 *   solved independently
 * - Result is a set of complex band-diagonal matrices to invert
 * 
 * Author: Ben Dudson, University of York, Oct 2011
 * 
 * Known issues:
 * ------------
 * 
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

#ifndef __INV_PAR_HYPRE_H__
#define __INV_PAR_HYPRE_H__

#include <invert_parderiv.hxx>
#include <dcomplex.hxx>
#include <mesh.hxx>
#include <field3d.hxx>
#include <field2d.hxx>

#include <HYPRE_struct_ls.h>

class InvertParHypre : public InvertPar {
public:
  InvertParHypre(Mesh *mesh);
  ~InvertParHypre();
  const Field3D solve(const Field3D &f);
  const Field3D solve(const Field3D &f, const Field3D &start);
  
  void setCoefA(const Field2D &f) {A = f;}
  void setCoefB(const Field2D &f) {B = f;}
  
private:
  Field2D A, B;
  
  Mesh *m;
  dcomplex **crhs; ///< Array for storing FFT of RHS
  dcomplex **cresult; ///< Result array
  
  HYPRE_StructGrid *grid; // Array of grids, one per surface
  HYPRE_StructStencil stencil; // Only one stencil needed
  HYPRE_StructMatrix **mat; // One per surface, per Z mode
  
};

#endif // __INV_PAR_HYPRE_H__
