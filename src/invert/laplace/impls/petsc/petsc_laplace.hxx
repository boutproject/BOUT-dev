/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using PETSc Solvers
 *
 **************************************************************************
 * Copyright 2013 J. Buchan
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
class LaplacePetsc;

#ifndef __PETSC_LAPLACE_H__
#define __PETSC_LAPLACE_H__

#ifndef BOUT_HAS_PETSC_3_3

#include <boutexception.hxx>
#include <invert_laplace.hxx>

class LaplacePetsc : public Laplacian {
public:
  LaplacePetsc(Options *opt = NULL) { throw BoutException("No PETSc solver available"); }
  
  void setCoefA(const Field2D &val) {}
  void setCoefB(const Field2D &val) {}
  void setCoefC(const Field2D &val) {}
  void setCoefD(const Field2D &val) {}
  void setCoefEx(const Field2D &val) {}
  void setCoefEz(const Field2D &val) {}
  
  const FieldPerp solve(const FieldPerp &b) {throw BoutException("PETSc not available");}
};

#else

#include <globals.hxx>
#include <output.hxx>
#include <petscksp.h>
#include <options.hxx>
#include <invert_laplace.hxx>
#include <bout/petsclib.hxx>
#include <boutexception.hxx>

class LaplacePetsc : public Laplacian {
public:
  LaplacePetsc(Options *opt = NULL);
  ~LaplacePetsc() {
    KSPDestroy( &ksp ); 
    VecDestroy( &xs );  
    VecDestroy( &bs ); 
    MatDestroy( &MatA );
  }
  
  void setCoefA(const Field2D &val) { A = val; /*Acoefchanged = true;*/}
  void setCoefC(const Field2D &val) { C = val; /*coefchanged = true;*/}
  void setCoefD(const Field2D &val) { D = val; /*coefchanged = true;*/}
  void setCoefEx(const Field2D &val) { Ex = val; /*coefchanged = true;*/}
  void setCoefEz(const Field2D &val) { Ez = val; /*coefchanged = true;*/}

  void setCoefA(const Field3D &val) { A = val; /*Acoefchanged = true;*/}
  void setCoefC(const Field3D &val) { C = val; /*coefchanged = true;*/}
  void setCoefD(const Field3D &val) { D = val; /*coefchanged = true;*/}
  void setCoefEx(const Field3D &val) { Ex = val; /*coefchanged = true;*/}
  void setCoefEz(const Field3D &val) { Ez = val; /*coefchanged = true;*/}
  
  const FieldPerp solve(const FieldPerp &b);
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0);

private:
  void Element(int i, int x, int z, int xshift, int zshift, PetscScalar ele, Mat &MatA );
  void Coeffs( int x, int y, int z, BoutReal &A1, BoutReal &A2, BoutReal &A3, BoutReal &A4, BoutReal &A5 );
  
  Field3D A, C, D, Ex, Ez;
// Metrics are not constant in y-direction, so matrix always changes as you loop over the grid
// Hence using coefchanged switch to avoid recomputing the mmatrix is not a useful thing to do (unless maybe in a cylindrical machine, but not worth implementing just for that)
//   bool coefchanged;           // Set to true when C, D, Ex or Ez coefficients are changed
//   bool Acoefchanged;	      // Set to true when A coefficient is changed
  int lastflag;               // The flag used to construct the matrix

  FieldPerp sol;              // solution Field
  
  // Istart is the first row of MatA owned by the process, Iend is 1 greater than the last row.
  int Istart, Iend; 

  int meshx, meshz, size, localN;
  MPI_Comm comm;
  Mat MatA;
  Vec xs, bs;                 // Solution and RHS vectors
  KSP ksp;
  
  Options *opts;              // Laplace Section Options Object
  KSPType ksptype;            // Solver Type;
  PCType pctype;	      // Preconditioner type

  // Values specific to particular solvers
  BoutReal richardson_damping_factor;  
  BoutReal chebyshev_max, chebyshev_min;
  int gmres_max_steps;

  // Convergence Parameters. Solution is considered converged if |r_k| < max( rtol * |b| , atol )
  // where r_k = b - Ax_k. The solution is considered diverged if |r_k| > dtol * |b|.
  BoutReal rtol, atol, dtol;
  int maxits; // Maximum number of iterations in solver.
  bool direct; //Use direct LU solver if true.
  bool fourth_order;

  PetscLib lib;
  
  #ifdef CHECK
    int implemented_flags;
  #endif
};

#endif //BOUT_HAS_PETSC_DEV

#endif //__PETSC_LAPLACE_H__
