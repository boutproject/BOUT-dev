/**************************************************************************
 * Perpendicular Laplacian inversion.
 *                           Using PETSc Solvers
 *
 * Equation solved is: d*\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + a x = b
 *
 **************************************************************************
 * Copyright 2013 J. Buchanan, J. Omotani
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

#ifndef BOUT_HAS_PETSC

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
    delete [] ksptype;
    delete [] pctype;
  }

  void setCoefA(const Field2D &val) { A = val; /*Acoefchanged = true;*/ if(pcsolve) pcsolve->setCoefA(val); }
  void setCoefC(const Field2D &val) { C1 = val; C2 = val; issetC = true; /*coefchanged = true;*/ if(pcsolve) pcsolve->setCoefC(val);  }
  void setCoefC1(const Field2D &val) { C1 = val; issetC = true;}
  void setCoefC2(const Field2D &val) { C2 = val; issetC = true;}
  void setCoefD(const Field2D &val) { D = val; issetD = true; /*coefchanged = true;*/ if(pcsolve) pcsolve->setCoefD(val); }
  void setCoefEx(const Field2D &val) { Ex = val; issetE = true; /*coefchanged = true;*/ if(pcsolve) pcsolve->setCoefEx(val); }
  void setCoefEz(const Field2D &val) { Ez = val; issetE = true; /*coefchanged = true;*/ if(pcsolve) pcsolve->setCoefEz(val); }

  void setCoefA(const Field3D &val) { A = val; /*Acoefchanged = true;*/ if(pcsolve) pcsolve->setCoefA(val);}
  void setCoefC(const Field3D &val) { C1 = val; C2 = val; issetC = true; /*coefchanged = true;*/ if(pcsolve) pcsolve->setCoefC(val); }
  void setCoefC1(const Field3D &val) { C1 = val; issetC = true;}
  void setCoefC2(const Field3D &val) { C2 = val; issetC = true;}
  void setCoefD(const Field3D &val) { D = val; issetD = true; /*coefchanged = true;*/ if(pcsolve) pcsolve->setCoefD(val); }
  void setCoefEx(const Field3D &val) { Ex = val; issetE = true; /*coefchanged = true;*/ if(pcsolve) pcsolve->setCoefEx(val); }
  void setCoefEz(const Field3D &val) { Ez = val; issetE = true; /*coefchanged = true;*/ if(pcsolve) pcsolve->setCoefEz(val); }

  const FieldPerp solve(const FieldPerp &b);
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0);

  int precon(Vec x, Vec y); ///< Preconditioner function

private:
  void Element(int i, int x, int z, int xshift, int zshift, PetscScalar ele, Mat &MatA );
  void Coeffs( int x, int y, int z, BoutReal &A1, BoutReal &A2, BoutReal &A3, BoutReal &A4, BoutReal &A5 );

  Field3D A, C1, C2, D, Ex, Ez;
// Metrics are not constant in y-direction, so matrix always changes as you loop over the grid
// Hence using coefchanged switch to avoid recomputing the mmatrix is not a useful thing to do (unless maybe in a cylindrical machine, but not worth implementing just for that)
//   bool coefchanged;           // Set to true when C, D, Ex or Ez coefficients are changed
//   bool Acoefchanged;       // Set to true when A coefficient is changed
  bool issetD;
  bool issetC;
  bool issetE;
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
  PCType pctype;          // Preconditioner type

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

  bool use_precon;  // Switch for preconditioning
  bool rightprec;   // Right preconditioning
  Laplacian *pcsolve; // Laplacian solver for preconditioning

  void vecToField(Vec x, FieldPerp &f);        // Copy a vector into a fieldperp
  void fieldToVec(const FieldPerp &f, Vec x);  // Copy a fieldperp into a vector

#ifdef CHECK
  int implemented_flags;
  int implemented_boundary_flags;
#endif
};

#endif //BOUT_HAS_PETSC_DEV

#endif //__PETSC_LAPLACE_H__
