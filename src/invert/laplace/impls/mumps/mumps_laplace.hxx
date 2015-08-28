/**************************************************************************
 * Perpendicular Laplacian inversion.
 *                           Using MUMPS Solver
 *
 **************************************************************************
 * Copyright 2013 Copyright 2013 J. Omotani (based on petsc_laplace (C) J. Buchan)
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
class LaplaceMumps;

#ifndef __MUMPS_LAPLACE_H__
#define __MUMPS_LAPLACE_H__

#include <invert_laplace.hxx>

#ifndef BOUT_HAS_MUMPS

#include <boutexception.hxx>

class LaplaceMumps : public Laplacian {
public:
  LaplaceMumps(Options *opt = NULL) { throw BoutException("Mumps library not available"); }

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
#include <options.hxx>
#include <boutexception.hxx>
#include "dmumps_c.h"

#define MUMPS_JOB_INIT -1
#define MUMPS_JOB_END -2
#define MUMPS_JOB_ANALYSIS 1
#define MUMPS_JOB_FACTORIZATION 2
#define MUMPS_JOB_SOLUTION 3
#define MUMPS_JOB_ANALYSIS_AND_FACTORIZATION 4 // combines job 1 and job 2
#define MUMPS_JOB_BOTH 5 // combines job 2 and job 3
#define MUMPS_JOB_ALL 6 // combines job 1, job 2 and job 3

class LaplaceMumps : public Laplacian {
public:
  LaplaceMumps(Options *opt = NULL);
  ~LaplaceMumps() {
    mumps_struc.job = -2;
    dmumps_c(&mumps_struc);
    delete [] mumps_struc.irn_loc;
    delete [] mumps_struc.jcn_loc;
    delete [] mumps_struc.a_loc;
    delete [] mumps_struc.isol_loc;
    delete [] rhs;
//     delete [] rhs_slice;
    delete [] localrhs;
    delete [] localrhs_size_array;
    delete [] rhs_positions;
  }

  void setCoefA(const Field2D &val) { A = val; }
  void setCoefC(const Field2D &val) { C1 = val; C2 = val; issetC = true; }
  void setCoefC1(const Field2D &val) { C1 = val; issetC = true; }
  void setCoefC2(const Field2D &val) { C2 = val; issetC = true; }
  void setCoefD(const Field2D &val) { D = val; issetD = true; }
  void setCoefEx(const Field2D &val) { Ex = val; issetE = true; }
  void setCoefEz(const Field2D &val) { Ez = val; issetE = true; }

  void setCoefA(const Field3D &val) { A = val; }
  void setCoefC(const Field3D &val) { C1 = val; C2 = val; issetC = true; }
  void setCoefC1(const Field3D &val) { C1 = val; issetC = true; }
  void setCoefC2(const Field3D &val) { C2 = val; issetC = true; }
  void setCoefD(const Field3D &val) { D = val; issetD = true; }
  void setCoefEx(const Field3D &val) { Ex = val; issetE = true; }
  void setCoefEz(const Field3D &val) { Ez = val; issetE = true; }

  void setFlags(int f) {throw BoutException("May not change the value of flags during run in LaplaceMumps as it might change the number of non-zero matrix elements: flags may only be set in the options file.");}

  const FieldPerp solve(const FieldPerp &b);
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0);
//   const Field3D solve(const Field3D &b);
//   const Field3D solve(const Field3D &b, const Field3D &x0);

private:
  void solve(BoutReal* rhs, int y);
  void Coeffs( int x, int y, int z, BoutReal &A1, BoutReal &A2, BoutReal &A3, BoutReal &A4, BoutReal &A5 );

  Field3D A, C1, C2, D, Ex, Ez;

  bool issetD;
  bool issetC;
  bool issetE;
//   int repeat_analysis; // Repeat analysis step after this many iterations
//   int iteration_count; // Use this to count the number of iterations since last analysis

  BoutReal* rhs; // Array to collect rhs field onto host processor
//   BoutReal* rhs_slice; // Array to pass xz-slice of rhs to solve
  BoutReal* localrhs;
  int localrhssize;
  int* localrhs_size_array;
  int* rhs_positions;
  FieldPerp sol;              // solution Field

  // Istart is the first row of MatA owned by the process, Iend is 1 greater than the last row.
  int Istart, Iend;

  int meshx, meshz, size, localN, nxguards;
  MPI_Comm comm;

  Options *opts;              // Laplace Section Options Object
  bool fourth_order;

#ifdef CHECK
  int implemented_flags;
  int implemented_boundary_flags;
#endif

  DMUMPS_STRUC_C mumps_struc;
};

#endif //BOUT_HAS_MUMPS

#endif //MUMPS_LAPLACE_H
