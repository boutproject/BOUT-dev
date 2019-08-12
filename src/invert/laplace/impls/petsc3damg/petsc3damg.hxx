/**************************************************************************
 * 3D Laplace inversion using PETSc solvers with algebraic multigrid
 *  preconditioning
 *
 * Equation solved is: \f$d\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + ex\nabla_x x + ez\nabla_z x + a x = b\f$
 *
 **************************************************************************
 * Copyright 2019 J. Buchanan, J. Omotani, C. MacMackin
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
class LaplacePetsc3dAmg;

#ifndef __PETSC_LAPLACE_3DAMG_H__
#define __PETSC_LAPLACE_3DAMG_H__

#ifndef BOUT_HAS_PETSC

#include <boutexception.hxx>
#include <invert_laplace.hxx>

class LaplacePetsc3dAmg : public Laplacian {
public:
  LaplacePetsc3dAmg(Options *UNUSED(opt) = nullptr, const CELL_LOC UNUSED(loc) = CELL_CENTRE, Mesh *UNUSED(mesh_in) = nullptr) {
    throw BoutException("No PETSc solver available");
  }

  using Laplacian::setCoefA;
  void setCoefA(const Field2D &UNUSED(val)) override {}
  using Laplacian::setCoefC;
  void setCoefC(const Field2D &UNUSED(val)) override {}
  using Laplacian::setCoefD;
  void setCoefD(const Field2D &UNUSED(val)) override {}
  using Laplacian::setCoefEx;
  void setCoefEx(const Field2D &UNUSED(val)) override {}
  using Laplacian::setCoefEz;
  void setCoefEz(const Field2D &UNUSED(val)) override {}

  using Laplacian::solve;
  const FieldPerp solve(const FieldPerp &UNUSED(b)) override {throw BoutException("PETSc not available");}
};

#else

#include <globals.hxx>
#include <output.hxx>
#include <petscksp.h>
#include <options.hxx>
#include <invert_laplace.hxx>
#include <bout/petsclib.hxx>
#include <boutexception.hxx>

class LaplacePetsc3dAmg : public Laplacian {
public:
  LaplacePetsc3dAmg(Options *opt = nullptr, const CELL_LOC loc = CELL_CENTRE, Mesh *mesh_in = nullptr);
  ~LaplacePetsc3dAmg(){};

  void setCoefA(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    A = val;
  }
  void setCoefC(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1 = val;
    C2 = val;
  }
  void setCoefC1(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1 = val;
    issetC = true;
  }
  void setCoefC2(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C2 = val;
    issetC = true;
  }
  void setCoefD(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    D = val;
    issetD = true;
  }
  void setCoefEx(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Ex = val;
    issetE = true;
  }
  void setCoefEz(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Ez = val;
    issetE = true;
  }

  void setCoefA(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    A = val;
  }
  void setCoefC(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1 = val;
    C2 = val;
    issetC = true;
  }
  void setCoefC1(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1 = val;
    issetC = true;
  }
  void setCoefC2(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C2 = val;
    issetC = true;
  }
  void setCoefD(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    D = val;
    issetD = true;
  }
  void setCoefEx(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Ex = val;
    issetE = true;
  }
  void setCoefEz(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Ez = val;
    issetE = true;
  }

  const Field3D solve(const Field3D &b) override {
    Field3D zero(b.getMesh());
    zero = 0.;
    return solve(b, zero);
  }
  const Field3D solve(const Field3D &b_in, const Field3D &x0) override;

  const FieldPerp solve(const FieldPerp& UNUSED(b)) override {
    throw BoutException("LaplacePetsc3DAmg cannot solve for FieldPerp");
  }

private:
 
  /* Ex and Ez
   * Additional 1st derivative terms to allow for solution field to be
   * components of a vector
   *
   * See LaplacePetsc::Coeffs for details an potential pit falls
   */
  Field3D A, C1, C2, D, Ex, Ez;
  bool issetD;
  bool issetC;
  bool issetE;
  int lastflag;               // The flag used to construct the matrix

  int meshx, meshz, size, localN; // Mesh sizes, total size, no of points on this processor
  MPI_Comm comm;

  Options *opts;              // Laplace Section Options Object
  std::string ksptype; ///< KSP solver type
  std::string pctype;  ///< Preconditioner type

  // Convergence Parameters. Solution is considered converged if |r_k| < max( rtol * |b| , atol )
  // where r_k = b - Ax_k. The solution is considered diverged if |r_k| > dtol * |b|.
  BoutReal rtol, atol, dtol;
  int maxits; // Maximum number of iterations in solver.
  bool direct; //Use direct LU solver if true.
  bool fourth_order;

  PetscLib lib;

  bool use_precon;  // Switch for preconditioning
  bool rightprec;   // Right preconditioning

  #if CHECK > 0
    int implemented_flags;
    int implemented_boundary_flags;
  #endif
};

#endif //BOUT_HAS_PETSC

#endif //__PETSC_LAPLACE_3DAMG_H__
