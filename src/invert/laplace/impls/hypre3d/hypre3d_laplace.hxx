/**************************************************************************
 * 3D Laplace inversion using Hypre solvers with algebraic multigrid
 *  preconditioning
 *
 * Equation solved is: \f$d\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + ex\nabla_x x + ez\nabla_z x + a x = b\f$
 *
 **************************************************************************
 * Copyright 2021 J. Omotani, C. MacMackin
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
class LaplaceHypre3d;

#include "bout/build_config.hxx"

#ifndef __LAPLACE_HYPRE3D_H__
#define __LAPLACE_HYPRE3D_H__

#if BOUT_HAS_HYPRE

#include <globals.hxx>
#include <output.hxx>
#include <options.hxx>
#include <invert_laplace.hxx>
#include <boutexception.hxx>
#include <bout/operatorstencil.hxx>
#include <bout/hypre_interface.hxx>

class LaplaceHypre3d;

namespace {
RegisterLaplace<LaplaceHypre3d> registerlaplacehypre3d(LAPLACE_HYPRE3D);
}

class LaplaceHypre3d : public Laplacian {
public:
  LaplaceHypre3d(Options *opt = nullptr, const CELL_LOC loc = CELL_CENTRE, Mesh *mesh_in = nullptr);
  ~LaplaceHypre3d() override;

  void setCoefA(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    A = val;
    updateRequired = true;
  }
  void setCoefC(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1 = val;
    C2 = val;
    updateRequired = issetC = true;
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
    updateRequired = issetC = true;
  }
  void setCoefD(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    D = val;
    updateRequired = issetD = true;
  }
  void setCoefEx(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Ex = val;
    updateRequired = issetE = true;
  }
  void setCoefEz(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Ez = val;
    updateRequired = issetE = true;
  }

  void setCoefA(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    A = val;
    updateRequired = true;
  }
  void setCoefC(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1 = val;
    C2 = val;
    updateRequired = issetC = true;
  }
  void setCoefC1(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1 = val;
    updateRequired = issetC = true;
  }
  void setCoefC2(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C2 = val;
    updateRequired = issetC = true;
  }
  void setCoefD(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    D = val;
    updateRequired = issetD = true;
  }
  void setCoefEx(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Ex = val;
    updateRequired = issetE = true;
  }
  void setCoefEz(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Ez = val;
    updateRequired = issetE = true;
  }

  
  // Return a reference to the matrix objects representing the Laplace
  // operator. These will be (re)construct if necessary.
  bout::HypreMatrix<Field3D>& getMatrix3D();
  IndexerPtr<Field3D> getIndexer() { return indexer; }

  virtual Field2D solve(const Field2D &b) override;

  virtual Field3D solve(const Field3D &b) override {
    Field3D zero = zeroFrom(b);
    return solve(b, zero);
  }
  virtual Field3D solve(const Field3D &b_in, const Field3D &x0) override;


  virtual FieldPerp solve(const FieldPerp& UNUSED(b)) override {
    throw BoutException("LaplaceHypre3d cannot solve for FieldPerp");
  }

//private:
  public:

  // (Re)compute the values of the matrix representing the Laplacian operator
  void updateMatrix3D();

  // Set up a stencil describing the structure of the operator. This
  // will be used to preallocate the correct ammount of memory in the
  // matrix.
  //
  // The interpolation done to get along-field values in y makes this
  // tricky. For now we will just assume that the footprint of cells
  // used for interpolation is the same everywhere.
  static OperatorStencil<Ind3D> getStencil(Mesh* localmesh,
                                           const RangeIterator &lowerYBound,
                                           const RangeIterator &upperYBound);
  
  /* Ex and Ez
   * Additional 1st derivative terms to allow for solution field to be
   * components of a vector
   *
   * See LaplaceHypre3d::Coeffs for details an potential pit falls
   */
  Field3D A, C1, C2, D, Ex, Ez;
  bool issetD = false;
  bool issetC = false;
  bool issetE = false;
  bool updateRequired = true;
  int lastflag;               // The flag used to construct the matrix
  int lower_boundary_flags;
  int upper_boundary_flags;

  int meshx, meshz, size, localN; // Mesh sizes, total size, no of points on this processor

  Options *opts;              // Laplace Section Options Object
  std::string ksptype; ///< KSP solver type
  std::string pctype;  ///< Preconditioner type

  // Convergence Parameters. Solution is considered converged if |r_k| < max( rtol * |b| , atol )
  // where r_k = b - Ax_k. The solution is considered diverged if |r_k| > dtol * |b|.
  BoutReal rtol, atol, dtol;
  int maxits; // Maximum number of iterations in solver.

  RangeIterator lowerY, upperY;

  IndexerPtr<Field3D> indexer;
  bout::HypreMatrix<Field3D> operator3D;
  bout::HypreVector<Field3D> solution;
  bout::HypreVector<Field3D> rhs;
  bout::HypreSystem<Field3D> linearSystem;

  bool use_precon;  // Switch for preconditioning
  bool rightprec;   // Right preconditioning

  // These are the implemented flags
  static constexpr int implemented_flags = INVERT_START_NEW,
    implemented_boundary_flags = INVERT_AC_GRAD + INVERT_SET + INVERT_RHS;
};

#endif //BOUT_HAS_HYPRE

#endif //__LAPLACE_HYPRE3D_H__
