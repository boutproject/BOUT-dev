/**************************************************************************
 * Perpendicular Laplacian inversion. Parallel code using FFT
 * and tridiagonal solver.
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

class LaplaceParallelTri;

#ifndef __PARALLEL_TRI_H__
#define __PARALLEL_TRI_H__

#include <invert_laplace.hxx>
#include <dcomplex.hxx>
#include <options.hxx>

class LaplaceParallelTri : public Laplacian {
public:
  LaplaceParallelTri(Options *opt = nullptr, const CELL_LOC loc = CELL_CENTRE, Mesh *mesh_in = nullptr);
  ~LaplaceParallelTri(){};

  using Laplacian::setCoefA;
  void setCoefA(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    A = val;
  }
  using Laplacian::setCoefC;
  void setCoefC(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C = val;
  }
  using Laplacian::setCoefD;
  void setCoefD(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    D = val;
  }
  using Laplacian::setCoefEx;
  void setCoefEx(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceParallelTri does not have Ex coefficient");
  }
  using Laplacian::setCoefEz;
  void setCoefEz(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceParallelTri does not have Ez coefficient");
  }

  using Laplacian::solve;
  FieldPerp solve(const FieldPerp &b) override;
  FieldPerp solve(const FieldPerp &b, const FieldPerp &x0) override;
  //FieldPerp solve(const FieldPerp &b, const FieldPerp &x0, const FieldPerp &b0 = 0.0);

  BoutReal getMeanIterations() const { return ipt_mean_its; }
  void resetMeanIterations() { ipt_mean_its = 0; }

  void ensure_stability(const int jy, const int kz, Array<dcomplex> &r,
      Tensor<dcomplex> &lowerGuardVector, Tensor<dcomplex> &upperGuardVector,
      bool &lowerUnstable, bool &upperUnstable);
  void swapHaloInteriorLower(Array<dcomplex> &x);
  void swapHaloInteriorUpper(Array<dcomplex> &x);
  void check_diagonal_dominance(const Array<dcomplex> &a, const Array<dcomplex> &b,
      const Array<dcomplex> &c, const int ncx, const int jy, const int kz);

  void resetSolver();

private:
  // The coefficents in
  // D*grad_perp^2(x) + (1/C)*(grad_perp(C))*grad_perp(x) + A*x = b
  Field2D A, C, D;

  BoutReal omega;
  BoutReal Borig;
  Field3D Bvals;
  Matrix<bool> first_call;
  // Save previous x in Fourier space
  Tensor<dcomplex> x0saved;

  /// Solver tolerances
  BoutReal rtol, atol;

  /// Maximum number of iterations
  int maxits;

  /// Kludge factor
  BoutReal B;

  /// Mean number of iterations taken by the solver
  BoutReal ipt_mean_its;

  /// Counter for the number of times the solver has been called
  int ncalls;

};

#endif // __PARALLEL_TRI_H__
