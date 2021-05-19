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

class LaplacePCR;

#ifndef __PCR_H__
#define __PCR_H__

#include <dcomplex.hxx>
#include <invert_laplace.hxx>
#include <options.hxx>
#include <utils.hxx>

namespace {
RegisterLaplace<LaplacePCR> registerlaplacepcr(LAPLACE_PCR);
}

class LaplacePCR : public Laplacian {
public:
  LaplacePCR(Options* opt = nullptr, const CELL_LOC loc = CELL_CENTRE,
             Mesh* mesh_in = nullptr);
  ~LaplacePCR() = default;

  using Laplacian::setCoefA;
  void setCoefA(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Acoef = val;
  }
  using Laplacian::setCoefC;
  void setCoefC(const Field2D &val) override {
    setCoefC1(val);
    setCoefC2(val);
  }
  using Laplacian::setCoefC1;
  void setCoefC1(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1coef = val;
  }
  using Laplacian::setCoefC2;
  void setCoefC2(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C2coef = val;
  }
  using Laplacian::setCoefD;
  void setCoefD(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Dcoef = val;
  }
  using Laplacian::setCoefEx;
  void setCoefEx(const Field2D& UNUSED(val)) override {
    throw BoutException("LaplaceParallelTriMG does not have Ex coefficient");
  }
  using Laplacian::setCoefEz;
  void setCoefEz(const Field2D& UNUSED(val)) override {
    throw BoutException("LaplaceParallelTriMG does not have Ez coefficient");
  }

  using Laplacian::solve;
  FieldPerp solve(const FieldPerp& b) override { return solve(b, b); }
  FieldPerp solve(const FieldPerp& b, const FieldPerp& x0) override;

  Field3D solve(const Field3D &b) override {return solve(b,b);}
  Field3D solve(const Field3D &b, const Field3D &x0) override;

        void setup(int n, int np_world, int rank_world);
        void cr_solver    (double *a_mpi, double *b_mpi, double *c_mpi, double *r_mpi, double *x_mpi);
        //void cr_pcr_solver(double *a_mpi, double *b_mpi, double *c_mpi, double *r_mpi, double *x_mpi);
        //void cr_pcr_solver(Tensor<dcomplex> &a_mpi, Tensor<dcomplex> &b_mpi, Tensor<dcomplex> &c_mpi, Matrix<dcomplex> &r_mpi, Matrix<dcomplex> &x_mpi, int jy);
        void cr_pcr_solver(Matrix<dcomplex> &a_mpi, Matrix<dcomplex> &b_mpi, Matrix<dcomplex> &c_mpi, Matrix<dcomplex> &r_mpi, Matrix<dcomplex> &x_mpi);
        void Thomas_pcr_solver(double *a_mpi, double *b_mpi, double *c_mpi, double *r_mpi, double *x_mpi);
        void verify_solution(double *a_ver, double *b_ver, double *c_ver, double *r_ver, double *x_sol);


private:
        Field2D Acoef, C1coef, C2coef, Dcoef;
        Matrix<dcomplex> bcmplx, xcmplx;

        /// Number of rows per MPI process and should be 2^n.
        int n_mpi;
        
        /// Number of MPI process and should be also 2^m.
        int nprocs;

        /// MPI process ID
        int myrank;

        /// Local private pointer for coefficient maxtix a
        Matrix<dcomplex> a, aa;
        /// Local private pointer for coefficient maxtix b
        Matrix<dcomplex> b, bb;
        /// Local private pointer for coefficient maxtix c
        Matrix<dcomplex> c, cc;
        /// Local private pointer for RHS vector r
        Matrix<dcomplex> r;
        /// Local private pointer for solution vector x
        Matrix<dcomplex> x;

        void cr_forward_multiple_row(Matrix<dcomplex> &a,Matrix<dcomplex> &b,Matrix<dcomplex> &c,Matrix<dcomplex> &r);
        void cr_backward_multiple_row(Matrix<dcomplex> &a,Matrix<dcomplex> &b,Matrix<dcomplex> &c,Matrix<dcomplex> &r,Matrix<dcomplex> &x);
        void apply_boundary_conditions(const Matrix<dcomplex> &a,const Matrix<dcomplex> &b,const Matrix<dcomplex> &c,const Matrix<dcomplex> &r,Matrix<dcomplex> &x);
        void eliminate_boundary_rows(const Matrix<dcomplex> &a,Matrix<dcomplex> &b,const Matrix<dcomplex> &c,Matrix<dcomplex> &r);
        void cr_forward_single_row();
        void cr_backward_single_row();
        void pcr_forward_single_row(Matrix<dcomplex> &a,Matrix<dcomplex> &b,Matrix<dcomplex> &c,Matrix<dcomplex> &r,Matrix<dcomplex> &x);
        void verify_solution(const Matrix<dcomplex> &a_ver, const Matrix<dcomplex> &b_ver, const Matrix<dcomplex> &c_ver, const Matrix<dcomplex> &r_ver, const Matrix<dcomplex> &x_sol);

        //void pThomas_forward_multiple_row();
        //void pcr_double_row_substitution();

  /// The coefficents in
  /// $D*grad_perp^2(x) + (1/C)*(grad_perp(C))*grad_perp(x) + A*x = b$
///  Field2D A, C, D;

  /// Number of unfiltered Fourier modes
  int nmode;

  /// Number of systems to solve = number of unfiltered Fourier modes times number of y points
  int nsys;

  /// Number of local x, y points
  int ncx, ny;

  /// Current y index
///  int jy;

  /// Lower-, on- and upper-diagonal terms of the operator matrix
  Tensor<dcomplex> avec, bvec, cvec;

  /// Counter for the number of times the solver has been called
  int ncalls{0};

  /// Neighbouring processors in the in and out directions
  int proc_in, proc_out;

  /// This processor's unique ID
  int myproc;

  /// Shorthand for localmesh->NXPE
  int nproc;

  /// First and last interior points xstart, xend
  int xs, xe;

  bool isGlobalFlagSet(int flag) const {
    return (global_flags & flag) != 0;
  }
  bool isInnerBoundaryFlagSet(int flag) const {
    return (inner_boundary_flags & flag) != 0;
  }
  bool isOuterBoundaryFlagSet(int flag) const {
    return (outer_boundary_flags & flag) != 0;
  }

  bool dst;
};

#endif // __PCR_H__
