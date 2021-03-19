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

class Laplace1DMG;

#ifndef __1DMG_H__
#define __1DMG_H__

#include <dcomplex.hxx>
#include <invert_laplace.hxx>
#include <options.hxx>
#include <utils.hxx>

namespace {
RegisterLaplace<Laplace1DMG> registerlaplace1dmg(LAPLACE_1DMG);
}

class Laplace1DMG : public Laplacian {
public:
  Laplace1DMG(Options *opt = nullptr, const CELL_LOC loc = CELL_CENTRE,
	  Mesh* mesh_in = nullptr);
  ~Laplace1DMG() = default;

  friend class Level;

  using Laplacian::setCoefA;
  void setCoefA(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    A = val;
  }
  using Laplacian::setCoefC;
  void setCoefC(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C = val;
  }
  using Laplacian::setCoefD;
  void setCoefD(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    D = val;
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

  BoutReal getMeanIterations() const { return ipt_mean_its; }
  void resetMeanIterations() { ipt_mean_its = 0; }

  void get_initial_guess(const int jy, const int kz, Matrix<dcomplex> &r,
      Tensor<dcomplex> &lowerGuardVector, Tensor<dcomplex> &upperGuardVector,
      Matrix<dcomplex> &xk1d);

  void resetSolver();

  bool all(const Array<bool>);

  class Level {

  public:
    // Constructor for zeroth level; needs to modify \p lap
    Level(Laplace1DMG& lap);
    // Constructor for all other levels; uses the previous level \p lup
    Level(const Laplace1DMG& lap, const Level& lup, std::size_t current_level);

    Matrix<dcomplex> xloc;
    Matrix<dcomplex> residual;
    Tensor<dcomplex> ar, br, cr, brinv;
    Matrix<dcomplex> rr;

    /// Unique ID
    int myproc;
    /// In-neighbour
    int proc_in;
    /// Out-neighbour
    int proc_out;
    /// Whether this processor is included in this grid level's calculation
    bool included;
    /// Colouring of processor for Gauss-Seidel
    bool red, black;
    /// Current grid level being solved
    std::size_t current_level;

    // indexing to remove branches from tight loops
    int index_start;
    int index_end;

    // Save some proc properties from the level above - this allows us
    // to NOT pass the level above as an argument in some functions

    /// Whether this processor is involved in the calculation on the grid one level more
    /// refined
    bool included_up;
    /// This processor's neighbours on the level above
    int proc_in_up, proc_out_up;

    void calculate_residual(const Laplace1DMG& lap);
    void calculate_total_residual(const Laplace1DMG& lap, Array<BoutReal>& error,
                                  Array<bool>& converged);
    void coarsen(const Laplace1DMG& lap, const Matrix<dcomplex>& fine_residual);
    void gauss_seidel_red_black(const Laplace1DMG& lap);
    void gauss_seidel_red_black_local(const Laplace1DMG& lap);
    void init_rhs(Laplace1DMG& lap, const Matrix<dcomplex>& bcmplx);
    bool is_diagonally_dominant(const Laplace1DMG& lap);
    void refine(const Laplace1DMG& lap, Matrix<dcomplex>& fine_error);
    void synchronize_reduced_field(const Laplace1DMG& lap, Matrix<dcomplex>& field);
    void update_solution(const Laplace1DMG& lap);
  };

  void transpose(Matrix<dcomplex> &matrix_transposed, const Matrix<dcomplex> &matrix);

private:
  /// Solver tolerances
  BoutReal rtol, atol;

  /// Maximum number of iterations
  int maxits;

  /// Maximum number of coarse grids
  int max_level;

  /// Maximum number of iterations per grid
  int max_cycle;

  /// Predict when convergence will be reached, and skip expensive convergence
  /// checks at earlier iterations.
  bool predict_exit;

  /// The coefficents in
  /// $D*grad_perp^2(x) + (1/C)*(grad_perp(C))*grad_perp(x) + A*x = b$
  Field2D A, C, D;

  /// Number of unfiltered Fourier modes
  int nmode;

  /// Number of local x, y points
  int ncx, ny;

  /// Information about the grids
  std::vector<Level> levels;

  /// Current y index
  int jy;

  /// Lower-, on- and upper-diagonal terms of the operator matrix
  Tensor<dcomplex> avec, bvec, cvec;


  /// Flag to state whether this is the first time the solver is called
  /// on the point (jy,kz).
  Array<bool> first_call;

  /// Save previous x in Fourier space
  Tensor<dcomplex> x0saved;

  /// Mean number of iterations taken by the solver
  BoutReal ipt_mean_its{0.0};

  /// Counter for the number of times the solver has been called
  int ncalls{0};

  /// Neighbouring processors in the in and out directions
  int proc_in, proc_out;

  /// This processor's unique ID
  int myproc;

  /// Shorthand for localmesh->NXPE
  int nproc;

  /// Array recording whether a kz mode is converged
  Array<bool> converged;

  /// Error interpolated onto the grid one finer than current grid
  Matrix<dcomplex> fine_error;

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
};

#endif // __1DMG_H__
