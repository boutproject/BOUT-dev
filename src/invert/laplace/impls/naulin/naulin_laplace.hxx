/// \file
/// Iterative solver to handle non-constant-in-z coefficients
///
/// Scheme suggested by Volker Naulin: solve
/// \f{eqnarray}
///   \nabla^2(\phi[i+1])
///      + 1/DC(C_1 D)\nabla_\perp(DC(C_2))\nabla_\perp(\phi[i+1])
///      + DC(A/D)\phi[i+1] \\
///   = rhs(\phi[i])
///        + 1/DC(C_1 D)\nabla_\perp(DC(C_2))\nabla_\perp(\phi[i])
///        + DC(A/D)\phi[i]
/// \f}
///
/// using standard FFT-based solver, iterating to include other terms
/// by evaluating them on ``rhs`` using \f$\phi\f$ from previous
/// iteration.  DC part (i.e. `Field2D` part) of \f$C_1 D\f$,
/// \f$C_2\f$ and \f$A,D\f$ is kept in the FFT inversion to improve
/// convergence by including as much as possible in the direct solve
/// and so that all Neumann boundary conditions can be used at least
/// when \f$DC(A/D)!=0\f$
///
/// Explanation of the procedure
/// ----------------------------
///
/// A way to invert the equation
/// \f$\Omega^D = \nabla\cdot(n\nabla_\perp \phi)\f$
/// invented by Naulin, V.
/// In an orthogonal system, we have that:
///
/// \f{eqnarray}{
/// \Omega^D &=& \nabla\cdot(n\nabla_\perp \phi)\                    \
///       &=& n \nabla_\perp^2 \phi + \nabla n\cdot\nabla_\perp \phi\\
///       &=& n \Omega + \nabla n\cdot\nabla_\perp \phi\\
///       &=& n \Omega + \nabla_\perp n\cdot\nabla_\perp \phi
/// \f}
///
/// Rearranging gives
///
/// \f{eqnarray}{
/// \Omega  &=& \frac{\Omega^D}{n} - \nabla_\perp \ln(n)\cdot\nabla_\perp \phi\ \
/// \nabla_\perp^2 \phi
/// &=& \frac{\Omega^D}{n} - \nabla_\perp \ln(n)\cdot\nabla_\perp \phi
/// \f}
///
/// In fact we allow for the slightly more general form
///
/// \f{eqnarray}{
/// \nabla_\perp^2 \phi + <\frac{A}{D}>\phi
/// &=& rhs/D - \frac{1}{D\,C1} \nabla_\perp C2\cdot\nabla_\perp \phi - (\frac{A}{D} - <\frac{A}{D}>)\phi
/// \f}
///
/// The iteration can be under-relaxed to help it converge. Amount of under-relaxation is
/// set by the parameter 'underrelax_factor'. 0<underrelax_factor<=1, with
/// underrelax_factor=1 corresponding to no under-relaxation. The amount of
/// under-relaxation is temporarily increased if the iteration starts diverging, the
/// starting value uof underrelax_factor can be set with the initial_underrelax_factor
/// option.
///
/// The iteration now works as follows:
///  1. Get the vorticity from
///     \code{.cpp}
///     vort = (vortD/n) - grad_perp(ln_n)*grad_perp(phiCur)
///     [Delp2(phiNext) + 1/DC(C2*D)*grad_perp(DC(C2))*grad_perp(phiNext) + DC(A/D)*phiNext
///      = b(phiCur)
///      = (rhs/D) - (1/C1/D*grad_perp(C2)*grad_perp(phiCur) - 1/DC(C2*D)*grad_perp(DC(C2))*grad_perp(phiCur)) - (A/D - DC(A/D))*phiCur]
///    \endcode
///    where phiCur is phi of the current iteration
///    [and DC(f) is the constant-in-z component of f]
/// 2. Invert \f$phi\f$ to find the voricity using
///    \code{.cpp}
///    phiNext = invert_laplace_perp(vort)
///    [set Acoef of laplace_perp solver to DC(A/D)
///     and C1coef of laplace_perp solver to DC(C1*D)
///     and C2coef of laplace_perp solver to DC(C2)
///     then phiNext = invert_laplace_perp(underrelax_factor*b(phiCur) - (1-underrelax_factor)*b(phiPrev))]
///     where b(phiPrev) is the previous rhs value, which (up to rounding errors) is
///     the same as the lhs of the direct solver applied to phiCur.
///    \endcode
///    where phiNext is the newly obtained \f$phi\f$
/// 3. Calculate the error at phi=phiNext
///    \code{.cpp}
///    error3D = Delp2(phiNext) + 1/C1*grad_perp(C2)*grad_perp(phiNext) + A/D*phiNext - rhs/D
///            = b(phiCur) - b(phiNext)
///    as b(phiCur) = Delp2(phiNext) + 1/DC(C2*D)*grad_perp(DC(C2))*grad_perp(phiNext) + DC(A/D)*phiNext
///    up to rounding errors
///    \endcode
/// 4. Calculate the infinity norms of the error
///    \code{.cpp}
///    EAbsLInf = max(error3D)
///    ERelLInf = EAbsLInf/sqrt( max((rhs/D)^2) )
///    \endcode
/// 5. Check whether
///    \code{.cpp}
///    EAbsLInf > atol
///    \endcode
///     * If yes
///         * Check whether
///           \code{.cpp}
///           ERelLInf > rtol
///           \endcode
///         * If yes
///             * Check whether
///             \code{.cpp}
///             EAbsLInf > EAbsLInf(previous step)
///             \endcode
///               * If yes
///                 \code{.cpp}
///                 underrelax_factor *= 0.9
///                 \endcode
///                 Restart iteration
///               * If no
///                 * Set
///                   \code{.cpp}
///                   phiCur = phiNext
///                   \endcode
///                   increase curCount and start from step 1
///                 * If number of iteration is above maxit, throw exception
///         * If no
///             * Stop: Function returns phiNext
///     * if no
///         * Stop: Function returns phiNext

/**************************************************************************
 * Copyright 2018 B.D.Dudson, M. Loiten, J. Omotani
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

class LaplaceNaulin;

#ifndef BOUT_LAP_NAULIN_H
#define BOUT_LAP_NAULIN_H

#include <bout/invert_laplace.hxx>
#include <bout/options.hxx>

namespace {
RegisterLaplace<LaplaceNaulin> registerlaplacenaulin(LAPLACE_NAULIN);
}

/// Solves the 2D Laplacian equation
/*!
 * 
 */
class LaplaceNaulin : public Laplacian {
public:
  LaplaceNaulin(Options* opt = NULL, const CELL_LOC loc = CELL_CENTRE,
                Mesh* mesh_in = nullptr, Solver* solver = nullptr);
  ~LaplaceNaulin() = default;

  using Laplacian::setCoefA;
  using Laplacian::setCoefC;
  using Laplacian::setCoefC1;
  using Laplacian::setCoefC2;
  using Laplacian::setCoefD;
  using Laplacian::setCoefEx;
  using Laplacian::setCoefEz;

  // ACoef is not implemented because the delp2solver that we use can probably
  // only handle a Field2D for Acoef, but we would have to pass Acoef/Dcoef,
  // where we allow Dcoef to be a Field3D
  void setCoefA(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Acoef = val;
  }
  void setCoefA(const Field3D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Acoef = val;
  }
  void setCoefC(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    setCoefC1(val);
    setCoefC2(val);
  }
  void setCoefC(const Field3D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    setCoefC1(val);
    setCoefC2(val);
  }
  void setCoefC1(const Field3D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1coef = val;
  }
  void setCoefC1(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1coef = val;
  }
  void setCoefC2(const Field3D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C2coef = val;
  }
  void setCoefC2(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C2coef = val;
  }
  void setCoefD(const Field3D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Dcoef = val;
  }
  void setCoefD(const Field2D& val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Dcoef = val;
  }
  void setCoefEx(const Field2D& UNUSED(val)) override {
    throw BoutException("LaplaceNaulin does not have Ex coefficient");
  }
  void setCoefEz(const Field2D& UNUSED(val)) override {
    throw BoutException("LaplaceNaulin does not have Ez coefficient");
  }

  bool uses3DCoefs() const override { return true; }

  using Laplacian::solve;

  FieldPerp solve(const FieldPerp& b) override { return solve(b, b); }
  FieldPerp solve(const FieldPerp& UNUSED(b), const FieldPerp& UNUSED(x0)) override {
    throw BoutException(
        "LaplaceNaulin has no solve(FieldPerp), must call solve(Field3D)");
  }
  Field3D solve(const Field3D& b, const Field3D& x0) override;
  Field3D solve(const Field3D& b) override { return solve(b, zeroFrom(b)); }

  // Override flag-setting methods to set delp2solver's flags as well
  void setGlobalFlags(int f) override {
    Laplacian::setGlobalFlags(f);
    delp2solver->setGlobalFlags(f);
  }
  void setInnerBoundaryFlags(int f) override {
    Laplacian::setInnerBoundaryFlags(f);
    delp2solver->setInnerBoundaryFlags(f);
  }
  void setOuterBoundaryFlags(int f) override {
    Laplacian::setOuterBoundaryFlags(f);
    delp2solver->setOuterBoundaryFlags(f);
  }

  BoutReal getMeanIterations() const { return naulinsolver_mean_its; }
  void resetMeanIterations() { naulinsolver_mean_its = 0; }

  void outputVars(Options& output_options,
                  const std::string& time_dimension) const override;

private:
  LaplaceNaulin(const LaplaceNaulin&);
  LaplaceNaulin& operator=(const LaplaceNaulin&);
  Field3D Acoef, C1coef, C2coef, Dcoef;

  /// Laplacian solver used to solve the equation with constant-in-z coefficients
  std::unique_ptr<Laplacian> delp2solver{nullptr};

  /// Solver tolerances
  BoutReal rtol, atol;

  /// Maximum number of iterations
  int maxits;

  /// Initial choice for under-relaxation factor, should be greater than 0 and
  /// less than or equal to 1. Value of 1 means no underrelaxation
  BoutReal initial_underrelax_factor{1.};

  /// Mean number of iterations taken by the solver
  BoutReal naulinsolver_mean_its;

  /// Mean number of times the underrelaxation factor is reduced
  BoutReal naulinsolver_mean_underrelax_counts{0.};

  /// Counter for the number of times the solver has been called
  int ncalls;

  /// Copy the boundary guard cells from the input 'initial guess' x0 into x.
  /// These may be used to set non-zero-value boundary conditions
  void copy_x_boundaries(Field3D& x, const Field3D& x0, Mesh* mesh);
};

#endif // BOUT_LAP_NAULIN_H
