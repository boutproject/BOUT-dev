// clang-format off
/*!
 * \file naulin_laplace.cxx
 *
 * \brief Iterative solver to handle non-constant-in-z coefficients
 *
 * Scheme suggested by Volker Naulin: solve
 * Delp2(phi[i+1]) + 1/DC(C1*D)*Grad_perp(DC(C2))*Grad_perp(phi[i+1]) + DC(A/D)*phi[i+1]
 *   = rhs(phi[i]) + 1/DC(C1*D)*Grad_perp(DC(C2))*Grad_perp(phi[i]) + DC(A/D)*phi[i]
 * using standard FFT-based solver, iterating to include other terms by
 * evaluating them on rhs using phi from previous iteration.
 * DC part (i.e. Field2D part) of C1*D, C2 and A/D is kept in the FFT inversion
 * to improve convergence by including as much as possible in the direct solve
 * and so that all Neumann boundary conditions can be used at least when
 * DC(A/D)!=0.
 *
 * CHANGELOG
 * =========
 *
 **************************************************************************
 * Copyright 2018 B.D.Dudson, M. Loiten, J. Omotani
 *
 * Contact: Ben Dudson, benjamin.dudson@york.ac.uk
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
 * ## Explanation of the procedure:
 * A way to invert the equation
 * \f$\Omega^D = \nabla\cdot(n\nabla_\perp \phi)\f$
 * invented by Naulin, V.
 * In an orthogonal system, we have that:
 *
 * \f{eqnarray}{
 * \Omega^D &=& \nabla\cdot(n\nabla_\perp \phi)\\
 *       &=& n \nabla_\perp^2 \phi + \nabla n\cdot\nabla_\perp \phi\\
 *       &=& n \Omega + \nabla n\cdot\nabla_\perp \phi\\
 *       &=& n \Omega + \nabla_\perp n\cdot\nabla_\perp \phi
 * \f}
 *
 * Rearranging gives
 *
 * \f{eqnarray}{
 * \Omega  &=& \frac{\Omega^D}{n} - \nabla_\perp \ln(n)\cdot\nabla_\perp \phi\\
 * \nabla_\perp^2 \phi
 * &=& \frac{\Omega^D}{n} - \nabla_\perp \ln(n)\cdot\nabla_\perp \phi
 * \f}
 *
 * In fact we allow for the slightly more general form
 *
 * \f{eqnarray}{
 * \nabla_\perp^2 \phi + <\frac{A}{D}>\phi
 * &=& rhs/D - \frac{1}{D\,C1} \nabla_\perp C2\cdot\nabla_\perp \phi - (\frac{A}{D} - <\frac{A}{D}>)*\phi
 * \f}
 *
 * The iteration can be under-relaxed to help it converge. Amount of under-relaxation is
 * set by the parameter 'underrelax_factor'. 0<underrelax_factor<=1, with
 * underrelax_factor=1 corresponding to no under-relaxation. The amount of
 * under-relaxation is temporarily increased if the iteration starts diverging, the
 * starting value uof underrelax_factor can be set with the initial_underrelax_factor
 * option.
 *
 * The iteration now works as follows:
 *      1. Get the vorticity from
 *         \code{.cpp}
 *         vort = (vortD/n) - grad_perp(ln_n)*grad_perp(phiCur)
 *         [Delp2(phiNext) + 1/DC(C2*D)*grad_perp(DC(C2))*grad_perp(phiNext) + DC(A/D)*phiNext
 *          = b(phiCur)
 *          = (rhs/D) - (1/C1/D*grad_perp(C2)*grad_perp(phiCur) - 1/DC(C2*D)*grad_perp(DC(C2))*grad_perp(phiCur)) - (A/D - DC(A/D))*phiCur]
 *         \endcode
 *         where phiCur is phi of the current iteration
 *         [and DC(f) is the constant-in-z component of f]
 *      2. Invert \f$phi\f$ to find the voricity using
 *         \code{.cpp}
 *         phiNext = invert_laplace_perp(vort)
 *         [set Acoef of laplace_perp solver to DC(A/D)
 *          and C1coef of laplace_perp solver to DC(C1*D)
 *          and C2coef of laplace_perp solver to DC(C2)
 *          then phiNext = invert_laplace_perp(underrelax_factor*b(phiCur) - (1-underrelax_factor)*b(phiPrev))]
 *          where b(phiPrev) is the previous rhs value, which (up to rounding errors) is
 *          the same as the lhs of the direct solver applied to phiCur.
 *         \endcode
 *         where phiNext is the newly obtained \f$phi\f$
 *      3. Calculate the error at phi=phiNext
 *         \code{.cpp}
 *         error3D = Delp2(phiNext) + 1/C1*grad_perp(C2)*grad_perp(phiNext) + A/D*phiNext - rhs/D
 *                 = b(phiCur) - b(phiNext)
 *         as b(phiCur) = Delp2(phiNext) + 1/DC(C2*D)*grad_perp(DC(C2))*grad_perp(phiNext) + DC(A/D)*phiNext
 *         up to rounding errors
 *         \endcode
 *      4. Calculate the infinity norms of the error
 *         \code{.cpp}
 *         EAbsLInf = max(error3D)
 *         ERelLInf = EAbsLInf/sqrt( max((rhs/D)^2) )
 *         \endcode
 *      5. Check whether
 *         \code{.cpp}
 *         EAbsLInf > atol
 *         \endcode
 *          * If yes
 *              * Check whether
 *                \code{.cpp}
 *                ERelLInf > rtol
 *                \endcode
 *              * If yes
 *                  * Check whether
 *                  \code{.cpp}
 *                  EAbsLInf > EAbsLInf(previous step)
 *                  \endcode
 *                    * If yes
 *                      \code{.cpp}
 *                      underrelax_factor *= 0.9
 *                      \endcode
 *                      Restart iteration
 *                    * If no
 *                      * Set
 *                        \code{.cpp}
 *                        phiCur = phiNext
 *                        \endcode
 *                        increase curCount and start from step 1
 *                      * If number of iteration is above maxit, throw exception
 *              * If no
 *                  * Stop: Function returns phiNext
 *          * if no
 *              * Stop: Function returns phiNext
 */
// clang-format on

#include <bout/boutexception.hxx>
#include <bout/coordinates.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/sys/timer.hxx>

#include "naulin_laplace.hxx"

LaplaceNaulin::LaplaceNaulin(Options* opt, const CELL_LOC loc, Mesh* mesh_in,
                             Solver* UNUSED(solver))
    : Laplacian(opt, loc, mesh_in), Acoef(0.0), C1coef(1.0), C2coef(0.0), Dcoef(1.0),
      delp2solver(nullptr), naulinsolver_mean_its(0.), ncalls(0) {

  ASSERT1(opt
          != nullptr); // An Options pointer should always be passed in by LaplaceFactory
  Options& options = *opt;

  Acoef.setLocation(location);
  C1coef.setLocation(location);
  C2coef.setLocation(location);
  Dcoef.setLocation(location);

  // Get options
  rtol = options["rtol"].doc("relative tolerance").withDefault(1.e-7);
  atol = options["atol"].doc("absolute tolerance").withDefault(1.e-20);
  rtol_accept =
      options["rtol_accept"].doc("Accept this rtol after maxits").withDefault(rtol);
  atol_accept =
      options["atol_accept"].doc("Accept this atol after maxits").withDefault(atol);

  maxits = options["maxits"].doc("maximum number of iterations").withDefault(100);
  initial_underrelax_factor =
      options["initial_underrelax_factor"]
          .doc("Initial underrelaxation factor for the fixed point iteration.")
          .withDefault(1.0);
  ASSERT0(initial_underrelax_factor > 0. and initial_underrelax_factor <= 1.);
  underrelax_threshold = options["underrelax_threshold"]
                             .doc("Threshold for relative increase of error in a step "
                                  "that triggers decrease of "
                                  "the underrelaxation factor.")
                             .withDefault(1.5);
  ASSERT0(underrelax_threshold >= 1.);
  underrelax_decrease_factor =
      options["underrelax_decrease_factor"]
          .doc("Factor to decrease underrelax_factor at each stage of the sub-loop if "
               "underrelax_threshold was crossed.")
          .withDefault(0.9);
  ASSERT0(underrelax_decrease_factor < 1. && underrelax_decrease_factor > 0.);
  underrelax_decrease_maxits =
      options["underrelax_decrease_maxits"]
          .doc(
              "Maximum number of iterations in the decreasing-underrelax_factor subcycle "
              "before trying to continue the main iteration loop.")
          .withDefault(10);
  underrelax_recovery =
      options["underrelax_recovery"]
          .doc("Factor to increase underrelax_factor by at the end of a successful "
               "iteration "
               "if it has been decreased below initial_underrelax_factor.")
          .withDefault(1.1);
  ASSERT0(underrelax_recovery >= 1.);
  delp2solver = create(opt->getSection("delp2solver"), location, localmesh);
  std::string delp2type;
  opt->getSection("delp2solver")->get("type", delp2type, "cyclic");
  // Check delp2solver is using an FFT scheme, otherwise it will not exactly
  // invert Delp2 and we will not converge
  ASSERT0(delp2type == "cyclic" || delp2type == "spt" || delp2type == "tri");
  // Use same flags for FFT solver as for NaulinSolver
  delp2solver->setGlobalFlags(getGlobalFlags());
  delp2solver->setInnerBoundaryFlags(getInnerBoundaryFlags());
  delp2solver->setOuterBoundaryFlags(getOuterBoundaryFlags());

  static int naulinsolver_count = 1;
  setPerformanceName(fmt::format("{}{}", "naulinsolver", ++naulinsolver_count));
}

Field3D LaplaceNaulin::solve(const Field3D& rhs, const Field3D& x0) {
  // Rearrange equation so first term is just Delp2(x):
  //   D*Delp2(x) + 1/C1*Grad_perp(C2).Grad_perp(phi) = rhs
  //   -> Delp2(x) + 1/(C1*D)*Grad_perp(C2).Grad_perp(phi) = rhs/D

  Timer timer("invert"); ///< Start timer

  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  ASSERT1(Dcoef.getLocation() == location);
  ASSERT1(C1coef.getLocation() == location);
  ASSERT1(C2coef.getLocation() == location);
  ASSERT1(Acoef.getLocation() == location);
  ASSERT1(localmesh == rhs.getMesh() && localmesh == x0.getMesh());

  Field3D rhsOverD = rhs / Dcoef;

  Field3D C1TimesD = C1coef * Dcoef; // This is needed several times

  // x-component of 1./(C1*D) * Grad_perp(C2)
  Field3D coef_x = DDX(C2coef, location, "C2") / C1TimesD;

  // y-component of 1./(C1*D) * Grad_perp(C2)
  Field3D coef_y = DDY(C2coef, location, "C2") / C1TimesD;

  // z-component of 1./(C1*D) * Grad_perp(C2)
  Field3D coef_z = DDZ(C2coef, location, "FFT") / C1TimesD;

  Field3D AOverD = Acoef / Dcoef;

  // Split coefficients into DC and AC parts so that delp2solver can use DC part.
  // This allows all-Neumann boundary conditions as long as AOverD_DC is non-zero

  Field2D C1coefTimesD_DC = DC(C1TimesD);
  Field2D C2coef_DC = DC(C2coef);

  // Our naming is slightly misleading here, as coef_x_AC may actually have a
  // DC component, as the AC components of C2coef and C1coefTimesD are not
  // necessarily in phase.
  // This is the piece that cannot be passed to an FFT-based Laplacian solver
  // (through our current interface).
  Field3D coef_x_AC = coef_x - DDX(C2coef_DC, location, "C2") / C1coefTimesD_DC;

  // coef_z is a z-derivative so must already have zero DC component

  Field2D AOverD_DC = DC(AOverD);
  Field3D AOverD_AC = AOverD - AOverD_DC;

  delp2solver->setCoefA(AOverD_DC);
  delp2solver->setCoefC1(C1coefTimesD_DC);
  delp2solver->setCoefC2(C2coef_DC);

  // Use this below to normalize error for relative error estimate
  BoutReal RMS_rhsOverD = sqrt(mean(
      SQ(rhsOverD), true,
      "RGN_NOBNDRY")); // use sqrt(mean(SQ)) to make sure we do not divide by zero at a point

  BoutReal error_rel = 1e20, error_abs = 1e20, last_error = error_abs;
  int count = 0;
  int underrelax_count = 0;
  BoutReal underrelax_factor = initial_underrelax_factor;

  auto calc_b_guess = [&](const Field3D& x_in) {
    // Derivatives of x
    Field3D ddx_x = DDX(x_in, location, "C2");
    Field3D ddz_x = DDZ(x_in, location, "FFT");
    return rhsOverD
           - (coords->g11 * coef_x_AC * ddx_x + coords->g33 * coef_z * ddz_x
              + coords->g13 * (coef_x_AC * ddz_x + coef_z * ddx_x))
           - AOverD_AC * x_in;
  };

  auto calc_b_x_pair = [&, this](Field3D b, Field3D x_guess) {
    // Note take a copy of the 'b' argument, because we want to return a copy of it in the
    // result

    if (isInnerBoundaryFlagSet(INVERT_SET) || isOuterBoundaryFlagSet(INVERT_SET)) {
      // This passes in the boundary conditions from x0's guard cells
      copy_x_boundaries(x_guess, x0, localmesh);
    }

    // NB need to pass x_guess in case boundary flags require 'x0', even if delp2solver is
    // not iterative and does not use an initial guess
    Field3D x = delp2solver->solve(b, x_guess);
    localmesh->communicate(x);

    return std::make_pair(b, x);
  };

  Field3D b = calc_b_guess(x0);
  // Need to make a copy of x0 here to make sure we don't change x0
  auto b_x_pair = calc_b_x_pair(b, x0);
  auto b_x_pair_old = b_x_pair;

  while (true) {
    Field3D bnew = calc_b_guess(b_x_pair.second);

    Field3D error3D = b_x_pair.first - bnew;
    error_abs = max(abs(error3D, "RGN_NOBNDRY"), true, "RGN_NOBNDRY");
    error_rel = error_abs / RMS_rhsOverD;

    if (error_rel < rtol or error_abs < atol) {
      break;
    }

    ++count;
    if (count > maxits) {
      // Perhaps accept a worse solution
      if (error_rel < rtol_accept or error_abs < atol_accept) {
        break;
      }
      throw BoutException(
          "LaplaceNaulin error: Not converged within maxits={:d} iterations.", maxits);
    }

    int local_count = 0;
    while (error_abs > last_error * underrelax_threshold) {
      // Iteration seems to be diverging... try underrelaxing and restart
      underrelax_factor *= underrelax_decrease_factor;
      ++underrelax_count;

      // Restart from b_x_pair_old - that was our best guess
      bnew = calc_b_guess(b_x_pair_old.second);
      b_x_pair = calc_b_x_pair(underrelax_factor * bnew
                                   + (1. - underrelax_factor) * b_x_pair_old.first,
                               b_x_pair_old.second);

      bnew = calc_b_guess(b_x_pair.second);

      error3D = b_x_pair.first - bnew;
      error_abs = max(abs(error3D, "RGN_NOBNDRY"), true, "RGN_NOBNDRY");
      error_rel = error_abs / RMS_rhsOverD;

      ++local_count;
      if (local_count > underrelax_decrease_maxits) {
        // Give up on trying to underrelax. Attempt to continue iteration anyway...
        break;
      }

      // effectively another iteration, so increment the counter
      ++count;
      if (count > maxits) {
        if (error_rel < rtol_accept or error_abs < atol_accept) {
          break;
        }
        throw BoutException(
            "LaplaceNaulin error: Not converged within maxits={:d} iterations.", maxits);
      }
    }

    // Might have met convergence criterion while in underrelaxation loop
    if (error_rel < rtol or error_abs < atol) {
      break;
    }

    if (underrelax_factor < initial_underrelax_factor) {
      underrelax_factor *= underrelax_recovery;
      if (underrelax_factor > initial_underrelax_factor) {
        underrelax_factor = initial_underrelax_factor;
      }
    }

    last_error = error_abs;
    b_x_pair_old = b_x_pair;
    b_x_pair = calc_b_x_pair(underrelax_factor * bnew
                                 + (1. - underrelax_factor) * b_x_pair.first,
                             b_x_pair.second);
  }

  ++ncalls;
  naulinsolver_mean_its =
      (naulinsolver_mean_its * BoutReal(ncalls - 1) + BoutReal(count)) / BoutReal(ncalls);
  naulinsolver_mean_underrelax_counts =
      (naulinsolver_mean_underrelax_counts * BoutReal(ncalls - 1)
       + BoutReal(underrelax_count))
      / BoutReal(ncalls);

  checkData(b_x_pair.second);

  return b_x_pair.second;
}

void LaplaceNaulin::copy_x_boundaries(Field3D& x, const Field3D& x0, Mesh* localmesh) {
  if (localmesh->firstX()) {
    for (int i = localmesh->xstart - 1; i >= 0; --i) {
      for (int j = localmesh->ystart; j <= localmesh->yend; ++j) {
        for (int k = 0; k < localmesh->LocalNz; ++k) {
          x(i, j, k) = x0(i, j, k);
        }
      }
    }
  }
  if (localmesh->lastX()) {
    for (int i = localmesh->xend + 1; i < localmesh->LocalNx; ++i) {
      for (int j = localmesh->ystart; j <= localmesh->yend; ++j) {
        for (int k = 0; k < localmesh->LocalNz; ++k) {
          x(i, j, k) = x0(i, j, k);
        }
      }
    }
  }
}

void LaplaceNaulin::outputVars(Options& output_options,
                               const std::string& time_dimension) const {
  output_options[fmt::format("{}_mean_its", getPerformanceName())].assignRepeat(
      naulinsolver_mean_its, time_dimension);
  output_options[fmt::format("{}_mean_underrelax_counts", getPerformanceName())]
      .assignRepeat(naulinsolver_mean_underrelax_counts, time_dimension);
}
