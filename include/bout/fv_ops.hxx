/*
  Finite-volume discretisation methods. Flux-conservative form
 */

#ifndef BOUT_FV_OPS_H
#define BOUT_FV_OPS_H

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/build_defines.hxx"
#include "bout/coordinates.hxx"
#include "bout/field.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"
#include "bout/vector2d.hxx"

#include <cmath>
#include <cstdlib>

namespace FV {
/*!
 * Div ( a Grad_perp(f) ) -- ∇⊥ ( a ⋅ ∇⊥ f) -- Vorticity
 */
Field3D Div_a_Grad_perp(const Field3D& a_coef, const Field3D& f);

[[deprecated("Please use Div_a_Grad_perp instead")]] inline Field3D
Div_a_Laplace_perp(const Field3D& a_coef, const Field3D& f) {
  return Div_a_Grad_perp(a_coef, f);
}

/*!
   * Divergence of a parallel diffusion Div( k * Grad_par(f) )
   */
Field3D Div_par_K_Grad_par(const Field3D& k, const Field3D& f, bool bndry_flux = true);

/*!
   * 4th-order derivative in Y, using derivatives
   * on cell boundaries.
   *
   * A one-sided 3rd-order derivative, given a value
   * at a boundary is:
   *
   * d3f/dx3 ~= 16/5 f_b - 6 f_0 + 4 f_1 - 6/5 f_2
   *
   * where f_b is the value on the boundary; f_0 is the cell
   * to the left of the boundary; f_1 to the left of f_0 and f_2
   * to the left of f_1
   *
   *    f_2 | f_1 | f_0 |
   *                   f_b
   *
   * NB: Uses to/from FieldAligned coordinates
   *
   * No fluxes through domain boundaries
   */
Field3D D4DY4(const Field3D& d_in, const Field3D& f);

/*!
   * 4th-order dissipation term
   *
   *
   * A one-sided 3rd-order derivative, given a value
   * at a boundary is:
   *
   * d3f/dx3 ~= 16/5 f_b - 6 f_0 + 4 f_1 - 6/5 f_2
   *
   * where f_b is the value on the boundary; f_0 is the cell
   * to the left of the boundary; f_1 to the left of f_0 and f_2
   * to the left of f_1
   *
   *    f_2 | f_1 | f_0 |
   *                   f_b
   */
Field3D D4DY4_index(const Field3D& f_in, bool bndry_flux = true);

/*!
   * Stencil used for Finite Volume calculations
   * which includes cell face values L and R
   */
struct Stencil1D {
  // Cell centre values
  BoutReal c = BoutNaN, m = BoutNaN, p = BoutNaN, mm = BoutNaN, pp = BoutNaN;

  // Left and right cell face values
  BoutReal L = BoutNaN, R = BoutNaN;
};

/*!
   * First order upwind for testing
   */
struct Upwind {
  void operator()(Stencil1D& n) { n.L = n.R = n.c; }
};

/*!
   * Fromm method
   */
struct Fromm {
  void operator()(Stencil1D& n) {
    n.L = n.c - 0.25 * (n.p - n.m);
    n.R = n.c + 0.25 * (n.p - n.m);
  }
};

/*!
   * Second order slope limiter method
   * 
   * Limits slope to minimum absolute value
   * of left and right gradients. If at a maximum
   * or minimum slope set to zero, i.e. reverts
   * to first order upwinding
   */
struct MinMod {
  void operator()(Stencil1D& n) {
    // Choose the gradient within the cell
    // as the minimum (smoothest) solution
    const BoutReal slope = _minmod(n.p - n.c, n.c - n.m);
    n.L = n.c - 0.5 * slope;
    n.R = n.c + 0.5 * slope;
  }

private:
  /*!
     * Internal helper function for minmod slope limiter
     *
     * If the inputs have different signs then
     * returns zero, otherwise chooses the value
     * with the minimum magnitude.
     */
  static BoutReal _minmod(BoutReal lhs, BoutReal rhs) {
    if (lhs * rhs <= 0.0) {
      return 0.0;
    }

    if (std::abs(lhs) < std::abs(rhs)) {
      return lhs;
    }
    return rhs;
  }
};

/*!
   * Monotonised Central (MC) second order slope limiter (Van Leer)
   * 
   * Limits the slope based on taking the slope with 
   * the minimum absolute value from central, 2*left and
   * 2*right. If any of these slopes have different signs
   * then the slope reverts to zero (i.e. 1st-order upwinding).
   */
struct MC {
  void operator()(Stencil1D& n) {
    const BoutReal slope = minmod(2. * (n.p - n.c),  // 2*right difference
                                  0.5 * (n.p - n.m), // Central difference
                                  2. * (n.c - n.m)); // 2*left difference
    n.L = n.c - 0.5 * slope;
    n.R = n.c + 0.5 * slope;
  }

private:
  // Return zero if any signs are different
  // otherwise return the value with the minimum magnitude
  static BoutReal minmod(BoutReal lhs, BoutReal centre, BoutReal rhs) {
    // if any of the signs are different, return zero gradient
    if ((lhs * centre <= 0.0) || (lhs * rhs <= 0.0)) {
      return 0.0;
    }

    // Return the minimum absolute value
    return SIGN(lhs) * BOUTMIN(std::abs(lhs), std::abs(centre), std::abs(rhs));
  }
};

/*!
   * Communicate fluxes between processors
   * Takes values in guard cells, and adds them to cells
   */
void communicateFluxes(Field3D& f);

/// Finite volume parallel divergence
///
/// Preserves the sum of f*J*dx*dy*dz over the domain
///
/// @param[in] f_in   The field being advected.
///                   This will be reconstructed at cell faces
///                   using the given CellEdges method
/// @param[in] v_in   The advection velocity.
///                   This will be interpolated to cell boundaries
///                   using linear interpolation
/// @param[in] wave_speed_in  Local maximum speed of all waves in the system at each
//                            point in space
/// @param[in] fixflux     Fix the flux at the boundary to be the value at the
///                        midpoint (for boundary conditions)
///
/// NB: Uses to/from FieldAligned coordinates
template <typename CellEdges = MC>
Field3D Div_par(const Field3D& f_in, const Field3D& v_in, const Field3D& wave_speed_in,
                bool fixflux = true) {

  ASSERT1_FIELDS_COMPATIBLE(f_in, v_in);
  ASSERT1_FIELDS_COMPATIBLE(f_in, wave_speed_in);

  Mesh* mesh = f_in.getMesh();

  CellEdges cellboundary;

  ASSERT2(f_in.getDirectionY() == v_in.getDirectionY());
  ASSERT2(f_in.getDirectionY() == wave_speed_in.getDirectionY());
  const bool are_unaligned =
      ((f_in.getDirectionY() == YDirectionType::Standard)
       and (v_in.getDirectionY() == YDirectionType::Standard)
       and (wave_speed_in.getDirectionY() == YDirectionType::Standard));

  const Field3D f = are_unaligned ? toFieldAligned(f_in, "RGN_NOX") : f_in;
  const Field3D v = are_unaligned ? toFieldAligned(v_in, "RGN_NOX") : v_in;
  Field3D wave_speed =
      are_unaligned ? toFieldAligned(wave_speed_in, "RGN_NOX") : wave_speed_in;

  Coordinates* coord = f_in.getCoordinates();

  Field3D result{zeroFrom(f)};

  // Only need one guard cell, so no need to communicate fluxes
  // Instead calculate in guard cells to preserve fluxes
  int y_start = mesh->ystart - 1;
  int y_end = mesh->yend + 1;

  for (int i = mesh->xstart; i <= mesh->xend; i++) {

    if (!mesh->firstY(i) || mesh->periodicY(i)) {
      // Calculate in guard cell to get fluxes consistent between processors
      y_start = mesh->ystart - 1;
    } else {
      // Don't include the boundary cell. Note that this implies special
      // handling of boundaries later
      y_start = mesh->ystart;
    }

    if (!mesh->lastY(i) || mesh->periodicY(i)) {
      // Calculate in guard cells
      y_end = mesh->yend + 1;
    } else {
      // Not in boundary cells
      y_end = mesh->yend;
    }

    using std::sqrt;

    for (int j = y_start; j <= y_end; j++) {
      // Pre-calculate factors which multiply fluxes
#if not(BOUT_USE_METRIC_3D)
      // For right cell boundaries
      const BoutReal common_factor_right =
          (coord->J(i, j) + coord->J(i, j + 1))
          / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j + 1)));

      const BoutReal flux_factor_rc =
          common_factor_right / (coord->dy(i, j) * coord->J(i, j));
      const BoutReal flux_factor_rp =
          common_factor_right / (coord->dy(i, j + 1) * coord->J(i, j + 1));

      // For left cell boundaries
      const BoutReal common_factor_left =
          (coord->J(i, j) + coord->J(i, j - 1))
          / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j - 1)));

      const BoutReal flux_factor_lc =
          common_factor_left / (coord->dy(i, j) * coord->J(i, j));
      const BoutReal flux_factor_lm =
          common_factor_left / (coord->dy(i, j - 1) * coord->J(i, j - 1));
#endif
      for (int k = 0; k < mesh->LocalNz; k++) {
#if BOUT_USE_METRIC_3D
        // For right cell boundaries
        const BoutReal common_factor_right =
            (coord->J(i, j, k) + coord->J(i, j + 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j + 1, k)));

        const BoutReal flux_factor_rc =
            common_factor_right / (coord->dy(i, j, k) * coord->J(i, j, k));
        const BoutReal flux_factor_rp =
            common_factor_right / (coord->dy(i, j + 1, k) * coord->J(i, j + 1, k));

        // For left cell boundaries
        const BoutReal common_factor_left =
            (coord->J(i, j, k) + coord->J(i, j - 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j - 1, k)));

        const BoutReal flux_factor_lc =
            common_factor_left / (coord->dy(i, j, k) * coord->J(i, j, k));
        const BoutReal flux_factor_lm =
            common_factor_left / (coord->dy(i, j - 1, k) * coord->J(i, j - 1, k));
#endif

        ////////////////////////////////////////////
        // Reconstruct f at the cell faces
        // This calculates s.R and s.L for the Right and Left
        // face values on this cell

        // Reconstruct f at the cell faces
        Stencil1D stencil;
        stencil.c = f(i, j, k);
        stencil.m = f(i, j - 1, k);
        stencil.p = f(i, j + 1, k);

        cellboundary(stencil); // Calculate s.R and s.L

        ////////////////////////////////////////////
        // Right boundary

        // Calculate velocity at right boundary (y+1/2)
        const BoutReal vpar_right = 0.5 * (v(i, j, k) + v(i, j + 1, k));
        BoutReal flux = BoutNaN;

        if (mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
          // Last point in domain

          const BoutReal bndryval = 0.5 * (stencil.c + stencil.p);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar_right;
          } else {
            // Add flux due to difference in boundary values
            flux = stencil.R * vpar_right + wave_speed(i, j, k) * (stencil.R - bndryval);
          }
        } else {

          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k));

          if (vpar_right > amax) {
            // Supersonic flow out of this cell
            flux = stencil.R * vpar_right;
          } else if (vpar_right < -amax) {
            // Supersonic flow into this cell
            flux = 0.0;
          } else {
            // Subsonic flow, so a mix of right and left fluxes
            flux = stencil.R * 0.5 * (vpar_right + amax);
          }
        }

        result(i, j, k) += flux * flux_factor_rc;
        result(i, j + 1, k) -= flux * flux_factor_rp;

        ////////////////////////////////////////////
        // Calculate at left boundary

        const BoutReal vpar_left = 0.5 * (v(i, j, k) + v(i, j - 1, k));

        if (mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
          // First point in domain
          const BoutReal bndryval = 0.5 * (stencil.c + stencil.m);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar_left;
          } else {
            // Add flux due to difference in boundary values
            flux = stencil.L * vpar_left - wave_speed(i, j, k) * (stencil.L - bndryval);
          }
        } else {

          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k));

          if (vpar_left < -amax) {
            // Supersonic out of this cell
            flux = stencil.L * vpar_left;
          } else if (vpar_left > amax) {
            // Supersonic into this cell
            flux = 0.0;
          } else {
            flux = stencil.L * 0.5 * (vpar_left - amax);
          }
        }

        result(i, j, k) -= flux * flux_factor_lc;
        result(i, j - 1, k) += flux * flux_factor_lm;
      }
    }
  }
  return are_unaligned ? fromFieldAligned(result, "RGN_NOBNDRY") : result;
}

/*!
   * Div ( n * v )  -- Magnetic drifts
   *
   * This uses the expression
   * 
   * Div( A ) = 1/J * d/di ( J * A^i )
   * 
   * Hence the input vector should be contravariant
   *
   * Note: Uses to/from FieldAligned
   *
   */
template <typename CellEdges = MC>
Field3D Div_f_v(const Field3D& n_in, const Vector3D& v, bool bndry_flux) {
  ASSERT1(n_in.getLocation() == v.getLocation());
  ASSERT1_FIELDS_COMPATIBLE(n_in, v.x);

  const Mesh* mesh = n_in.getMesh();

  CellEdges cellboundary;

  const Coordinates* coord = n_in.getCoordinates();

  if (v.covariant) {
    // Got a covariant vector instead
    throw BoutException("Div_f_v passed a covariant v");
  }

  Field3D result{zeroFrom(n_in)};
  const auto& J = coord->J;

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    // Calculate velocities
    const BoutReal v_U = 0.25 * (v.z[i.zp()] + v.z[i]) * (J[i.zp()] + J[i]);
    const BoutReal v_D = 0.25 * (v.z[i.zm()] + v.z[i]) * (J[i.zm()] + J[i]);
    const BoutReal v_L = 0.25 * (v.x[i.xm()] + v.x[i]) * (J[i.xm()] + J[i]);
    const BoutReal v_R = 0.25 * (v.x[i.xp()] + v.x[i]) * (J[i.xp()] + J[i]);

    // X direction
    Stencil1D stencil;
    stencil.c = n_in[i];
    stencil.m = n_in[i.xm()];
    stencil.mm = n_in[i.xmm()];
    stencil.p = n_in[i.xp()];
    stencil.pp = n_in[i.xpp()];

    cellboundary(stencil);

    if ((i.x() == mesh->xend) && (mesh->lastX())) {
      // At right boundary in X
      if (bndry_flux) {
        BoutReal flux = BoutNaN;
        if (v_R > 0.0) {
          // Flux to boundary
          flux = v_R * stencil.R;
        } else {
          // Flux in from boundary
          flux = v_R * 0.5 * (n_in[i.xp()] + n_in[i]);
        }
        result[i] += flux / (coord->dx[i] * coord->J[i]);
        result[i.xp()] -= flux / (coord->dx[i.xp()] * coord->J[i.xp()]);
      }
    } else {
      // Not at a boundary
      if (v_R > 0.0) {
        // Flux out into next cell
        const BoutReal flux = v_R * stencil.R;
        result[i] += flux / (coord->dx[i] * coord->J[i]);
        result[i.xp()] -= flux / (coord->dx[i.xp()] * coord->J[i.xp()]);
      }
    }

    // Left side

    if ((i.x() == mesh->xstart) && (mesh->firstX())) {
      // At left boundary in X

      if (bndry_flux) {
        BoutReal flux = BoutNaN;
        if (v_L < 0.0) {
          // Flux to boundary
          flux = v_L * stencil.L;
        } else {
          // Flux in from boundary
          flux = v_L * 0.5 * (n_in[i.xm()] + n_in[i]);
        }
        result[i] -= flux / (coord->dx[i] * coord->J[i]);
        result[i.xm()] += flux / (coord->dx[i.xm()] * coord->J[i.xm()]);
      }
    } else {
      // Not at a boundary
      if (v_L < 0.0) {
        const BoutReal flux = v_L * stencil.L;
        result[i] -= flux / (coord->dx[i] * coord->J[i]);
        result[i.xm()] += flux / (coord->dx[i.xm()] * coord->J[i.xm()]);
      }
    }

    /// NOTE: Need to communicate fluxes

    // Z direction
    stencil.m = n_in[i.zm()];
    stencil.mm = n_in[i.zmm()];
    stencil.p = n_in[i.zp()];
    stencil.pp = n_in[i.zpp()];

    cellboundary(stencil);

    if (v_U > 0.0) {
      const BoutReal flux = v_U * stencil.R;
      result[i] += flux / (coord->J[i] * coord->dz[i]);
      result[i.zp()] -= flux / (coord->J[i.zp()] * coord->dz[i.zp()]);
    }
    if (v_D < 0.0) {
      const BoutReal flux = v_D * stencil.L;
      result[i] -= flux / (coord->J[i] * coord->dz[i]);
      result[i.zm()] += flux / (coord->J[i.zm()] * coord->dz[i.zm()]);
    }
  }

  communicateFluxes(result);

  // Y advection
  // Currently just using simple centered differences
  // so no fluxes need to be exchanged

  const Field3D n_aligned = toFieldAligned(n_in, "RGN_NOX");
  const Field3D v_y = toFieldAligned(v.y, "RGN_NOX");

  Field3D yresult = 0.0;
  yresult.setDirectionY(YDirectionType::Aligned);

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    // Y velocities on y boundaries
    const BoutReal v_U = 0.25 * (v_y[i] + v_y[i.yp()]) * (coord->J[i] + coord->J[i.yp()]);
    const BoutReal v_D = 0.25 * (v_y[i] + v_y[i.ym()]) * (coord->J[i] + coord->J[i.ym()]);

    // n (advected quantity) on y boundaries
    // Note: Use unshifted n_in variable
    const BoutReal n_U = 0.5 * (n_aligned[i] + n_aligned[i.yp()]);
    const BoutReal n_D = 0.5 * (n_aligned[i] + n_aligned[i.ym()]);

    yresult[i] = (n_U * v_U - n_D * v_D) / (coord->J[i] * coord->dy[i]);
  }
  return result + fromFieldAligned(yresult, "RGN_NOBNDRY");
}

/*!
   * X-Z Finite Volume diffusion operator
   */
Field3D Div_Perp_Lap(const Field3D& a_coeff, const Field3D& f,
                     CELL_LOC outloc = CELL_DEFAULT);
} // namespace FV
#endif // BOUT_FV_OPS_H
