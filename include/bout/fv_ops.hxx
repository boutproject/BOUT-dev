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
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/output_bout_types.hxx" // NOLINT(unused-includes, misc-include-cleaner)
#include "bout/region.hxx"
#include "bout/utils.hxx"
#include "bout/vector2d.hxx"

#include <cmath>

namespace FV {
/*!
 * Div ( a Grad_perp(f) ) -- ∇⊥ ( a ⋅ ∇⊥ f) -- Vorticity
 */
Field3D Div_a_Grad_perp(const Field3D& a, const Field3D& f);

[[deprecated("Please use Div_a_Grad_perp instead")]] inline Field3D
Div_a_Laplace_perp(const Field3D& a, const Field3D& f) {
  return Div_a_Grad_perp(a, f);
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
Field3D D4DY4(const Field3D& d, const Field3D& f);

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
Field3D D4DY4_Index(const Field3D& f, bool bndry_flux = true);

/*!
   * Stencil used for Finite Volume calculations
   * which includes cell face values L and R
   */
struct Stencil1D {
  /// Cell centre values
  BoutReal c;
  BoutReal m;
  BoutReal p;
  BoutReal mm = BoutNaN;
  BoutReal pp = BoutNaN;

  /// Left cell face value
  BoutReal L = BoutNaN;
  /// Right cell face value
  BoutReal R = BoutNaN;
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
    n.L = n.c - (0.25 * (n.p - n.m));
    n.R = n.c + (0.25 * (n.p - n.m));
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
    n.L = n.c - (0.5 * slope);
    n.R = n.c + (0.5 * slope);
  }

private:
  /*!
     * Internal helper function for minmod slope limiter
     *
     * If the inputs have different signs then
     * returns zero, otherwise chooses the value
     * with the minimum magnitude.
     */
  static BoutReal _minmod(BoutReal a, BoutReal b) {
    if (a * b <= 0.0) {
      return 0.0;
    }

    if (fabs(a) < fabs(b)) {
      return a;
    }
    return b;
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
    n.L = n.c - (0.5 * slope);
    n.R = n.c + (0.5 * slope);
  }

private:
  // Return zero if any signs are different
  // otherwise return the value with the minimum magnitude
  static BoutReal minmod(BoutReal a, BoutReal b, BoutReal c) {
    // if any of the signs are different, return zero gradient
    if ((a * b <= 0.0) || (a * c <= 0.0)) {
      return 0.0;
    }

    // Return the minimum absolute value
    return SIGN(a) * BOUTMIN(fabs(a), fabs(b), fabs(c));
  }
};

/// Superbee limiter
///
/// This corresponds to the limiter function
///    φ(r) = max(0, min(2r, 1), min(r,2)
///
/// The value at cell right (i.e. i + 1/2) is:
///
///   n.R = n.c - φ(r) (n.c - (n.p + n.c)/2)
///       = n.c + φ(r) (n.p - n.c)/2
///
/// Four regimes:
///  a) r < 1/2 -> φ(r) = 2r
///     n.R = n.c + gL
///  b) 1/2 < r < 1 -> φ(r) = 1
///     n.R = n.c + gR/2
///  c) 1 < r < 2 -> φ(r) = r
///     n.R = n.c + gL/2
///  d) 2 < r  -> φ(r) = 2
///     n.R = n.c + gR
///
///  where the left and right gradients are:
///   gL = n.c - n.m
///   gR = n.p - n.c
///
struct Superbee {
  void operator()(Stencil1D& n) {
    const BoutReal gL = n.c - n.m;
    const BoutReal gR = n.p - n.c;

    // r = gL / gR
    // Limiter is φ(r)
    if (gL * gR < 0) {
      // Different signs => Zero gradient
      n.L = n.R = n.c;
    } else {
      const BoutReal sign = SIGN(gL);
      const BoutReal abs_gL = fabs(gL);
      const BoutReal abs_gR = fabs(gR);
      const BoutReal half_slope =
          sign * BOUTMAX(BOUTMIN(abs_gL, 0.5 * abs_gR), BOUTMIN(abs_gR, 0.5 * abs_gL));
      n.L = n.c - half_slope;
      n.R = n.c + half_slope;
    }
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

  Mesh const* mesh = f_in.getMesh();

  CellEdges cellboundary;

  ASSERT2(f_in.getDirectionY() == v_in.getDirectionY());
  ASSERT2(f_in.getDirectionY() == wave_speed_in.getDirectionY());
  const bool are_unaligned =
      ((f_in.getDirectionY() == YDirectionType::Standard)
       and (v_in.getDirectionY() == YDirectionType::Standard)
       and (wave_speed_in.getDirectionY() == YDirectionType::Standard));

  Field3D f = are_unaligned ? toFieldAligned(f_in, "RGN_NOX") : f_in;
  Field3D v = are_unaligned ? toFieldAligned(v_in, "RGN_NOX") : v_in;
  Field3D wave_speed =
      are_unaligned ? toFieldAligned(wave_speed_in, "RGN_NOX") : wave_speed_in;

  Coordinates* coord = f_in.getCoordinates();

  Field3D result{zeroFrom(f)};

  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    const bool is_periodic_y = mesh->periodicY(i);
    const bool is_first_y = mesh->firstY(i);
    const bool is_last_y = mesh->lastY(i);

    // Only need one guard cell, so no need to communicate fluxes Instead
    // calculate in guard cells to get fluxes consistent between processors, but
    // don't include the boundary cell. Note that this implies special handling
    // of boundaries later
    const int ys = (!is_first_y || is_periodic_y) ? mesh->ystart - 1 : mesh->ystart;
    const int ye = (!is_last_y || is_periodic_y) ? mesh->yend + 1 : mesh->yend;

    for (int j = ys; j <= ye; j++) {
      // Pre-calculate factors which multiply fluxes
#if not(BOUT_USE_METRIC_3D)
      // For right cell boundaries
      BoutReal common_factor = (coord->J(i, j) + coord->J(i, j + 1))
                               / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j + 1)));

      const BoutReal flux_factor_rc = common_factor / (coord->dy(i, j) * coord->J(i, j));
      const BoutReal flux_factor_rp =
          common_factor / (coord->dy(i, j + 1) * coord->J(i, j + 1));

      // For left cell boundaries
      common_factor = (coord->J(i, j) + coord->J(i, j - 1))
                      / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j - 1)));

      const BoutReal flux_factor_lc = common_factor / (coord->dy(i, j) * coord->J(i, j));
      const BoutReal flux_factor_lm =
          common_factor / (coord->dy(i, j - 1) * coord->J(i, j - 1));
#endif
      for (int k = 0; k < mesh->LocalNz; k++) {
#if BOUT_USE_METRIC_3D
        // For right cell boundaries
        BoutReal common_factor =
            (coord->J(i, j, k) + coord->J(i, j + 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j + 1, k)));

        BoutReal flux_factor_rc =
            common_factor / (coord->dy(i, j, k) * coord->J(i, j, k));
        BoutReal flux_factor_rp =
            common_factor / (coord->dy(i, j + 1, k) * coord->J(i, j + 1, k));

        // For left cell boundaries
        common_factor = (coord->J(i, j, k) + coord->J(i, j - 1, k))
                        / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j - 1, k)));

        BoutReal flux_factor_lc =
            common_factor / (coord->dy(i, j, k) * coord->J(i, j, k));
        BoutReal flux_factor_lm =
            common_factor / (coord->dy(i, j - 1, k) * coord->J(i, j - 1, k));
#endif

        ////////////////////////////////////////////
        // Reconstruct f at the cell faces
        // This calculates s.R and s.L for the Right and Left
        // face values on this cell

        // Reconstruct f at the cell faces
        Stencil1D s;
        s.c = f(i, j, k);
        s.m = f(i, j - 1, k);
        s.p = f(i, j + 1, k);

        cellboundary(s); // Calculate s.R and s.L

        ////////////////////////////////////////////
        // Right boundary

        // Calculate velocity at right boundary (y+1/2)
        BoutReal vpar = 0.5 * (v(i, j, k) + v(i, j + 1, k));
        BoutReal flux = NAN;

        if (is_last_y && (j == mesh->yend) && !is_periodic_y) {
          // Last point in domain

          const BoutReal bndryval = 0.5 * (s.c + s.p);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux = (s.R * vpar) + (wave_speed(i, j, k) * (s.R - bndryval));
          }
        } else {

          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k));

          if (vpar > amax) {
            // Supersonic flow out of this cell
            flux = s.R * vpar;
          } else if (vpar < -amax) {
            // Supersonic flow into this cell
            flux = 0.0;
          } else {
            // Subsonic flow, so a mix of right and left fluxes
            flux = s.R * 0.5 * (vpar + amax);
          }
        }

        result(i, j, k) += flux * flux_factor_rc;
        result(i, j + 1, k) -= flux * flux_factor_rp;

        ////////////////////////////////////////////
        // Calculate at left boundary

        vpar = 0.5 * (v(i, j, k) + v(i, j - 1, k));

        if (is_first_y && (j == mesh->ystart) && !is_periodic_y) {
          // First point in domain
          const BoutReal bndryval = 0.5 * (s.c + s.m);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux = (s.L * vpar) - (wave_speed(i, j, k) * (s.L - bndryval));
          }
        } else {

          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k));

          if (vpar < -amax) {
            // Supersonic out of this cell
            flux = s.L * vpar;
          } else if (vpar > amax) {
            // Supersonic into this cell
            flux = 0.0;
          } else {
            flux = s.L * 0.5 * (vpar - amax);
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

  Coordinates* coord = n_in.getCoordinates();

  if (v.covariant) {
    // Got a covariant vector instead
    throw BoutException("Div_f_v passed a covariant v");
  }

  Field3D result{zeroFrom(n_in)};

  Field3D vx = v.x;
  Field3D vz = v.z;
  Field3D n = n_in;

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    // Calculate velocities
    const BoutReal vU = 0.25 * (vz[i.zp()] + vz[i]) * (coord->J[i.zp()] + coord->J[i]);
    const BoutReal vD = 0.25 * (vz[i.zm()] + vz[i]) * (coord->J[i.zm()] + coord->J[i]);
    const BoutReal vL = 0.25 * (vx[i.xm()] + vx[i]) * (coord->J[i.xm()] + coord->J[i]);
    const BoutReal vR = 0.25 * (vx[i.xp()] + vx[i]) * (coord->J[i.xp()] + coord->J[i]);

    // X direction
    Stencil1D s;
    s.c = n[i];
    s.m = n[i.xm()];
    s.mm = n[i.xmm()];
    s.p = n[i.xp()];
    s.pp = n[i.xpp()];

    cellboundary(s);

    if ((i.x() == mesh->xend) && (mesh->lastX())) {
      // At right boundary in X
      if (bndry_flux) {
        BoutReal flux = NAN;
        if (vR > 0.0) {
          // Flux to boundary
          flux = vR * s.R;
        } else {
          // Flux in from boundary
          flux = vR * 0.5 * (n[i.xp()] + n[i]);
        }
        result[i] += flux / (coord->dx[i] * coord->J[i]);
        result[i.xp()] -= flux / (coord->dx[i.xp()] * coord->J[i.xp()]);
      }
    } else {
      // Not at a boundary
      if (vR > 0.0) {
        // Flux out into next cell
        const BoutReal flux = vR * s.R;
        result[i] += flux / (coord->dx[i] * coord->J[i]);
        result[i.xp()] -= flux / (coord->dx[i.xp()] * coord->J[i.xp()]);
      }
    }

    // Left side

    if ((i.x() == mesh->xstart) && (mesh->firstX())) {
      // At left boundary in X

      if (bndry_flux) {
        BoutReal flux = NAN;
        if (vL < 0.0) {
          // Flux to boundary
          flux = vL * s.L;
        } else {
          // Flux in from boundary
          flux = vL * 0.5 * (n[i.xm()] + n[i]);
        }
        result[i] -= flux / (coord->dx[i] * coord->J[i]);
        result[i.xm()] += flux / (coord->dx[i.xm()] * coord->J[i.xm()]);
      }
    } else {
      // Not at a boundary
      if (vL < 0.0) {
        const BoutReal flux = vL * s.L;
        result[i] -= flux / (coord->dx[i] * coord->J[i]);
        result[i.xm()] += flux / (coord->dx[i.xm()] * coord->J[i.xm()]);
      }
    }

    /// NOTE: Need to communicate fluxes

    // Z direction
    s.m = n[i.zm()];
    s.mm = n[i.zmm()];
    s.p = n[i.zp()];
    s.pp = n[i.zpp()];

    cellboundary(s);

    if (vU > 0.0) {
      const BoutReal flux = vU * s.R;
      result[i] += flux / (coord->J[i] * coord->dz[i]);
      result[i.zp()] -= flux / (coord->J[i.zp()] * coord->dz[i.zp()]);
    }
    if (vD < 0.0) {
      const BoutReal flux = vD * s.L;
      result[i] -= flux / (coord->J[i] * coord->dz[i]);
      result[i.zm()] += flux / (coord->J[i.zm()] * coord->dz[i.zm()]);
    }
  }

  communicateFluxes(result);

  // Y advection
  // Currently just using simple centered differences
  // so no fluxes need to be exchanged

  n = toFieldAligned(n_in, "RGN_NOX");
  Field3D vy = toFieldAligned(v.y, "RGN_NOX");

  Field3D yresult = 0.0;
  yresult.setDirectionY(YDirectionType::Aligned);

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    // Y velocities on y boundaries
    const BoutReal vU = 0.25 * (vy[i] + vy[i.yp()]) * (coord->J[i] + coord->J[i.yp()]);
    const BoutReal vD = 0.25 * (vy[i] + vy[i.ym()]) * (coord->J[i] + coord->J[i.ym()]);

    // n (advected quantity) on y boundaries
    // Note: Use unshifted n_in variable
    const BoutReal nU = 0.5 * (n[i] + n[i.yp()]);
    const BoutReal nD = 0.5 * (n[i] + n[i.ym()]);

    yresult[i] = (nU * vU - nD * vD) / (coord->J[i] * coord->dy[i]);
  }
  return result + fromFieldAligned(yresult, "RGN_NOBNDRY");
}

/*!
   * X-Z Finite Volume diffusion operator
   */
Field3D Div_Perp_Lap(const Field3D& a, const Field3D& f, CELL_LOC outloc = CELL_DEFAULT);

/// Finite volume parallel divergence
///
/// NOTE: Modified version, applies limiter to velocity and field
///       Performs better (smaller overshoots) than Div_par
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
/// @param[out] flow_ylow    Flow at the lower Y cell boundary
///                          Already includes area factor * flux
template <typename CellEdges = MC>
Field3D Div_par_mod(const Field3D& f_in, const Field3D& v_in,
                    const Field3D& wave_speed_in, Field3D& flow_ylow,
                    bool fixflux = true) {

  Coordinates* coord = f_in.getCoordinates();

  if (f_in.isFci()) {
    // Use mid-point (cell boundary) averages
    if (flow_ylow.isAllocated()) {
      flow_ylow = emptyFrom(flow_ylow);
    }

    ASSERT1(f_in.hasParallelSlices());
    ASSERT1(v_in.hasParallelSlices());

    const auto& f_up = f_in.yup();
    const auto& f_down = f_in.ydown();

    const auto& v_up = v_in.yup();
    const auto& v_down = v_in.ydown();

    Field3D result{emptyFrom(f_in)};
    BOUT_FOR(i, f_in.getRegion("RGN_NOBNDRY")) {
      const auto iyp = i.yp();
      const auto iym = i.ym();

      result[i] = (0.25 * (f_in[i] + f_up[iyp]) * (v_in[i] + v_up[iyp])
                       * (coord->J[i] + coord->J.yup()[iyp])
                       / (sqrt(coord->g_22[i]) + sqrt(coord->g_22.yup()[iyp]))
                   - 0.25 * (f_in[i] + f_down[iym]) * (v_in[i] + v_down[iym])
                         * (coord->J[i] + coord->J.ydown()[iym])
                         / (sqrt(coord->g_22[i]) + sqrt(coord->g_22.ydown()[iym])))
                  / (coord->dy[i] * coord->J[i]);
    }
    return result;
  }
  ASSERT1_FIELDS_COMPATIBLE(f_in, v_in);
  ASSERT1_FIELDS_COMPATIBLE(f_in, wave_speed_in);

  const Mesh* mesh = f_in.getMesh();

  CellEdges cellboundary;

  ASSERT2(f_in.getDirectionY() == v_in.getDirectionY());
  ASSERT2(f_in.getDirectionY() == wave_speed_in.getDirectionY());
  const bool are_unaligned =
      ((f_in.getDirectionY() == YDirectionType::Standard)
       and (v_in.getDirectionY() == YDirectionType::Standard)
       and (wave_speed_in.getDirectionY() == YDirectionType::Standard));

  const Field3D f = are_unaligned ? toFieldAligned(f_in, "RGN_NOX") : f_in;
  const Field3D v = are_unaligned ? toFieldAligned(v_in, "RGN_NOX") : v_in;
  const Field3D wave_speed =
      are_unaligned ? toFieldAligned(wave_speed_in, "RGN_NOX") : wave_speed_in;

  Field3D result{zeroFrom(f)};
  flow_ylow = zeroFrom(f);

  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    const bool is_periodic_y = mesh->periodicY(i);
    const bool is_first_y = mesh->firstY(i);
    const bool is_last_y = mesh->lastY(i);

    // Only need one guard cell, so no need to communicate fluxes Instead
    // calculate in guard cells to get fluxes consistent between processors, but
    // don't include the boundary cell. Note that this implies special handling
    // of boundaries later
    const int ys = (!is_first_y || is_periodic_y) ? mesh->ystart - 1 : mesh->ystart;
    const int ye = (!is_last_y || is_periodic_y) ? mesh->yend + 1 : mesh->yend;

    for (int j = ys; j <= ye; j++) {
      // Pre-calculate factors which multiply fluxes
#if not(BOUT_USE_METRIC_3D)
      // For right cell boundaries
      const BoutReal common_factor_r =
          (coord->J(i, j) + coord->J(i, j + 1))
          / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j + 1)));

      const BoutReal flux_factor_rc =
          common_factor_r / (coord->dy(i, j) * coord->J(i, j));
      const BoutReal flux_factor_rp =
          common_factor_r / (coord->dy(i, j + 1) * coord->J(i, j + 1));

      const BoutReal area_rp =
          common_factor_r * coord->dx(i, j + 1) * coord->dz(i, j + 1);

      // For left cell boundaries
      const BoutReal common_factor_l =
          (coord->J(i, j) + coord->J(i, j - 1))
          / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j - 1)));

      const BoutReal flux_factor_lc =
          common_factor_l / (coord->dy(i, j) * coord->J(i, j));
      const BoutReal flux_factor_lm =
          common_factor_l / (coord->dy(i, j - 1) * coord->J(i, j - 1));

      const BoutReal area_lc = common_factor_l * coord->dx(i, j) * coord->dz(i, j);
#endif
      for (int k = 0; k < mesh->LocalNz; k++) {
#if BOUT_USE_METRIC_3D
        // For right cell boundaries
        const BoutReal common_factor_r =
            (coord->J(i, j, k) + coord->J(i, j + 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j + 1, k)));

        const BoutReal flux_factor_rc =
            common_factor_r / (coord->dy(i, j, k) * coord->J(i, j, k));
        const BoutReal flux_factor_rp =
            common_factor_r / (coord->dy(i, j + 1, k) * coord->J(i, j + 1, k));

        const BoutReal area_rp =
            common_factor_r * coord->dx(i, j + 1, k) * coord->dz(i, j + 1, k);

        // For left cell boundaries
        const BoutReal common_factor_l =
            (coord->J(i, j, k) + coord->J(i, j - 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j - 1, k)));

        const BoutReal flux_factor_lc =
            common_factor_l / (coord->dy(i, j, k) * coord->J(i, j, k));
        const BoutReal flux_factor_lm =
            common_factor_l / (coord->dy(i, j - 1, k) * coord->J(i, j - 1, k));

        const BoutReal area_lc =
            common_factor_l * coord->dx(i, j, k) * coord->dz(i, j, k);
#endif

        ////////////////////////////////////////////
        // Reconstruct f at the cell faces
        // This calculates s.R and s.L for the Right and Left
        // face values on this cell

        // Reconstruct f at the cell faces
        // TODO(peter): We can remove this #ifdef guard after switching to C++20
#if __cpp_designated_initializers >= 201707L
        Stencil1D s{.c = f(i, j, k), .m = f(i, j - 1, k), .p = f(i, j + 1, k)};
#else
        Stencil1D s{f(i, j, k), f(i, j - 1, k), f(i, j + 1, k), BoutNaN,
                    BoutNaN,    BoutNaN,        BoutNaN};
#endif
        cellboundary(s); // Calculate s.R and s.L

        ////////////////////////////////////////////
        // Reconstruct v at the cell faces
        // TODO(peter): We can remove this #ifdef guard after switching to C++20
#if __cpp_designated_initializers >= 201707L
        Stencil1D sv{.c = v(i, j, k), .m = v(i, j - 1, k), .p = v(i, j + 1, k)};
#else
        Stencil1D sv{v(i, j, k), v(i, j - 1, k), v(i, j + 1, k), BoutNaN,
                     BoutNaN,    BoutNaN,        BoutNaN};
#endif
        cellboundary(sv); // Calculate sv.R and sv.L

        ////////////////////////////////////////////
        // Right boundary

        BoutReal flux = BoutNaN;

        if (is_last_y && (j == mesh->yend) && !is_periodic_y) {
          // Last point in domain

          // Calculate velocity at right boundary (y+1/2)
          const BoutReal vpar = 0.5 * (v(i, j, k) + v(i, j + 1, k));

          const BoutReal bndryval = 0.5 * (s.c + s.p);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux = (s.R * vpar) + (wave_speed(i, j, k) * (s.R - bndryval));
          }

        } else {
          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k),
                                        fabs(v(i, j, k)), fabs(v(i, j + 1, k)));

          flux = s.R * 0.5 * (sv.R + amax);
        }

        result(i, j, k) += flux * flux_factor_rc;
        result(i, j + 1, k) -= flux * flux_factor_rp;

        flow_ylow(i, j + 1, k) += flux * area_rp;

        ////////////////////////////////////////////
        // Calculate at left boundary

        if (is_first_y && (j == mesh->ystart) && !is_periodic_y) {
          // First point in domain
          const BoutReal bndryval = 0.5 * (s.c + s.m);
          const BoutReal vpar = 0.5 * (v(i, j, k) + v(i, j - 1, k));
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux = (s.L * vpar) - (wave_speed(i, j, k) * (s.L - bndryval));
          }
        } else {

          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k),
                                        fabs(v(i, j, k)), fabs(v(i, j - 1, k)));

          flux = s.L * 0.5 * (sv.L - amax);
        }

        result(i, j, k) -= flux * flux_factor_lc;
        result(i, j - 1, k) += flux * flux_factor_lm;

        flow_ylow(i, j, k) += flux * area_lc;
      }
    }
  }
  if (are_unaligned) {
    flow_ylow = fromFieldAligned(flow_ylow, "RGN_NOBNDRY");
  }
  return are_unaligned ? fromFieldAligned(result, "RGN_NOBNDRY") : result;
}

/// This operator calculates Div_par(f v v)
/// It is used primarily (only?) in the parallel momentum equation.
///
/// This operator is used rather than Div(f fv) so that the values of
/// f and v are consistent with other advection equations: The product
/// fv is not interpolated to cell boundaries.
template <typename CellEdges = MC>
Field3D Div_par_fvv(const Field3D& f_in, const Field3D& v_in,
                    const Field3D& wave_speed_in, bool fixflux = true) {
  ASSERT1_FIELDS_COMPATIBLE(f_in, v_in);
  const Mesh* mesh = f_in.getMesh();
  const Coordinates* coord = f_in.getCoordinates();
  CellEdges cellboundary;

  if (f_in.isFci()) {
    // FCI version, using yup/down fields
    ASSERT1(f_in.hasParallelSlices());
    ASSERT1(v_in.hasParallelSlices());

    const auto& B = coord->Bxy;
    const auto& B_up = coord->Bxy.yup();
    const auto& B_down = coord->Bxy.ydown();

    const auto& f_up = f_in.yup();
    const auto& f_down = f_in.ydown();

    const auto& v_up = v_in.yup();
    const auto& v_down = v_in.ydown();

    const auto& g_22 = coord->g_22;
    const auto& dy = coord->dy;

    Field3D result{emptyFrom(f_in)};
    BOUT_FOR(i, f_in.getRegion("RGN_NOBNDRY")) {
      const auto iyp = i.yp();
      const auto iym = i.ym();

      // Maximum local wave speed
      const BoutReal amax =
          BOUTMAX(wave_speed_in[i], fabs(v_in[i]), fabs(v_up[iyp]), fabs(v_down[iym]));

      const BoutReal term = (f_up[iyp] * v_up[iyp] * v_up[iyp] / B_up[iyp])
                            - (f_down[iym] * v_down[iym] * v_down[iym] / B_down[iym]);

      // Penalty terms. This implementation is very dissipative.
      BoutReal penalty =
          (amax * (f_in[i] * v_in[i] - f_up[iyp] * v_up[iyp]) / (B[i] + B_up[iyp]))
          + (amax * (f_in[i] * v_in[i] - f_down[iym] * v_down[iym])
             / (B[i] + B_down[iym]));

      if (fabs(penalty) > fabs(term) and penalty * v_in[i] > 0) {
        if (term * penalty > 0) {
          penalty = term;
        } else {
          penalty = -term;
        }
      }

      result[i] = B[i] * (term + penalty) / (2 * dy[i] * sqrt(g_22[i]));

#if CHECK > 0
      if (!std::isfinite(result[i])) {
        throw BoutException("Non-finite value in Div_par_fvv at {}\n"
                            "fup {} vup {} fdown {} vdown {} amax {}\n",
                            "B {} Bup {} Bdown {} dy {} sqrt(g_22} {}", i, f_up[i],
                            v_up[i], f_down[i], v_down[i], amax, B[i], B_up[i], B_down[i],
                            dy[i], sqrt(g_22[i]));
      }
#endif
    }
    return result;
  }

  ASSERT1(areFieldsCompatible(f_in, wave_speed_in));

  /// Ensure that f, v and wave_speed are field aligned
  Field3D f = toFieldAligned(f_in, "RGN_NOX");
  Field3D v = toFieldAligned(v_in, "RGN_NOX");
  Field3D wave_speed = toFieldAligned(wave_speed_in, "RGN_NOX");

  Field3D result{zeroFrom(f)};

  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    const bool is_periodic_y = mesh->periodicY(i);
    const bool is_first_y = mesh->firstY(i);
    const bool is_last_y = mesh->lastY(i);

    // Only need one guard cell, so no need to communicate fluxes Instead
    // calculate in guard cells to get fluxes consistent between processors, but
    // don't include the boundary cell. Note that this implies special handling
    // of boundaries later
    const int ys = (!is_first_y || is_periodic_y) ? mesh->ystart - 1 : mesh->ystart;
    const int ye = (!is_last_y || is_periodic_y) ? mesh->yend + 1 : mesh->yend;

    for (int j = ys; j <= ye; j++) {
      // Pre-calculate factors which multiply fluxes

      for (int k = 0; k < mesh->LocalNz; k++) {
        // For right cell boundaries
        const BoutReal common_factor_r =
            (coord->J(i, j, k) + coord->J(i, j + 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j + 1, k)));

        const BoutReal flux_factor_rc =
            common_factor_r / (coord->dy(i, j, k) * coord->J(i, j, k));
        const BoutReal flux_factor_rp =
            common_factor_r / (coord->dy(i, j + 1, k) * coord->J(i, j + 1, k));

        // For left cell boundaries
        const BoutReal common_factor_l =
            (coord->J(i, j, k) + coord->J(i, j - 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j - 1, k)));

        const BoutReal flux_factor_lc =
            common_factor_l / (coord->dy(i, j, k) * coord->J(i, j, k));
        const BoutReal flux_factor_lm =
            common_factor_l / (coord->dy(i, j - 1, k) * coord->J(i, j - 1, k));

        ////////////////////////////////////////////
        // Reconstruct f at the cell faces
        // This calculates s.R and s.L for the Right and Left
        // face values on this cell

        // Reconstruct f at the cell faces
#if __cpp_designated_initializers >= 201707L
        Stencil1D s{.c = f(i, j, k), .m = f(i, j - 1, k), .p = f(i, j + 1, k)};
#else
        Stencil1D s{f(i, j, k), f(i, j - 1, k), f(i, j + 1, k), BoutNaN,
                    BoutNaN,    BoutNaN,        BoutNaN};
#endif
        cellboundary(s); // Calculate s.R and s.L

        ////////////////////////////////////////////
        // Reconstruct v at the cell faces
        // TODO(peter): We can remove this #ifdef guard after switching to C++20
#if __cpp_designated_initializers >= 201707L
        Stencil1D sv{.c = v(i, j, k), .m = v(i, j - 1, k), .p = v(i, j + 1, k)};
#else
        Stencil1D sv{v(i, j, k), v(i, j - 1, k), v(i, j + 1, k), BoutNaN,
                     BoutNaN,    BoutNaN,        BoutNaN};
#endif
        cellboundary(sv);

        ////////////////////////////////////////////
        // Right boundary

        // Calculate velocity at right boundary (y+1/2)
        const BoutReal v_mid_r = 0.5 * (sv.c + sv.p);
        // And mid-point density at right boundary
        const BoutReal n_mid_r = 0.5 * (s.c + s.p);
        BoutReal flux = NAN;

        if (mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
          // Last point in domain

          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = n_mid_r * v_mid_r * v_mid_r;
          } else {
            // Add flux due to difference in boundary values
            flux = (s.R * sv.R * sv.R) // Use right cell edge values
                   + (BOUTMAX(wave_speed(i, j, k), fabs(sv.c), fabs(sv.p)) * n_mid_r
                      * (sv.R - v_mid_r)); // Damp differences in velocity, not flux
          }
        } else {
          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k),
                                        fabs(sv.c), fabs(sv.p));

          flux = s.R * 0.5 * (sv.R + amax) * sv.R;
        }

        result(i, j, k) += flux * flux_factor_rc;
        result(i, j + 1, k) -= flux * flux_factor_rp;

        ////////////////////////////////////////////
        // Calculate at left boundary

        const BoutReal v_mid_l = 0.5 * (sv.c + sv.m);
        const BoutReal n_mid_l = 0.5 * (s.c + s.m);

        if (mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
          // First point in domain
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = n_mid_l * v_mid_l * v_mid_l;
          } else {
            // Add flux due to difference in boundary values
            flux = (s.L * sv.L * sv.L)
                   - (BOUTMAX(wave_speed(i, j, k), fabs(sv.c), fabs(sv.m)) * n_mid_l
                      * (sv.L - v_mid_l));
          }
        } else {
          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k),
                                        fabs(sv.c), fabs(sv.m));

          flux = s.L * 0.5 * (sv.L - amax) * sv.L;
        }

        result(i, j, k) -= flux * flux_factor_lc;
        result(i, j - 1, k) += flux * flux_factor_lm;
      }
    }
  }
  return fromFieldAligned(result, "RGN_NOBNDRY");
}
} // namespace FV
#endif // BOUT_FV_OPS_H
