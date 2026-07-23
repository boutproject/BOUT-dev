/*
  Finite-volume discretisation methods. Flux-conservative form
 */

#ifndef BOUT_FV_OPS_IMPL_H
#define BOUT_FV_OPS_IMPL_H

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/build_defines.hxx"
#include "bout/coordinates.hxx"
#include "bout/field.hxx"
#include "bout/field3d.hxx"
#include "bout/fv_ops.hxx" // NOLINT(unused-includes, misc-include-cleaner)
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/output_bout_types.hxx" // NOLINT(unused-includes, misc-include-cleaner)
#include "bout/region.hxx"
#include "bout/utils.hxx"
#include "bout/vector2d.hxx"

#include <cmath>

namespace FV {
/*!
   * Stencil used for Finite Volume calculations
   * which includes cell face values L and R
   */
struct Stencil1D {
  /// Cell centre values
  BoutReal c{};
  BoutReal m{};
  BoutReal p{};
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
// NB: The templates need to be explicitly instantiated in fv_ops.cxx
struct Upwind {
  void operator()(Stencil1D& n) { n.L = n.R = n.c; }
};

/*!
   * Fromm method
   */
// NB: The templates need to be explicitly instantiated in fv_ops.cxx
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
// NB: The templates need to be explicitly instantiated in fv_ops.cxx
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

    if (std::abs(a) < std::abs(b)) {
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
// NB: The templates need to be explicitly instantiated in fv_ops.cxx
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
    return SIGN(a) * BOUTMIN(std::abs(a), std::abs(b), std::abs(c));
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
// NB: The templates need to be explicitly instantiated in fv_ops.cxx
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
      const BoutReal abs_gL = std::abs(gL);
      const BoutReal abs_gR = std::abs(gR);
      const BoutReal half_slope =
          sign * BOUTMAX(BOUTMIN(abs_gL, 0.5 * abs_gR), BOUTMIN(abs_gR, 0.5 * abs_gL));
      n.L = n.c - half_slope;
      n.R = n.c + half_slope;
    }
  }
};

/*!
   * Symmetric Van Albada second order slope limiter
   *
   * Uses a smooth (differentiable) approximation to `max(a*b, 0)` to avoid
   * introducing a kink at extrema, which can be helpful for nonlinear solvers
   * and finite-difference Jacobian calculations.
   *
   * The limited slope is calculated from the left and right differences
   * `dl = c - m` and `dr = p - c` as
   *
   *   slope = (pos(dl*dr) * (dl + dr)) / (dl^2 + dr^2)
   *
   * where `pos(x)` is a smooth approximation to `max(x, 0)`.
   */
// NB: The templates need to be explicitly instantiated in fv_ops.cxx
struct VanAlbada {
  void operator()(Stencil1D& n) {
    const BoutReal dl = n.c - n.m;
    const BoutReal dr = n.p - n.c;

    const BoutReal denom = (dl * dl) + (dr * dr);

    // Smoothness parameters:
    // - keep division well-defined when dl=dr=0
    // - provide a differentiable approximation to max(dl*dr, 0)
    const BoutReal eps = (1e-12 * denom) + 1e-30;

    const BoutReal ab = dl * dr;
    const BoutReal ab_pos = 0.5 * (ab + sqrt((ab * ab) + (eps * eps)));

    const BoutReal slope = (ab_pos * (dl + dr)) / (denom + eps);

    n.L = n.c - (0.5 * slope);
    n.R = n.c + (0.5 * slope);
  }
};

/*!
   * WENO3-JS (Jiang-Shu) reconstruction to cell faces
   *
   * This is a third-order essentially non-oscillatory reconstruction using two
   * candidate second-order polynomials and smoothness-weighted blending.
   *
   * Unlike TVD slope limiters (e.g. ``MC``), WENO reconstruction is generally
   * smooth (differentiable) for all inputs, but it does not enforce strict
   * monotonicity.
   *
   * Uses only the three-point stencil (`m`, `c`, `p`), so it is a drop-in
   * replacement anywhere `Stencil1D` is populated with those values.
   */
// NB: The templates need to be explicitly instantiated in fv_ops.cxx
struct WENO3 {
  void operator()(Stencil1D& n) {
    // Right face (between c and p): value from cell c (left state at i+1/2)
    const BoutReal p0_r = 0.5 * (-n.m + 3.0 * n.c);
    const BoutReal p1_r = 0.5 * (n.c + n.p);

    const BoutReal beta0_r = SQ(n.c - n.m);
    const BoutReal beta1_r = SQ(n.p - n.c);

    // Left face (between m and c): value from cell c (right state at i-1/2)
    const BoutReal p0_l = 0.5 * (-n.p + 3.0 * n.c);
    const BoutReal p1_l = 0.5 * (n.m + n.c);

    const BoutReal beta0_l = beta1_r;
    const BoutReal beta1_l = beta0_r;

    // Smoothness parameter (scaled to local variation)
    const BoutReal eps = (1e-12 * (beta0_r + beta1_r)) + 1e-30;

    // Linear weights for WENO3-JS
    constexpr BoutReal d0 = 1.0 / 3.0;
    constexpr BoutReal d1 = 2.0 / 3.0;

    // Right face weights
    const BoutReal a0_r = d0 / SQ(eps + beta0_r);
    const BoutReal a1_r = d1 / SQ(eps + beta1_r);
    const BoutReal wsum_r = a0_r + a1_r;
    const BoutReal w0_r = a0_r / wsum_r;
    const BoutReal w1_r = a1_r / wsum_r;

    // Left face weights (mirrored)
    const BoutReal a0_l = d0 / SQ(eps + beta0_l);
    const BoutReal a1_l = d1 / SQ(eps + beta1_l);
    const BoutReal wsum_l = a0_l + a1_l;
    const BoutReal w0_l = a0_l / wsum_l;
    const BoutReal w1_l = a1_l / wsum_l;

    n.R = (w0_r * p0_r) + (w1_r * p1_r);
    n.L = (w0_l * p0_l) + (w1_l * p1_l);
  }
};

template <class T>
struct Slices {
  T c;
  T up;
  T down;

  Slices(bool use_slices, const T& field)
      : c(use_slices ? field : toFieldAligned(field)), up(use_slices ? field.yup() : c),
        down(use_slices ? field.ydown() : c) {}
};

template <class T>
Slices<T> makeslices(bool use_slices, const T& field) {
  return Slices<T>(use_slices, field);
}

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
template <typename CellEdges>
Field3D Div_par(const Field3D& f_in, const Field3D& v_in, const Field3D& wave_speed_in,
                bool fixflux) {

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
      const BoutReal area_r = coord->cell_area_yhigh()(i, j);
      const BoutReal flux_factor_rc = area_r / coord->cell_volume()(i, j);
      const BoutReal flux_factor_rp = area_r / coord->cell_volume()(i, j + 1);
      // For left cell boundaries
      const BoutReal area_l = coord->cell_area_ylow()(i, j);
      const BoutReal flux_factor_lc = area_l / coord->cell_volume()(i, j);
      const BoutReal flux_factor_lm = area_l / coord->cell_volume()(i, j - 1);
#endif
      for (int k = mesh->zstart; k <= mesh->zend; k++) {
#if BOUT_USE_METRIC_3D
        // For right cell boundaries
        const BoutReal area_r = coord->cell_area_yhigh()(i, j, k);
        const BoutReal flux_factor_rc = area_r / coord->cell_volume()(i, j, k);
        const BoutReal flux_factor_rp = area_r / coord->cell_volume()(i, j + 1, k);
        // For left cell boundaries
        const BoutReal area_l = coord->cell_area_ylow()(i, j, k);
        const BoutReal flux_factor_lc = area_l / coord->cell_volume()(i, j, k);
        const BoutReal flux_factor_lm = area_l / coord->cell_volume()(i, j - 1, k);
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
        BoutReal flux = BoutNaN;

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
template <typename CellEdges>
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
        BoutReal flux = BoutNaN;
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
        BoutReal flux = BoutNaN;
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
template <typename CellEdges>
Field3D Div_par_mod(const Field3D& f_in, const Field3D& v_in,
                    const Field3D& wave_speed_in, Field3D& flow_ylow, bool fixflux) {

  Coordinates* coord = f_in.getCoordinates();
  ASSERT1_FIELDS_COMPATIBLE(f_in, v_in);

  if (f_in.isFci()) {
    // Use mid-point (cell boundary) averages

    ASSERT1(f_in.hasParallelSlices());
    ASSERT1(v_in.hasParallelSlices());

    const auto& f_up = f_in.yup();
    const auto& f_down = f_in.ydown();

    const auto& v_up = v_in.yup();
    const auto& v_down = v_in.ydown();

    Field3D result{emptyFrom(f_in)};
    flow_ylow = zeroFrom(f_in);
    BOUT_FOR(i, f_in.getRegion("RGN_NOBNDRY")) {
      const auto iyp = i.yp();
      const auto iym = i.ym();

      result[i] = (0.25 * (f_in[i] + f_up[iyp]) * (v_in[i] + v_up[iyp])
                       * coord->cell_area_yhigh()[i]
                   - 0.25 * (f_in[i] + f_down[iym]) * (v_in[i] + v_down[iym])
                         * coord->cell_area_ylow()[i])
                  / coord->cell_volume()[i];
    }
    return result;
  }
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
      const BoutReal area_rp = coord->cell_area_yhigh()(i, j);

      const BoutReal flux_factor_rc = area_rp / coord->cell_volume()(i, j);
      const BoutReal flux_factor_rp = area_rp / coord->cell_volume()(i, j + 1);

      // For left cell boundaries
      const BoutReal area_lc = coord->cell_area_ylow()(i, j);

      const BoutReal flux_factor_lc = area_lc / coord->cell_volume()(i, j);
      const BoutReal flux_factor_lm = area_lc / coord->cell_volume()(i, j - 1);
#endif
      for (int k = 0; k < mesh->LocalNz; k++) {
#if BOUT_USE_METRIC_3D
        // For right cell boundaries
        const BoutReal area_rp = coord->cell_area_yhigh()(i, j, k);

        const BoutReal flux_factor_rc = area_rp / coord->cell_volume()(i, j, k);
        const BoutReal flux_factor_rp = area_rp / coord->cell_volume()(i, j + 1, k);

        // For left cell boundaries
        const BoutReal area_lc = coord->cell_area_ylow()(i, j, k);

        const BoutReal flux_factor_lc = area_lc / coord->cell_volume()(i, j, k);
        const BoutReal flux_factor_lm = area_lc / coord->cell_volume()(i, j - 1, k);
#endif

        ////////////////////////////////////////////
        // Reconstruct f at the cell faces
        // This calculates s.R and s.L for the Right and Left
        // face values on this cell

        // Reconstruct f at the cell faces
        Stencil1D s{.c = f(i, j, k), .m = f(i, j - 1, k), .p = f(i, j + 1, k)};
        cellboundary(s); // Calculate s.R and s.L

        ////////////////////////////////////////////
        // Reconstruct v at the cell faces
        Stencil1D sv{.c = v(i, j, k), .m = v(i, j - 1, k), .p = v(i, j + 1, k)};
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
                                        std::abs(v(i, j, k)), std::abs(v(i, j + 1, k)));

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
                                        std::abs(v(i, j, k)), std::abs(v(i, j - 1, k)));

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
template <typename CellEdges>
Field3D Div_par_fvv(const Field3D& f_in, const Field3D& v_in,
                    const Field3D& wave_speed_in, bool fixflux) {
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
      const BoutReal amax = BOUTMAX(wave_speed_in[i], std::abs(v_in[i]),
                                    std::abs(v_up[iyp]), std::abs(v_down[iym]));

      const BoutReal term = (f_up[iyp] * v_up[iyp] * v_up[iyp] / B_up[iyp])
                            - (f_down[iym] * v_down[iym] * v_down[iym] / B_down[iym]);

      // Penalty terms. This implementation is very dissipative.
      BoutReal penalty =
          (amax * (f_in[i] * v_in[i] - f_up[iyp] * v_up[iyp]) / (B[i] + B_up[iyp]))
          + (amax * (f_in[i] * v_in[i] - f_down[iym] * v_down[iym])
             / (B[i] + B_down[iym]));

      if (std::abs(penalty) > std::abs(term) and penalty * v_in[i] > 0) {
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
        const BoutReal area_r = coord->cell_area_yhigh()(i, j, k);

        const BoutReal flux_factor_rc = area_r / coord->cell_volume()(i, j, k);
        const BoutReal flux_factor_rp = area_r / coord->cell_volume()(i, j + 1, k);

        // For left cell boundaries
        const BoutReal area_l = coord->cell_area_ylow()(i, j, k);

        const BoutReal flux_factor_lc = area_l / coord->cell_volume()(i, j, k);
        const BoutReal flux_factor_lm = area_l / coord->cell_volume()(i, j - 1, k);

        ////////////////////////////////////////////
        // Reconstruct f at the cell faces
        // This calculates s.R and s.L for the Right and Left
        // face values on this cell

        // Reconstruct f at the cell faces
        Stencil1D s{.c = f(i, j, k), .m = f(i, j - 1, k), .p = f(i, j + 1, k)};
        cellboundary(s); // Calculate s.R and s.L

        ////////////////////////////////////////////
        // Reconstruct v at the cell faces
        Stencil1D sv{.c = v(i, j, k), .m = v(i, j - 1, k), .p = v(i, j + 1, k)};
        cellboundary(sv);

        ////////////////////////////////////////////
        // Right boundary

        // Calculate velocity at right boundary (y+1/2)
        const BoutReal v_mid_r = 0.5 * (sv.c + sv.p);
        // And mid-point density at right boundary
        const BoutReal n_mid_r = 0.5 * (s.c + s.p);
        BoutReal flux = BoutNaN;

        if (mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
          // Last point in domain

          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = n_mid_r * v_mid_r * v_mid_r;
          } else {
            // Add flux due to difference in boundary values
            flux =
                (s.R * sv.R * sv.R) // Use right cell edge values
                + (BOUTMAX(wave_speed(i, j, k), std::abs(sv.c), std::abs(sv.p)) * n_mid_r
                   * (sv.R - v_mid_r)); // Damp differences in velocity, not flux
          }
        } else {
          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k),
                                        std::abs(sv.c), std::abs(sv.p));

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
                   - (BOUTMAX(wave_speed(i, j, k), std::abs(sv.c), std::abs(sv.m))
                      * n_mid_l * (sv.L - v_mid_l));
          }
        } else {
          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k),
                                        std::abs(sv.c), std::abs(sv.m));

          flux = s.L * 0.5 * (sv.L - amax) * sv.L;
        }

        result(i, j, k) -= flux * flux_factor_lc;
        result(i, j - 1, k) += flux * flux_factor_lm;
      }
    }
  }
  return fromFieldAligned(result, "RGN_NOBNDRY");
}

// Calculates viscous heating due to numerical momentum fluxes
// and flow of kinetic energy (in flow_ylow)
template <typename CellEdges>
Field3D Div_par_fvv_heating(const Field3D& f_in, const Field3D& v_in,
                            const Field3D& wave_speed_in, Field3D& flow_ylow,
                            bool fixflux) {

  ASSERT1(areFieldsCompatible(f_in, v_in));
  ASSERT1(areFieldsCompatible(f_in, wave_speed_in));

  const Mesh* mesh = f_in.getMesh();
  Coordinates* coord = f_in.getCoordinates();
  CellEdges cellboundary;

  if (f_in.isFci()) {
    // FCI version, using yup/down fields
    ASSERT1(f_in.hasParallelSlices());
    ASSERT1(v_in.hasParallelSlices());

    const auto B = coord->Bxy;
    const auto B_up = coord->Bxy.yup();
    const auto B_down = coord->Bxy.ydown();

    const auto& f_up = f_in.yup();
    const auto& f_down = f_in.ydown();

    const auto& v_up = v_in.yup();
    const auto& v_down = v_in.ydown();

    const auto g_22 = coord->g_22;
    const auto dy = coord->dy;

    Field3D result{emptyFrom(f_in)};
    flow_ylow = zeroFrom(f_in);

    BOUT_FOR(i, f_in.getRegion("RGN_NOBNDRY")) {
      const auto iyp = i.yp();
      const auto iym = i.ym();

      //Maximum local wave speed
      const BoutReal amax = BOUTMAX(wave_speed_in[i], std::abs(v_in[i]),
                                    std::abs(v_up[iyp]), std::abs(v_down[iym]));

      result[i] =
          B[i]
          * ((f_up[iyp] * v_up[iyp] * v_up[iyp] / B_up[iyp])
             - (f_down[iym] * v_down[iym] * v_down[iym] / B_down[iym])
             // Penalty terms. This implementation is very dissipative.
             // Note: This version adds a viscosity that damps gradients of velocity
             + amax * (f_in[i] + f_up[iyp]) * (v_in[i] - v_up[iyp]) / (B[i] + B_up[iyp])
             + amax * (f_in[i] + f_down[iym]) * (v_in[i] - v_down[iym])
                   / (B[i] + B_down[iym]))
          / (2 * dy[i] * sqrt(g_22[i]));

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

  /// Ensure that f, v and wave_speed are field aligned
  Field3D f = toFieldAligned(f_in, "RGN_NOX");
  Field3D v = toFieldAligned(v_in, "RGN_NOX");
  Field3D wave_speed = toFieldAligned(wave_speed_in, "RGN_NOX");

  // result and flow_ylow are field-aligned.
  // Will be converted to non-aligned before return.
  Field3D result{zeroFrom(f)};
  flow_ylow = zeroFrom(f);

  // Only need one guard cell, so no need to communicate fluxes
  // Instead calculate in guard cells to preserve fluxes
  int ys = mesh->ystart - 1;
  int ye = mesh->yend + 1;

  for (int i = mesh->xstart; i <= mesh->xend; i++) {

    if (!mesh->firstY(i) || mesh->periodicY(i)) {
      // Calculate in guard cell to get fluxes consistent between processors
      ys = mesh->ystart - 1;
    } else {
      // Don't include the boundary cell. Note that this implies special
      // handling of boundaries later
      ys = mesh->ystart;
    }

    if (!mesh->lastY(i) || mesh->periodicY(i)) {
      // Calculate in guard cells
      ye = mesh->yend + 1;
    } else {
      // Not in boundary cells
      ye = mesh->yend;
    }

    for (int j = ys; j <= ye; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        // Pre-calculate factors which multiply fluxes
        // Note: In 3D metric geometries these quantities can depend on (i,j,k),
        // so calculate inside the k loop.

        // For right cell boundaries
        BoutReal common_factor =
            (coord->J(i, j, k) + coord->J(i, j + 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j + 1, k)));

        const BoutReal flux_factor_rc =
            common_factor / (coord->dy(i, j, k) * coord->J(i, j, k));
        const BoutReal area_rp =
            common_factor * coord->dx(i, j + 1, k) * coord->dz(i, j + 1, k);

        // For left cell boundaries
        common_factor = (coord->J(i, j, k) + coord->J(i, j - 1, k))
                        / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j - 1, k)));

        const BoutReal flux_factor_lc =
            common_factor / (coord->dy(i, j, k) * coord->J(i, j, k));
        const BoutReal area_lc = common_factor * coord->dx(i, j, k) * coord->dz(i, j, k);

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

        // Reconstruct v at the cell faces
        Stencil1D sv;
        sv.c = v(i, j, k);
        sv.m = v(i, j - 1, k);
        sv.p = v(i, j + 1, k);

        cellboundary(sv);

        ////////////////////////////////////////////
        // Right boundary

        // Calculate velocity at right boundary (y+1/2)
        BoutReal v_mid = 0.5 * (sv.c + sv.p);
        // And mid-point density at right boundary
        BoutReal n_mid = 0.5 * (s.c + s.p);

        if (mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
          // Last point in domain

          // Expected loss of kinetic energy into boundary
          // This is used in the sheath boundary condition to calculate
          // energy losses.
          const BoutReal expected_ke = 0.5 * n_mid * v_mid * v_mid * v_mid;

          BoutReal flux_mom = BoutNaN;
          if (fixflux) {
            // Mid-point consistent with boundary conditions
            // but kinetic energy loss will not match expected
            // -> Adjust energy balance in pressure equation
            flux_mom = n_mid * v_mid * v_mid;
          } else {
            flux_mom = (s.R * sv.R * sv.R)
                       + (BOUTMAX(wave_speed(i, j, k), std::abs(sv.c), std::abs(sv.p))
                          * (s.R * sv.R - n_mid * v_mid));
          }

          // Assume that particle flux is fixed to boundary value
          const BoutReal flux_part = n_mid * v_mid;

          // d/dt(1/2 m n v^2) = v * d/dt(mnv) - 1/2 m v^2 * dn/dt
          const BoutReal actual_ke = (sv.c * flux_mom) - (0.5 * sv.c * sv.c * flux_part);

          // Note: If the actual loss was higher than expected, then
          //       plasma heating is needed to compensate
          result(i, j, k) += (actual_ke - expected_ke) * flux_factor_rc;

          // Final flow through boundary is the expected value
          flow_ylow(i, j + 1, k) += expected_ke * area_rp; //expected_ke * area_rp;

        } else {
          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k),
                                        std::abs(sv.c), std::abs(sv.p));

          // Viscous heating due to relaxation of velocity towards midpoint
          result(i, j, k) +=
              (amax + 0.5 * sv.R) * s.R * (sv.c - sv.p) * (sv.R - v_mid) * flux_factor_rc;

          // Kinetic energy flow into next cell.
          // Note: Different from flow out of this cell; the difference
          //       is in the viscous heating.
          const BoutReal flux_part = s.R * 0.5 * (sv.R + amax);
          const BoutReal flux_mom = flux_part * sv.R;

          flow_ylow(i, j + 1, k) +=
              (sv.p * flux_mom - 0.5 * SQ(sv.p) * flux_part) * area_rp;
        }

        ////////////////////////////////////////////
        // Calculate at left boundary

        v_mid = 0.5 * (sv.c + sv.m);
        n_mid = 0.5 * (s.c + s.m);

        // Expected KE loss. Note minus sign because negative v into boundary
        const BoutReal expected_ke = -0.5 * n_mid * v_mid * v_mid * v_mid;

        if (mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
          // First point in domain
          BoutReal flux_mom = BoutNaN;
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux_mom = n_mid * v_mid * v_mid;
          } else {
            // Add flux due to difference in boundary values
            flux_mom = (s.L * sv.L * sv.L)
                       - (BOUTMAX(wave_speed(i, j, k), std::abs(sv.c), std::abs(sv.m))
                          * (s.L * sv.L - n_mid * v_mid));
          }

          // Assume that density flux is fixed to boundary value
          const BoutReal flux_part = n_mid * v_mid;

          // d/dt(1/2 m n v^2) = v * d/dt(mnv) - 1/2 m v^2 * dn/dt
          const BoutReal actual_ke = (-sv.c * flux_mom) + (0.5 * sv.c * sv.c * flux_part);

          result(i, j, k) += (actual_ke - expected_ke) * flux_factor_lc;

          flow_ylow(i, j, k) -= expected_ke * area_lc;
        } else {
          // Maximum wave speed in the two cells
          const BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k),
                                        std::abs(sv.c), std::abs(sv.m));

          // Viscous heating due to relaxation
          result(i, j, k) +=
              (amax - 0.5 * sv.L) * s.L * (sv.c - sv.m) * (sv.L - v_mid) * flux_factor_lc;

          // Kinetic energy flow into this cell.
          // Note: Different from flow out of left cell; the difference
          //       is in the viscous heating.
          const BoutReal flux_part = s.L * 0.5 * (sv.L - amax);
          const BoutReal flux_mom = flux_part * sv.L;

          flow_ylow(i, j, k) += (sv.c * flux_mom - 0.5 * SQ(sv.c) * flux_part) * area_lc;
        }
      }
    }
  }
  flow_ylow = fromFieldAligned(flow_ylow, "RGN_NOBNDRY");
  return fromFieldAligned(result, "RGN_NOBNDRY");
}

/// Div ( a g Grad_perp(f) )  -- Perpendicular gradient-driven advection
///
/// This version uses a slope limiter to calculate cell edge values of g in X,
/// the advects the upwind cell edge.
///
/// 1st order upwinding is used in Y.
template <typename CellEdges>
Field3D Div_a_Grad_perp_limit(const Field3D& a, const Field3D& g, const Field3D& f) {
  ASSERT2(a.getLocation() == f.getLocation());

  Mesh* mesh = a.getMesh();

  // Requires at least 2 communication guard cells in X, 1 in Y
  ASSERT1(mesh->xstart >= 2);
  ASSERT1(mesh->ystart >= 1);

  CellEdges cellboundary;

  Field3D result{zeroFrom(f)};

  Coordinates* coord = f.getCoordinates();

  // Flux in x

  for (int i = mesh->xstart - 1; i <= mesh->xend; i++) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        // Calculate flux from i to i+1

        const BoutReal gradient = f(i + 1, j, k) - f(i, j, k);

        // Mid-point average boundary value of 'a'
        const BoutReal aedge = 0.5 * (a(i + 1, j, k) + a(i, j, k));
        BoutReal gedge = BoutNaN;
        if (((i == mesh->xstart - 1) and mesh->firstX())
            or ((i == mesh->xend) and mesh->lastX())) {
          // Mid-point average boundary value of 'g'
          gedge = 0.5 * (g(i + 1, j, k) + g(i, j, k));
        } else if (gradient > 0) {
          // Flux is from (i+1) to (i)
          // Reconstruct `g` at left of (i+1, j, k)

          Stencil1D sg;
          sg.m = g(i, j, k);
          sg.c = g(i + 1, j, k);
          sg.p = g(i + 2, j, k);
          cellboundary(sg); // Calculate sg.R and sg.L

          gedge = sg.L;
        } else {
          // Flux is from (i) to (i+1)
          // Reconstruct `g` at right of (i, j, k)

          Stencil1D sg;
          sg.m = g(i - 1, j, k);
          sg.c = g(i, j, k);
          sg.p = g(i + 1, j, k);
          cellboundary(sg); // Calculate sg.R and sg.L

          gedge = sg.R;
        }

        // Flux across cell edge
        const BoutReal fout = gradient * aedge * gedge
                              * (coord->J(i, j, k) * coord->g11(i, j, k)
                                 + coord->J(i + 1, j, k) * coord->g11(i + 1, j, k))
                              / (coord->dx(i, j, k) + coord->dx(i + 1, j, k));

        result(i, j, k) += fout / (coord->dx(i, j, k) * coord->J(i, j, k));
        result(i + 1, j, k) -= fout / (coord->dx(i + 1, j, k) * coord->J(i + 1, j, k));
      }
    }
  }

  const bool fci =
      f.hasParallelSlices() && a.hasParallelSlices() && g.hasParallelSlices();

#if BOUT_USE_METRIC_3D
  if (fci) {
    // 3D Metric, need yup/ydown fields.
    // Requires previous communication of metrics.
    if (!coord->g23.hasParallelSlices() || !coord->g_23.hasParallelSlices()
        || !coord->dy.hasParallelSlices() || !coord->dz.hasParallelSlices()
        || !coord->Bxy.hasParallelSlices() || !coord->J.hasParallelSlices()) {
      throw BoutException("metrics have no yup/down: Maybe communicate in init?");
    }
  }
#endif

  // Y and Z fluxes require Y derivatives

  // Values on this y slice (centre).
  // This is needed because toFieldAligned may modify the field
  const auto f_slice = makeslices(fci, f);
  const auto a_slice = makeslices(fci, a);
  const auto g_slice = makeslices(fci, g);

#if BOUT_USE_METRIC_3D
  const bool metric_fci = fci;
#else
  constexpr bool metric_fci = false;
#endif
  const auto g23 = makeslices(metric_fci, coord->g23);
  const auto g_23 = makeslices(metric_fci, coord->g_23);
  const auto J = makeslices(metric_fci, coord->J);
  const auto dy = makeslices(metric_fci, coord->dy);
  const auto dz = makeslices(metric_fci, coord->dz);
  const auto Bxy = makeslices(metric_fci, coord->Bxy);

  // Result of the Y and Z fluxes
  Field3D yzresult(0.0, mesh);
  if (!fci) {
    yzresult.setDirectionY(YDirectionType::Aligned);
  }

  // Y flux

  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
#if BOUT_USE_METRIC_3D
      for (int k = 0; k < mesh->LocalNz; k++)
#else
      const int k = 0;
#endif
      {
        const BoutReal coef_u =
            0.5
            * (g_23.c(i, j, k) / SQ(J.c(i, j, k) * Bxy.c(i, j, k))
               + g_23.up(i, j + 1, k) / SQ(J.up(i, j + 1, k) * Bxy.up(i, j + 1, k)));

        const BoutReal coef_d =
            0.5
            * (g_23.c(i, j, k) / SQ(J.c(i, j, k) * Bxy.c(i, j, k))
               + g_23.down(i, j - 1, k)
                     / SQ(J.down(i, j - 1, k) * Bxy.down(i, j - 1, k)));

#if not BOUT_USE_METRIC_3D
        for (int k = 0; k < mesh->LocalNz; k++)
#endif
        {
          // Calculate flux between j and j+1
          const int kp = (k + 1) % mesh->LocalNz;
          const int km = (k - 1 + mesh->LocalNz) % mesh->LocalNz;

          // Calculate Z derivative at y boundary
          BoutReal dfdz = 0.5
                          * (f_slice.c(i, j, kp) - f_slice.c(i, j, km)
                             + f_slice.up(i, j + 1, kp) - f_slice.up(i, j + 1, km))
                          / (dz.c(i, j, k) + dz.up(i, j + 1, k));

          // Y derivative
          BoutReal dfdy = 2. * (f_slice.up(i, j + 1, k) - f_slice.c(i, j, k))
                          / (dy.up(i, j + 1, k) + dy.c(i, j, k));

          BoutReal aedge = 0.5 * (a_slice.c(i, j, k) + a_slice.up(i, j + 1, k));
          BoutReal gedge = BoutNaN;
          if ((j == mesh->yend) and mesh->lastY(i)) {
            // Midpoint boundary value
            gedge = 0.5 * (g_slice.c(i, j, k) + g_slice.up(i, j + 1, k));
          } else if (dfdy > 0) {
            // Flux from (j+1) to (j)
            gedge = g_slice.up(i, j + 1, k);
          } else {
            // Flux from (j) to (j+1)
            gedge = g_slice.c(i, j, k);
          }

          BoutReal fout =
              aedge * gedge * 0.5
              * (J.c(i, j, k) * g23.c(i, j, k) + J.up(i, j + 1, k) * g23.up(i, j + 1, k))
              * (dfdz - coef_u * dfdy);

          yzresult(i, j, k) += fout / (dy.c(i, j, k) * J.c(i, j, k));

          // Calculate flux between j and j-1
          dfdz = 0.5
                 * (f_slice.c(i, j, kp) - f_slice.c(i, j, km) + f_slice.down(i, j - 1, kp)
                    - f_slice.down(i, j - 1, km))
                 / (dz.c(i, j, k) + dz.down(i, j - 1, k));

          dfdy = 2. * (f_slice.c(i, j, k) - f_slice.down(i, j - 1, k))
                 / (dy.c(i, j, k) + dy.down(i, j - 1, k));

          aedge = 0.5 * (a_slice.c(i, j, k) + a_slice.down(i, j - 1, k));
          if ((j == mesh->ystart) and mesh->firstY(i)) {
            gedge = 0.5 * (g_slice.c(i, j, k) + g_slice.down(i, j - 1, k));
          } else if (dfdy > 0) {
            gedge = g_slice.c(i, j, k);
          } else {
            gedge = g_slice.down(i, j - 1, k);
          }

          fout = aedge * gedge * 0.5
                 * (J.c(i, j, k) * g23.c(i, j, k)
                    + J.down(i, j - 1, k) * g23.down(i, j - 1, k))
                 * (dfdz - coef_d * dfdy);

          yzresult(i, j, k) -= fout / (dy.c(i, j, k) * J.c(i, j, k));
        }
      }
    }
  }

  // Z flux
  // Easier since all metrics constant in Z

  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
#if BOUT_USE_METRIC_3D
      for (int k = 0; k < mesh->LocalNz; k++)
#else
      const int k = 0;
#endif
      {
        // Coefficient in front of df/dy term
        const BoutReal coef =
            g_23.c(i, j, k)
            / (dy.up(i, j + 1, k) + 2. * dy.c(i, j, k) + dy.down(i, j - 1, k))
            / SQ(J.c(i, j, k) * Bxy.c(i, j, k));
#if not BOUT_USE_METRIC_3D
        for (int k = 0; k < mesh->LocalNz; k++)
#endif
        {
          // Calculate flux between k and k+1
          const int kp = (k + 1) % mesh->LocalNz;

          const BoutReal gradient =
              // df/dz
              ((f_slice.c(i, j, kp) - f_slice.c(i, j, k)) / dz.c(i, j, k))

              // - g_yz * df/dy / SQ(J*B)
              - (coef
                 * (f_slice.up(i, j + 1, k) + f_slice.up(i, j + 1, kp)
                    - f_slice.down(i, j - 1, k) - f_slice.down(i, j - 1, kp)));

          const BoutReal fout =
              gradient * 0.5 * (a_slice.c(i, j, kp) + a_slice.c(i, j, k))
              * ((gradient > 0) ? g_slice.c(i, j, kp) : g_slice.c(i, j, k));

          yzresult(i, j, k) += fout / dz.c(i, j, k);
          yzresult(i, j, kp) -= fout / dz.c(i, j, kp);
        }
      }
    }
  }
  // Check if we need to transform back
  if (fci) {
    result += yzresult;
  } else {
    result += fromFieldAligned(yzresult);
  }
  return result;
}

} // namespace FV
#endif // BOUT_FV_OPS_H
