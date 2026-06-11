/*
  Finite-volume discretisation methods. Flux-conservative form
 */

#ifndef BOUT_FV_OPS_H
#define BOUT_FV_OPS_H

#include "boutexception.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/build_defines.hxx"
#include "bout/coordinates.hxx"
#include "bout/field.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"
#include "bout/vector2d.hxx"
#include <bout/mesh.hxx>

namespace FV {
/*!
 * Div ( a Grad_perp(f) ) -- ∇⊥ ( a ⋅ ∇⊥ f) -- Vorticity
 */
Field3D Div_a_Grad_perp(const Field3D& a, const Field3D& x);

[[deprecated("Please use Div_a_Grad_perp instead")]] inline Field3D
Div_a_Laplace_perp(const Field3D& a, const Field3D& x) {
  return Div_a_Grad_perp(a, x);
}

/*!
   * Divergence of a parallel diffusion Div( k * Grad_par(f) )
   */
const Field3D Div_par_K_Grad_par(const Field3D& k, const Field3D& f,
                                 bool bndry_flux = true);

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
const Field3D D4DY4(const Field3D& d, const Field3D& f);

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
const Field3D D4DY4_Index(const Field3D& f, bool bndry_flux = true);

/*!
   * Stencil used for Finite Volume calculations
   * which includes cell face values L and R
   */
struct Stencil1D {
  // Cell centre values
  BoutReal c, m, p, mm, pp;

  // Left and right cell face values
  BoutReal L, R;
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
  BoutReal _minmod(BoutReal a, BoutReal b) {
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
    n.L = n.c - 0.5 * slope;
    n.R = n.c + 0.5 * slope;
  }

private:
  // Return zero if any signs are different
  // otherwise return the value with the minimum magnitude
  BoutReal minmod(BoutReal a, BoutReal b, BoutReal c) {
    // if any of the signs are different, return zero gradient
    if ((a * b <= 0.0) || (a * c <= 0.0)) {
      return 0.0;
    }

    // Return the minimum absolute value
    return SIGN(a) * BOUTMIN(fabs(a), fabs(b), fabs(c));
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
struct VanAlbada {
  void operator()(Stencil1D& n) {
    const BoutReal dl = n.c - n.m;
    const BoutReal dr = n.p - n.c;

    const BoutReal denom = dl * dl + dr * dr;

    // Smoothness parameters:
    // - keep division well-defined when dl=dr=0
    // - provide a differentiable approximation to max(dl*dr, 0)
    const BoutReal eps = 1e-12 * denom + 1e-30;

    const BoutReal ab = dl * dr;
    const BoutReal ab_pos = 0.5 * (ab + sqrt(ab * ab + eps * eps));

    const BoutReal slope = (ab_pos * (dl + dr)) / (denom + eps);

    n.L = n.c - 0.5 * slope;
    n.R = n.c + 0.5 * slope;
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
    const BoutReal eps = 1e-12 * (beta0_r + beta1_r) + 1e-30;

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

    n.R = w0_r * p0_r + w1_r * p1_r;
    n.L = w0_l * p0_l + w1_l * p1_l;
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

  Field3D f = are_unaligned ? toFieldAligned(f_in, "RGN_NOX") : f_in;
  Field3D v = are_unaligned ? toFieldAligned(v_in, "RGN_NOX") : v_in;
  Field3D wave_speed =
      are_unaligned ? toFieldAligned(wave_speed_in, "RGN_NOX") : wave_speed_in;

  Coordinates* coord = f_in.getCoordinates();

  Field3D result{zeroFrom(f)};

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
      for (int k = mesh->zstart; k <= mesh->zend; k++) {
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
        BoutReal flux;

        if (mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
          // Last point in domain

          const BoutReal bndryval = 0.5 * (s.c + s.p);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux = s.R * vpar + wave_speed(i, j, k) * (s.R - bndryval);
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

        if (mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
          // First point in domain
          const BoutReal bndryval = 0.5 * (s.c + s.m);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux = s.L * vpar - wave_speed(i, j, k) * (s.L - bndryval);
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

  Mesh* mesh = n_in.getMesh();

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
        BoutReal flux;
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
        BoutReal flux;
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
} // namespace FV
#endif // BOUT_FV_OPS_H
