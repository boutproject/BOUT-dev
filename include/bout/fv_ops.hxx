/*
  Finite-volume discretisation methods. Flux-conservative form
 */

#ifndef __FV_OPS_H__
#define __FV_OPS_H__

#include "../globals.hxx"
#include "../field3d.hxx"
#include "../vector2d.hxx"

#include "../utils.hxx"
#include <bout/mesh.hxx>

namespace FV {
  /*!
   * Div ( a Laplace_perp(x) )  -- Vorticity
   */
  const Field3D Div_a_Laplace_perp(const Field3D &a, const Field3D &x);
  
  /*!
   * Divergence of a parallel diffusion Div( k * Grad_par(f) )
   */
  const Field3D Div_par_K_Grad_par(const Field3D &k, const Field3D &f, bool bndry_flux=true);
  
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
   */
  const Field3D D4DY4(const Field3D &d, const Field3D &f);

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
  const Field3D D4DY4_Index(const Field3D &f, bool bndry_flux = true);
  
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
    void operator()(Stencil1D &n) {
      n.L = n.R = n.c;
    }
  };

  /*!
   * Fromm method
   */
  struct Fromm {
    void operator()(Stencil1D &n) {
      n.L = n.c - 0.25*(n.p - n.m);
      n.R = n.c + 0.25*(n.p - n.m);
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
    void operator()(Stencil1D &n) {
      // Choose the gradient within the cell
      // as the minimum (smoothest) solution
      BoutReal slope = _minmod(n.p - n.c, n.c - n.m);
      n.L = n.c - 0.5*slope;
      n.R = n.c + 0.5*slope;
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
      if( a*b <= 0.0 )
        return 0.0;
      
      if(fabs(a) < fabs(b))
        return a;
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
    void operator()(Stencil1D &n) {
      BoutReal slope = minmod(2. * (n.p - n.c),  // 2*right difference
                              0.5 * (n.p - n.m), // Central difference
                              2. * (n.c - n.m) ); // 2*left difference
      n.L = n.c - 0.5*slope;
      n.R = n.c + 0.5*slope;
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
   * Communicate fluxes between processors
   * Takes values in guard cells, and adds them to cells
   */
  void communicateFluxes(Field3D &f);

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
  /// @param[in] wave_speed  Maximum wave speed in the system
  /// @param[in] fixflux     Fix the flux at the boundary to be the value at the
  ///                        midpoint (for boundary conditions)
  ///
  /// NB: Uses to/from FieldAligned coordinates
  template<typename CellEdges = MC>
  const Field3D Div_par(const Field3D &f_in, const Field3D &v_in,
                        const Field3D &wave_speed, bool fixflux=true) {

    ASSERT1(areFieldsCompatible(f_in, v_in));
    ASSERT1(areFieldsCompatible(f_in, wave_speed));

    Mesh *localmesh = f_in.getMesh();
    
    CellEdges cellboundary;
    
    Field3D f = toFieldAligned(f_in, RGN_NOX);
    Field3D v = toFieldAligned(v_in, RGN_NOX);

    Coordinates *coord = f_in.getCoordinates();

    Field3D result{zeroFrom(f_in)};

    // Only need one guard cell, so no need to communicate fluxes
    // Instead calculate in guard cells to preserve fluxes
    int ys = localmesh->ystart-1;
    int ye = localmesh->yend+1;

    for (int i = localmesh->xstart; i <= localmesh->xend; i++) {

      if (!localmesh->firstY(i) || localmesh->periodicY(i)) {
        // Calculate in guard cell to get fluxes consistent between processors
        ys = localmesh->ystart - 1;
      } else {
        // Don't include the boundary cell. Note that this implies special
        // handling of boundaries later
        ys = localmesh->ystart;
      }

      if (!localmesh->lastY(i) || localmesh->periodicY(i)) {
        // Calculate in guard cells
        ye = localmesh->yend + 1;
      } else {
        // Not in boundary cells
        ye = localmesh->yend;
      }

      for (int j = ys; j <= ye; j++) {
        // Pre-calculate factors which multiply fluxes

        // For right cell boundaries
        BoutReal common_factor = (coord->J(i, j) + coord->J(i, j + 1)) /
          (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j + 1)));
        
        BoutReal flux_factor_rc = common_factor / (coord->dy(i, j) * coord->J(i, j));
        BoutReal flux_factor_rp = common_factor / (coord->dy(i, j + 1) * coord->J(i, j + 1));

        // For left cell boundaries
        common_factor = (coord->J(i, j) + coord->J(i, j - 1)) /
          (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j - 1)));

        BoutReal flux_factor_lc = common_factor / (coord->dy(i, j) * coord->J(i, j));
        BoutReal flux_factor_lm = common_factor / (coord->dy(i, j - 1) * coord->J(i, j - 1));
        
        for (int k = localmesh->zstart; k <= localmesh->zend; k++) {

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

          if (localmesh->lastY(i) && (j == localmesh->yend) && !localmesh->periodicY(i)) {
            // Last point in domain

            BoutReal bndryval = 0.5 * (s.c + s.p);
            if (fixflux) {
              // Use mid-point to be consistent with boundary conditions
              flux = bndryval * vpar;
            } else {
              // Add flux due to difference in boundary values
              flux = s.R * vpar + wave_speed(i, j, k) * (s.R - bndryval);
            }
          } else {
            
            // Maximum wave speed in the two cells
            BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k));

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

          if (localmesh->firstY(i) && (j == localmesh->ystart) && !localmesh->periodicY(i)) {
            // First point in domain
            BoutReal bndryval = 0.5 * (s.c + s.m);
            if (fixflux) {
              // Use mid-point to be consistent with boundary conditions
              flux = bndryval * vpar;
            } else {
              // Add flux due to difference in boundary values
              flux = s.L * vpar - wave_speed(i, j, k) * (s.L - bndryval);
            }
          } else {
            
            // Maximum wave speed in the two cells
            BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k));

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
    return fromFieldAligned(result, RGN_NOBNDRY);
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
  template<typename CellEdges = MC>
  const Field3D Div_f_v(const Field3D &n_in, const Vector3D &v, bool bndry_flux) {
    ASSERT1(n_in.getLocation() == v.getLocation());
    ASSERT1(areFieldsCompatible(n_in, v.x));

    Mesh* localmesh = n_in.getMesh();
    
    CellEdges cellboundary;
    
    Coordinates *coord = n_in.getCoordinates();

    if (v.covariant) {
      // Got a covariant vector instead
      throw BoutException("Div_f_v_XPPM passed a covariant v");
    }
    
    Field3D result{zeroFrom(n_in)};
    
    Field3D vx = v.x;
    Field3D vz = v.z;
    Field3D n  = n_in;
    
    // X-Z advection
    BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {
      // Calculate offset indices
      auto izp = i.zp();
      auto izm = i.zm();
      auto izpp = i.zpp();
      auto izmm = i.zmm();

      // Calculate velocities
      BoutReal vU = 0.5 * (vz[izp] + vz[i]) * coord->J[i];
      BoutReal vD = 0.5 * (vz[izm] + vz[i]) * coord->J[i];
      BoutReal vL = 0.25 * (vx[i.xm()] + vx[i]) * (coord->J[i.xm()] + coord->J[i]);
      BoutReal vR = 0.25 * (vx[i.xp()] + vx[i]) * (coord->J[i.xp()] + coord->J[i]);

      // X direction
      Stencil1D s;
      s.c  = n[i];
      s.m  = n[i.xm()];
      s.mm = n[i.xmm()];
      s.p  = n[i.xp()];
      s.pp = n[i.xpp()];
      
      cellboundary(s);

      if ((i.x() == localmesh->xend) && (localmesh->lastX())) {
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
          // Note: This statement might update another thread's data
          BOUT_OMP(atomic)
          result[i.xp()] -= flux / (coord->dx[i.xp()] * coord->J[i.xp()]);
        }
      } else {
        // Not at a boundary
        if (vR > 0.0) {
          // Flux out into next cell
          BoutReal flux = vR * s.R;
          result[i] += flux / (coord->dx[i] * coord->J[i]);
          // Note: This statement might update another thread's data
          BOUT_OMP(atomic)
          result[i.xp()] -= flux / (coord->dx[i.xp()] * coord->J[i.xp()]);
        }
      }
      
      // Left side

      if ((i.x() == localmesh->xstart) && (localmesh->firstX())) {
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
          BOUT_OMP(atomic)
          result[i.xm()] += flux / (coord->dx[i.xm()] * coord->J[i.xm()]);
        }
      } else {
        // Not at a boundary

        if (vL < 0.0) {
          BoutReal flux = vL * s.L;
          result[i] -= flux / (coord->dx[i] * coord->J[i]);
          BOUT_OMP(atomic)
          result[i.xm()] += flux / (coord->dx[i.xm()] * coord->J[i.xm()]);
        }
      }

      /// NOTE: Need to communicate fluxes
          
      // Z direction
      s.m  = n[izm];
      s.mm = n[izmm];
      s.p  = n[izp];
      s.pp = n[izpp];
      
      cellboundary(s);

      if (vU > 0.0) {
        BoutReal flux = vU * s.R / (coord->J[i] * coord->dz);
        result[i] += flux;
        BOUT_OMP(atomic)
        result[izp] -= flux;
      }
      if (vD < 0.0) {
        BoutReal flux = vD * s.L / (coord->J[i] * coord->dz);
        result[i] -= flux;
        BOUT_OMP(atomic)
        result[izm] += flux;
      }
    }
    communicateFluxes(result);
    
    // Y advection
    // Currently just using simple centered differences
    // so no fluxes need to be exchanged
    
    n = toFieldAligned(n_in, RGN_NOX);
    Field3D vy = toFieldAligned(v.y, RGN_NOX);

    Field3D yresult{0.0, localmesh};

    BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {
      // Y velocities on y boundaries
      BoutReal vU = 0.25 * (vy[i] + vy[i.yp()]) * (coord->J[i] + coord->J[i.yp()]);
      BoutReal vD = 0.25 * (vy[i] + vy[i.ym()]) * (coord->J[i] + coord->J[i.ym()]);

      // n (advected quantity) on y boundaries
      // Note: Use unshifted n_in variable
      BoutReal nU = 0.5 * (n[i] + n[i.yp()]);
      BoutReal nD = 0.5 * (n[i] + n[i.ym()]);

      yresult[i] = (nU * vU - nD * vD) / (coord->J[i] * coord->dy[i]);
    }

    return result + fromFieldAligned(yresult, RGN_NOBNDRY);
  }
}

#endif // __FV_OPS_H__
