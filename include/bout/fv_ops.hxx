/*
  Finite-volume discretisation methods. Flux-conservative form
 */

#ifndef __FV_OPS_H__
#define __FV_OPS_H__

#include "../globals.hxx"
#include "../field3d.hxx"
#include "../vector2d.hxx"

#include "../utils.hxx"

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

    ASSERT2(f_in.getLocation() == v_in.getLocation());

    CellEdges cellboundary;
    
    Field3D f = mesh->toFieldAligned(f_in);
    Field3D v = mesh->toFieldAligned(v_in);

    Coordinates *coord = f_in.getCoordinates();

    Field3D result = 0.0;
    
    // Only need one guard cell, so no need to communicate fluxes
    // Instead calculate in guard cells to preserve fluxes
    int ys = mesh->ystart-1;
    int ye = mesh->yend+1;

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
        
        for (int k = 0; k < mesh->LocalNz; k++) {

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

          if (mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
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
    return mesh->fromFieldAligned(result);
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
    ASSERT2(n_in.getLocation() == v.getLocation());

    CellEdges cellboundary;
    
    Coordinates *coord = n_in.getCoordinates();
    
    if(v.covariant) {
      // Got a covariant vector instead
      throw BoutException("Div_f_v_XPPM passed a covariant v");
    }
    
    Field3D result = 0.0;
    
    Field3D vx = v.x;
    Field3D vz = v.z;
    Field3D n  = n_in;
    
    // X-Z advection
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++)
        for(int k=0;k<mesh->LocalNz;k++) {
          int kp = (k+1) % (mesh->LocalNz);
          int kpp = (kp+1) % (mesh->LocalNz);
          int km = (k-1+mesh->LocalNz) % (mesh->LocalNz);
          int kmm = (km-1+mesh->LocalNz) % (mesh->LocalNz);
          
          // Calculate velocities
          BoutReal vU = 0.5*(vz(i,j,kp) + vz(i,j,k))*coord->J(i,j);
          BoutReal vD = 0.5*(vz(i,j,km) + vz(i,j,k))*coord->J(i,j);
          BoutReal vL = 0.25*(vx(i-1,j,k) + vx(i,j,k))*(coord->J(i-1,j) + coord->J(i,j));
          BoutReal vR = 0.25*(vx(i+1,j,k) + vx(i,j,k))*(coord->J(i+1,j) + coord->J(i,j));
          
          // X direction
          Stencil1D s;
          s.c  = n(i,  j,k);
          s.m  = n(i-1,j,k);
          s.mm = n(i-2,j,k);
          s.p  = n(i+1,j,k);
          s.pp = n(i+2,j,k);
          
          cellboundary(s);
        
          if((i==mesh->xend) && (mesh->lastX())) {
            // At right boundary in X
          
            if(bndry_flux) {
              BoutReal flux;
              if(vR > 0.0) {
                // Flux to boundary
                flux = vR * s.R;
              }else {
                // Flux in from boundary
                flux = vR * 0.5*(n(i+1,j,k) + n(i,j,k));
              }
              result(i,j,k)   += flux / (coord->dx(i,j) * coord->J(i,j));
              result(i+1,j,k) -= flux / (coord->dx(i+1,j) * coord->J(i+1,j));
            }
          }else {
            // Not at a boundary
            if(vR > 0.0) {
              // Flux out into next cell
              BoutReal flux = vR * s.R;
              result(i,j,k)   += flux / (coord->dx(i,j) * coord->J(i,j));
              result(i+1,j,k) -= flux / (coord->dx(i+1,j) * coord->J(i+1,j));
            }
          }
        
          // Left side
          
          if((i==mesh->xstart) && (mesh->firstX())) {
            // At left boundary in X
            
            if(bndry_flux) {
              BoutReal flux;
              
              if(vL < 0.0) {
                // Flux to boundary
                flux = vL * s.L;
                
              }else {
                // Flux in from boundary
                flux = vL * 0.5*(n(i-1,j,k) + n(i,j,k));
              }
              result(i,j,k)   -= flux / (coord->dx(i,j) * coord->J(i,j));
              result(i-1,j,k) += flux / (coord->dx(i-1,j) * coord->J(i-1,j));
            }
          }else {
            // Not at a boundary
            
            if(vL < 0.0) {
              BoutReal flux = vL * s.L;
              result(i,j,k)   -= flux / (coord->dx(i,j) * coord->J(i,j));
              result(i-1,j,k) += flux / (coord->dx(i-1,j) * coord->J(i-1,j));
            }
          }
          
          /// NOTE: Need to communicate fluxes
          
          // Z direction
          s.m  = n(i,j,km);
          s.mm = n(i,j,kmm);
          s.p  = n(i,j,kp);
          s.pp = n(i,j,kpp);
          
          cellboundary(s);
          
          if(vU > 0.0) {
            BoutReal flux = vU * s.R / (coord->J(i,j)*coord->dz);
            result(i,j,k)   += flux;
            result(i,j,kp)  -= flux;
          }
          if(vD < 0.0) {
            BoutReal flux = vD * s.L / (coord->J(i,j)*coord->dz);
            result(i,j,k)   -= flux;
            result(i,j,km)  += flux;
          }
          
        }
    communicateFluxes(result);
    
    // Y advection
    // Currently just using simple centered differences
    // so no fluxes need to be exchanged
    
    n = mesh->toFieldAligned(n_in);
    Field3D vy = mesh->toFieldAligned(v.y);
    
    Field3D yresult = 0.0;    
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++)
        for(int k=0;k<mesh->LocalNz;k++) {
          
          // Y velocities on y boundaries
          BoutReal vU = 0.25*(vy(i,j,k) + vy(i,j+1,k))*(coord->J(i,j) + coord->J(i,j+1));
          BoutReal vD = 0.25*(vy(i,j,k) + vy(i,j-1,k))*(coord->J(i,j) + coord->J(i,j-1));
          
          // n (advected quantity) on y boundaries
          // Note: Use unshifted n_in variable
          BoutReal nU = 0.5*(n(i,j,k) + n(i,j+1,k));
          BoutReal nD = 0.5*(n(i,j,k) + n(i,j-1,k));
          
          yresult(i,j,k) = (nU*vU - nD*vD) / (coord->J(i,j)*coord->dy(i,j));
        }
    
    return result + mesh->fromFieldAligned(yresult);
  }
}

#endif // __FV_OPS_H__
