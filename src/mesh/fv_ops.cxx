
#include <bout/fv_ops.hxx>
#include <globals.hxx>
#include <utils.hxx>
#include <msg_stack.hxx>

#include <output.hxx>

namespace FV {

  // Div ( a Laplace_perp(f) )  -- Vorticity
  const Field3D Div_a_Laplace_perp(const Field3D &a, const Field3D &f) {
    ASSERT1(a.getMesh() == f.getMesh());
    ASSERT2(a.getLocation() == f.getLocation());

    Mesh *localmesh = a.getMesh();

    Field3D result{0.0, localmesh};

    Coordinates *coord = f.getCoordinates();
    
    // Flux in x
  
    int xs = localmesh->xstart-1;
    int xe = localmesh->xend;
    
    for (int i = xs; i <= xe; i++)
      for (int j = localmesh->ystart; j <= localmesh->yend; j++) {
        for (int k = localmesh->zstart; k <= localmesh->zend; k++) {
          
          // Calculate flux from i to i+1

          BoutReal fout = (coord->J(i, j) * a(i, j, k) * coord->g11(i, j)
                           + coord->J(i + 1, j) * a(i + 1, j, k) * coord->g11(i + 1, j))
                          * (f(i + 1, j, k) - f(i, j, k))
                          / (coord->dx(i, j) + coord->dx(i + 1, j));

          result(i, j, k) += fout / (coord->dx(i, j) * coord->J(i, j));
          result(i + 1, j, k) -= fout / (coord->dx(i + 1, j) * coord->J(i + 1, j));
        }
      }

    // Y flux
    BOUT_FOR(i, f.getRegion("RGN_NOBNDRY")) {
    
      BoutReal coef = 0.5*(coord->g_23[i]/SQ(coord->J[i]*coord->Bxy[i]) + coord->g_23[i]/SQ(coord->J[i]*coord->Bxy[i]));
        
      // Calculate flux between j and j+1
      
      // Calculate Z derivative at y boundary
      BoutReal dfdz = 0.25
                      * (f[i.zp()] - f[i.zm()] + f.yup()[i.offset(0, 1, 1)]
                         - f.yup()[i.offset(0, 1, -1)])
                      / coord->dz;

      // Y derivative
      BoutReal dfdy = 2. * (f.yup()[i.yp()] - f[i]) / (coord->dy[i.yp()] + coord->dy[i]);

      BoutReal fout = 0.5
                      * (coord->J[i] * a[i] * coord->g23[i]
                         + coord->J[i.yp()] * a.yup()[i.yp()] * coord->g23[i.yp()])
                      * (dfdz - coef * dfdy);

      result[i] += fout / (coord->dy[i] * coord->J[i]);

      // Calculate flux between j and j-1
      dfdz = 0.25
             * (f[i.zp()] - f[i.zm()] + f.ydown()[i.offset(0, -1, 1)]
                - f.ydown()[i.offset(0, -1, -1)])
             / coord->dz;

      dfdy = 2. * (f[i] - f.ydown()[i.ym()]) / (coord->dy[i] + coord->dy[i.ym()]);

      fout = 0.5
             * (coord->J[i] * a[i] * coord->g23[i]
                + coord->J[i.ym()] * a.ydown()[i.ym()] * coord->g23[i.ym()])
             * (dfdz - coef * dfdy);

      result[i] -= fout / (coord->dy[i] * coord->J[i]);
    }

    // Z flux
    // Easier since all metrics constant in Z
    BOUT_FOR(i, f.getRegion("RGN_NOBNDRY")) {
      
	// Coefficient in front of df/dy term
        BoutReal coef = coord->g_23[i]
                        / (coord->dy[i.yp()] + 2. * coord->dy[i] + coord->dy[i.ym()])
                        / SQ(coord->J[i] * coord->Bxy[i]);

        // Calculate flux between k and k+1

        BoutReal fout =
            0.5 * (a[i] + a[i.xp()]) * coord->g33[i]
            * (
                  // df/dz
                  (f[i.zp()] - f[i]) / coord->dz

                  // - g_yz * df/dy / SQ(J*B)
                  - coef
                        * (f.yup()[i.yp()] + f.yup()[i.offset(0, 1, 1)]
                           - f.ydown()[i.ym()] - f.ydown()[i.offset(0, -1, 1)]));

        result[i] += fout / coord->dz;
        // Note: This statement might update another thread's data
        BOUT_OMP(atomic)
        result[i.zp()] -= fout / coord->dz;
    }
  
    return result;
  }

  const Field3D Div_par_K_Grad_par(const Field3D &Kin, const Field3D &fin, bool bndry_flux) {
    TRACE("FV::Div_par_K_Grad_par");

    ASSERT1(Kin.getMesh() == fin.getMesh());
    ASSERT2(Kin.getLocation() == fin.getLocation());

    Mesh *localmesh = Kin.getMesh();
    Field3D result{0.0, localmesh};

    bool use_yup_ydown = (Kin.hasYupYdown() && fin.hasYupYdown());

    const auto& K = use_yup_ydown ? Kin : localmesh->toFieldAligned(Kin);
    const auto& f = use_yup_ydown ? fin : localmesh->toFieldAligned(fin);

    // K and f fields in yup and ydown directions
    const auto& Kup = use_yup_ydown ? Kin.yup() : K;
    const auto& Kdown = use_yup_ydown ? Kin.ydown() : K;
    const auto& fup = use_yup_ydown ? fin.yup() : f;
    const auto& fdown = use_yup_ydown ? fin.ydown() : f;
    
    Coordinates *coord = fin.getCoordinates();

    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      // Calculate flux at upper surface
      
      const auto iyp = i.yp();
      const auto iym = i.ym();

      if (bndry_flux || !mesh->lastY() || (i.y() != mesh->yend)) {

        BoutReal c = 0.5*(K[i] + Kup[iyp]); // K at the upper boundary
        BoutReal J = 0.5*(coord->J[i] + coord->J[iyp]); // Jacobian at boundary
        BoutReal g_22 = 0.5*(coord->g_22[i] + coord->g_22[iyp]);
        
        BoutReal gradient = 2.*(fup[iyp] - f[i]) / (coord->dy[i] + coord->dy[iyp]);
        
        BoutReal flux = c * J * gradient / g_22;
            
        result[i] += flux / (coord->dy[i] * coord->J[i]);
      }
      
      // Calculate flux at lower surface
      if (bndry_flux || !mesh->firstY() || (i.y() != mesh->ystart)) {
        BoutReal c = 0.5*(K[i] + Kdown[iym]); // K at the lower boundary
        BoutReal J = 0.5*(coord->J[i] + coord->J[iym]); // Jacobian at boundary
        
        BoutReal g_22 = 0.5*(coord->g_22[i] + coord->g_22[i]);
        
        BoutReal gradient = 2.*(f[i] - fdown[iym]) / (coord->dy[i] + coord->dy[iym]);
        
        BoutReal flux = c * J * gradient / g_22;
        
        result[i] -= flux / (coord->dy[i] * coord->J[i]);
      }
    }
    
    if (!use_yup_ydown) {
      // Shifted to field aligned coordinates, so need to shift back
      result = localmesh->fromFieldAligned(result);
    }
    
    return result;
  }

  const Field3D D4DY4(const Field3D &d_in, const Field3D &f_in) {
    ASSERT1(d_in.getMesh() == f_in.getMesh());
    ASSERT2(d_in.getLocation() == f_in.getLocation());

    Mesh* localmesh = f_in.getMesh();
    
    Field3D result{0.0, localmesh};
    result.setLocation(f_in.getLocation());
    
    Coordinates *coord = f_in.getCoordinates();
    
    // Convert to field aligned coordinates
    Field3D d = localmesh->toFieldAligned(d_in);
    Field3D f = localmesh->toFieldAligned(f_in);

    for (int i = localmesh->xstart; i <= localmesh->xend; i++)
      for (int j = localmesh->ystart; j <= localmesh->yend; j++) {
        BoutReal dy3 = SQ(coord->dy(i, j)) * coord->dy(i, j);
        for (int k = localmesh->zstart; k <= localmesh->zend; k++) {

          // 3rd derivative at right boundary
          
          BoutReal d3fdx3 = (
                             f(i,j+2,k)
                             - 3.*f(i,j+1,k)
                             + 3.*f(i,j,  k)
                             -    f(i,j-1,k)
                             ) / dy3;
          
          BoutReal flux = 0.5*(d(i,j,k) + d(i,j+1,k))*(coord->J(i,j) + coord->J(i,j+1)) * d3fdx3;
          
          result(i,j,  k) += flux / (coord->J(i,j) * coord->dy(i,j));
          result(i,j+1,k) -= flux / (coord->J(i,j+1) * coord->dy(i,j+1));
          
          if(j == localmesh->ystart && (!localmesh->firstY())) {
            // Left cell boundary, no flux through boundaries
            d3fdx3 = (
                      f(i,j+1,k)
                      - 3.*f(i,j,  k)
                      + 3.*f(i,j-1,k)
                      -    f(i,j-2,k)
                      ) / dy3;
            
            flux = 0.5*(d(i,j,k) + d(i,j-1,k))*(coord->J(i,j) + coord->J(i,j-1)) * d3fdx3;
            
            result(i,j,  k) -= flux / (coord->J(i,j) * coord->dy(i,j));
            result(i,j-1,k) += flux / (coord->J(i,j-1) * coord->dy(i,j-1));
          }
        }
      }

    // Convert result back to non-aligned coordinates
    return localmesh->fromFieldAligned(result);
  }

  const Field3D D4DY4_Index(const Field3D &f_in, bool bndry_flux) {
    Mesh* localmesh = f_in.getMesh();

    Field3D result{0.0, localmesh};
    result.setLocation(f_in.getLocation());
    
    // Convert to field aligned coordinates
    Field3D f = localmesh->toFieldAligned(f_in);

    Coordinates *coord = f_in.getCoordinates();

    for (int i = localmesh->xstart; i <= localmesh->xend; i++) {
      bool yperiodic = localmesh->periodicY(i);
      
      bool has_upper_boundary = !yperiodic && localmesh->lastY(i);
      bool has_lower_boundary = !yperiodic && localmesh->firstY(i);

      for (int j = localmesh->ystart; j <= localmesh->yend; j++) {

        // Right boundary
        
        if (bndry_flux || (j != localmesh->yend) || !has_upper_boundary) {
          // Calculate the fluxes
          
          // Right boundary common factors
          BoutReal common_factor = 0.25 *
            (coord->dy(i, j) + coord->dy(i, j + 1)) *
            (coord->J(i, j) + coord->J(i, j + 1));
          
          BoutReal factor_rc = common_factor / (coord->J(i,j) * coord->dy(i,j));
          BoutReal factor_rp = common_factor / (coord->J(i,j+1) * coord->dy(i,j+1));
          if ( j != localmesh->yend || !has_upper_boundary ) {
            for (int k = localmesh->zstart; k <= localmesh->zend; k++) {
              // Not on domain boundary
              // 3rd derivative at right cell boundary
              
              BoutReal d3fdx3 = (
                                 f(i,j+2,k)
                                 - 3.*f(i,j+1,k)
                                 + 3.*f(i,j,  k)
                                 -    f(i,j-1,k)
                                 );
              
              result(i,j,  k) += d3fdx3 * factor_rc; 
              result(i,j+1,k) -= d3fdx3 * factor_rp;
            }
          } else {
            // At a domain boundary
            // Use a one-sided difference formula

            for (int k = localmesh->zstart; k <= localmesh->zend; k++) {

              BoutReal d3fdx3 = -((16. / 5) * 0.5 *
                                  (f(i, j + 1, k) + f(i, j, k)) // Boundary value f_b
                                  - 6. * f(i, j, k)                 // f_0
                                  + 4. * f(i, j - 1, k)             // f_1
                                  - (6. / 5) * f(i, j - 2, k)       // f_2
                                  );
              
              result(i,j,  k) += d3fdx3 * factor_rc; 
              result(i,j+1,k) -= d3fdx3 * factor_rp;
            }
          }
        }

        // Left cell boundary
        
        if (bndry_flux || (j != localmesh->ystart) || !has_lower_boundary) {
          // Calculate the fluxes

          BoutReal common_factor = 0.25 *
            (coord->dy(i, j) + coord->dy(i, j + 1)) *
            (coord->J(i, j) + coord->J(i, j - 1));

          BoutReal factor_lc = common_factor / (coord->J(i, j) * coord->dy(i, j));
          BoutReal factor_lm = common_factor / (coord->J(i, j - 1) * coord->dy(i, j - 1));
            
          if ( j != localmesh->ystart || !has_lower_boundary ) {
            for (int k = localmesh->zstart; k <= localmesh->zend; k++) {

              // Not on a domain boundary
              BoutReal d3fdx3 = (f(i, j + 1, k)
                                 - 3. * f(i, j, k)
                                 + 3. * f(i, j - 1, k)
                                 - f(i, j - 2, k));
              
              result(i, j    , k) -= d3fdx3 * factor_lc; 
              result(i, j - 1, k) += d3fdx3 * factor_lm;
            }
          } else {
            // On a domain (Y) boundary
            for (int k = localmesh->zstart; k <= localmesh->zend; k++) {
              BoutReal d3fdx3 =
                  -(-(16. / 5) * 0.5 * (f(i, j - 1, k) + f(i, j, k)) // Boundary value f_b
                    + 6. * f(i, j, k)                                // f_0
                    - 4. * f(i, j + 1, k)                            // f_1
                    + (6. / 5) * f(i, j + 2, k)                      // f_2
                  );

              result(i, j    , k) -= d3fdx3 * factor_lc; 
              result(i, j - 1, k) += d3fdx3 * factor_lm;
            }
          }
        }
      }
    }

    // Convert result back to non-aligned coordinates
    return localmesh->fromFieldAligned(result);
  }

  void communicateFluxes(Field3D &f) {
    Mesh* localmesh = f.getMesh();
    
    // Use X=0 as temporary buffer
    if (localmesh->xstart != 2)
      throw BoutException("communicateFluxes: Sorry!");

    int size = localmesh->LocalNy * localmesh->LocalNz;
    comm_handle xin, xout;
    // Cache results to silence spurious compiler warning about xin,
    // xout possibly being uninitialised when used
    bool not_first = !localmesh->firstX();
    bool not_last = !localmesh->lastX();
    if (not_first) {
      xin = localmesh->irecvXIn(f(0, 0), size, 0);
    }
    if (not_last) {
      xout = localmesh->irecvXOut(f(localmesh->LocalNx - 1, 0), size, 1);
    }
    // Send X=1 values
    if (not_first) {
      localmesh->sendXIn(f(1, 0), size, 1);
    }
    if (not_last) {
      localmesh->sendXOut(f(localmesh->LocalNx - 2, 0), size, 0);
    }
    // Wait
    if (not_first) {
      localmesh->wait(xin);
      // Add to cells
      for (int y = localmesh->ystart; y <= localmesh->yend; y++)
        for (int z = 0; z < localmesh->LocalNz; z++) {
          f(2, y, z) += f(0, y, z);
        }
    }
    if (not_last) {
      mesh->wait(xout);
      // Add to cells
      for (int y = localmesh->ystart; y <= localmesh->yend; y++)
        for (int z = 0; z < localmesh->LocalNz; z++) {
          f(localmesh->LocalNx - 3, y, z) += f(localmesh->LocalNx - 1, y, z);
        }
    }
  }

} // Namespace FV
