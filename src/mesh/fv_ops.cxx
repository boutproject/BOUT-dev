
#include <bout/fv_ops.hxx>
#include <globals.hxx>
#include <utils.hxx>
#include <msg_stack.hxx>

#include <output.hxx>

namespace FV {

  // Div ( a Laplace_perp(f) )  -- Vorticity
  const Field3D Div_a_Laplace_perp(const Field3D &a, const Field3D &f) {
    ASSERT2(a.getLocation() == f.getLocation());

    Mesh *mesh = a.getMesh();

    Field3D result{zeroFrom(f)};

    Coordinates *coord = f.getCoordinates();
    
    // Flux in x
  
    int xs = mesh->xstart-1;
    int xe = mesh->xend;
  
    /*
      if(mesh->firstX())
      xs += 1;
    */
    /*
      if(mesh->lastX())
      xe -= 1;
    */

    for(int i=xs;i<=xe;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
	for(int k=0;k<mesh->LocalNz;k++) {
	  // Calculate flux from i to i+1

          BoutReal fout = 0.5 * (a(i, j, k) + a(i + 1, j, k))
                          * (coord->J(i, j, k) * coord->g11(i, j, k)
                             + coord->J(i + 1, j, k) * coord->g11(i + 1, j, k))
                          * (f(i + 1, j, k) - f(i, j, k))
                          / (coord->dx(i, j, k) + coord->dx(i + 1, j, k));

          result(i, j, k) += fout / (coord->dx(i, j, k) * coord->J(i, j, k));
          result(i + 1, j, k) -= fout / (coord->dx(i + 1, j, k) * coord->J(i + 1, j, k));
        }
      }


    // Y and Z fluxes require Y derivatives

    // Fields containing values along the magnetic field
    Field3D fup(mesh), fdown(mesh);
    Field3D aup(mesh), adown(mesh);

    // Values on this y slice (centre).
    // This is needed because toFieldAligned may modify the field
    Field3D fc = f;
    Field3D ac = a;

    // Result of the Y and Z fluxes
    Field3D yzresult(mesh);
    yzresult.allocate();

    if (f.hasParallelSlices() && a.hasParallelSlices()) {
      // Both inputs have yup and ydown

      fup = f.yup();
      fdown = f.ydown();

      aup = a.yup();
      adown = a.ydown();
    } else {
      // At least one input doesn't have yup/ydown fields.
      // Need to shift to/from field aligned coordinates

      fup = fdown = fc = toFieldAligned(f);
      aup = adown = ac = toFieldAligned(a);
      yzresult.setDirectionY(YDirectionType::Aligned);
    }

    // Y flux

    for (int i = mesh->xstart; i <= mesh->xend; i++) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        for (int k = 0; k < mesh->LocalNz; k++) {
          // Calculate flux between j and j+1
          int kp = (k + 1) % mesh->LocalNz;
          int km = (k - 1 + mesh->LocalNz) % mesh->LocalNz;

          BoutReal coef =
              0.5
              * (coord->g_23(i, j, k) / SQ(coord->J(i, j, k) * coord->Bxy(i, j, k))
                 + coord->g_23(i, j + 1, k)
                       / SQ(coord->J(i, j + 1, k) * coord->Bxy(i, j + 1, k)));

          // Calculate Z derivative at y boundary
          BoutReal dfdz =
              0.25 * (fc(i, j, kp) - fc(i, j, km) + fup(i, j + 1, kp) - fup(i, j + 1, km))
              / coord->dz(i, j, k); // is this at right location?

          // Y derivative
          BoutReal dfdy = 2. * (fup(i, j + 1, k) - fc(i, j, k))
                          / (coord->dy(i, j + 1, k) + coord->dy(i, j, k));

          BoutReal fout = 0.25 * (ac(i, j, k) + aup(i, j + 1, k))
                          * (coord->J(i, j, k) * coord->g23(i, j, k)
                             + coord->J(i, j + 1, k) * coord->g23(i, j + 1, k))
                          * (dfdz - coef * dfdy);

          yzresult(i, j, k) = fout / (coord->dy(i, j, k) * coord->J(i, j, k));

          // Calculate flux between j and j-1
          dfdz =
              0.25
              * (fc(i, j, kp) - fc(i, j, km) + fdown(i, j - 1, kp) - fdown(i, j - 1, km))
              / coord->dz(i, j, k); // is this at right location?

          dfdy = 2. * (fc(i, j, k) - fdown(i, j - 1, k))
                 / (coord->dy(i, j, k) + coord->dy(i, j - 1, k));

          fout = 0.25 * (ac(i, j, k) + adown(i, j - 1, k))
                 * (coord->J(i, j, k) * coord->g23(i, j, k)
                    + coord->J(i, j - 1, k) * coord->g23(i, j - 1, k))
                 * (dfdz - coef * dfdy);

          yzresult(i, j, k) -= fout / (coord->dy(i, j, k) * coord->J(i, j, k));
        }
      }
    }

    // Z flux
    // Easier since all metrics constant in Z

    for (int i = mesh->xstart; i <= mesh->xend; i++) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        for (int k = 0; k < mesh->LocalNz; k++) {
          // Calculate flux between k and k+1
          int kp = (k + 1) % mesh->LocalNz;

          // Coefficient in front of df/dy term
          BoutReal coef = coord->g_23(i, j, k)
                          / (coord->dy(i, j + 1, k) + 2. * coord->dy(i, j, k)
                             + coord->dy(i, j - 1, k))
                          / SQ(coord->J(i, j, k) * coord->Bxy(i, j, k));

          BoutReal fout = 0.5 * (ac(i, j, k) + ac(i, j, kp)) * coord->g33(i, j, k)
                          * ( // df/dz
                                (fc(i, j, kp) - fc(i, j, k)) / coord->dz(i, j, k)
                                // - g_yz * df/dy / SQ(J*B)
                                - coef
                                      * (fup(i, j + 1, k) + fup(i, j + 1, kp)
                                         - fdown(i, j - 1, k) - fdown(i, j - 1, kp)));

          yzresult(i, j, k) += fout / coord->dz(i, j, k);
          yzresult(i, j, kp) -= fout / coord->dz(i, j, kp);
        }
      }
    }
    // Check if we need to transform back
    if (f.hasParallelSlices() && a.hasParallelSlices()) {
      result += yzresult;
    } else {
      result += fromFieldAligned(yzresult);
    }
    
    return result;
  }

  const Field3D Div_par_K_Grad_par(const Field3D &Kin, const Field3D &fin, bool bndry_flux) {
    TRACE("FV::Div_par_K_Grad_par");

    ASSERT2(Kin.getLocation() == fin.getLocation());

    Mesh *mesh = Kin.getMesh();

    bool use_parallel_slices = (Kin.hasParallelSlices() && fin.hasParallelSlices());

    const auto& K = use_parallel_slices ? Kin : toFieldAligned(Kin, "RGN_NOX");
    const auto& f = use_parallel_slices ? fin : toFieldAligned(fin, "RGN_NOX");

    Field3D result{zeroFrom(f)};

    // K and f fields in yup and ydown directions
    const auto& Kup = use_parallel_slices ? Kin.yup() : K;
    const auto& Kdown = use_parallel_slices ? Kin.ydown() : K;
    const auto& fup = use_parallel_slices ? fin.yup() : f;
    const auto& fdown = use_parallel_slices ? fin.ydown() : f;
    
    Coordinates *coord = fin.getCoordinates();

    BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
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
        
        BoutReal g_22 = 0.5*(coord->g_22[i] + coord->g_22[iym]);
        
        BoutReal gradient = 2.*(f[i] - fdown[iym]) / (coord->dy[i] + coord->dy[iym]);
        
        BoutReal flux = c * J * gradient / g_22;
        
        result[i] -= flux / (coord->dy[i] * coord->J[i]);
      }
    }
    
    if (!use_parallel_slices) {
      // Shifted to field aligned coordinates, so need to shift back
      result = fromFieldAligned(result, "RGN_NOBNDRY");
    }
    
    return result;
  }

  const Field3D D4DY4(const Field3D &d_in, const Field3D &f_in) {
    ASSERT1_FIELDS_COMPATIBLE(d_in, f_in);

    Mesh* mesh = d_in.getMesh();

    Coordinates *coord = f_in.getCoordinates();
    
    ASSERT2(d_in.getDirectionY() == f_in.getDirectionY());
    const bool are_unaligned = ((d_in.getDirectionY() == YDirectionType::Standard)
                                and (f_in.getDirectionY() == YDirectionType::Standard));

    // Convert to field aligned coordinates
    Field3D d = are_unaligned ? toFieldAligned(d_in, "RGN_NOX") : d_in;
    Field3D f = are_unaligned ? toFieldAligned(f_in, "RGN_NOX") : f_in;

    Field3D result{zeroFrom(f)};
    
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
        for(int k=0;k<mesh->LocalNz;k++) {
          BoutReal dy3 = SQ(coord->dy(i, j, k)) * coord->dy(i, j, k);
          // 3rd derivative at right boundary
          
          BoutReal d3fdx3 = (
                             f(i,j+2,k)
                             - 3.*f(i,j+1,k)
                             + 3.*f(i,j,  k)
                             -    f(i,j-1,k)
                             ) / dy3;

          BoutReal flux = 0.5 * (d(i, j, k) + d(i, j + 1, k))
                          * (coord->J(i, j, k) + coord->J(i, j + 1, k)) * d3fdx3;

          result(i, j, k) += flux / (coord->J(i, j, k) * coord->dy(i, j, k));
          result(i, j + 1, k) -= flux / (coord->J(i, j + 1, k) * coord->dy(i, j + 1, k));

          if(j == mesh->ystart && (!mesh->firstY())) {
            // Left cell boundary, no flux through boundaries
            d3fdx3 = (
                      f(i,j+1,k)
                      - 3.*f(i,j,  k)
                      + 3.*f(i,j-1,k)
                      -    f(i,j-2,k)
                      ) / dy3;

            flux = 0.5 * (d(i, j, k) + d(i, j - 1, k))
                   * (coord->J(i, j, k) + coord->J(i, j - 1, k)) * d3fdx3;

            result(i, j, k) -= flux / (coord->J(i, j, k) * coord->dy(i, j, k));
            result(i, j - 1, k) +=
                flux / (coord->J(i, j - 1, k) * coord->dy(i, j - 1, k));
          }
        }
      }
    
    // Convert result back to non-aligned coordinates
    return are_unaligned ? fromFieldAligned(result, "RGN_NOBNDRY") : result;
  }

  const Field3D D4DY4_Index(const Field3D &f_in, bool bndry_flux) {
    Mesh* mesh = f_in.getMesh();

    // Convert to field aligned coordinates
    const bool is_unaligned = (f_in.getDirectionY() == YDirectionType::Standard);
    Field3D f = is_unaligned ? toFieldAligned(f_in, "RGN_NOX") : f_in;

    Field3D result{zeroFrom(f)};

    Coordinates *coord = f_in.getCoordinates();
    
    for(int i=mesh->xstart;i<=mesh->xend;i++) {
      bool yperiodic = mesh->periodicY(i);
      
      bool has_upper_boundary = !yperiodic && mesh->lastY(i);
      bool has_lower_boundary = !yperiodic && mesh->firstY(i);
      
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
        
        // Right boundary
        
        if (bndry_flux || (j != mesh->yend) || !has_upper_boundary) {
          // Calculate the fluxes

          if (j != mesh->yend || !has_upper_boundary) {

            for (int k = 0; k < mesh->LocalNz; k++) {
              // Right boundary common factors
              const BoutReal common_factor =
                  0.25 * (coord->dy(i, j, k) + coord->dy(i, j + 1, k))
                  * (coord->J(i, j, j) + coord->J(i, j + 1, k));

              const BoutReal factor_rc =
                  common_factor / (coord->J(i, j, k) * coord->dy(i, j, k));
              const BoutReal factor_rp =
                  common_factor / (coord->J(i, j + 1, k) * coord->dy(i, j + 1, k));

              // Not on domain boundary
              // 3rd derivative at right cell boundary

              const BoutReal d3fdx3 = (f(i, j + 2, k) - 3. * f(i, j + 1, k)
                                       + 3. * f(i, j, k) - f(i, j - 1, k));

              result(i,j,  k) += d3fdx3 * factor_rc; 
              result(i,j+1,k) -= d3fdx3 * factor_rp; 
            }
          } else {
            // At a domain boundary
            // Use a one-sided difference formula
            
            for(int k=0;k<mesh->LocalNz;k++) {
              // Right boundary common factors
              const BoutReal common_factor =
                  0.25 * (coord->dy(i, j, k) + coord->dy(i, j + 1, k))
                  * (coord->J(i, j, j) + coord->J(i, j + 1, k));

              const BoutReal factor_rc =
                  common_factor / (coord->J(i, j, k) * coord->dy(i, j, k));
              const BoutReal factor_rp =
                  common_factor / (coord->J(i, j + 1, k) * coord->dy(i, j + 1, k));

              const BoutReal d3fdx3 =
                  -((16. / 5) * 0.5 * (f(i, j + 1, k) + f(i, j, k)) // Boundary value f_b
                    - 6. * f(i, j, k)                               // f_0
                    + 4. * f(i, j - 1, k)                           // f_1
                    - (6. / 5) * f(i, j - 2, k)                     // f_2
                  );

              result(i,j,  k) += d3fdx3 * factor_rc; 
              result(i,j+1,k) -= d3fdx3 * factor_rp;
            }
          }
        }

        // Left cell boundary
        
        if (bndry_flux || (j != mesh->ystart) || !has_lower_boundary) {
          // Calculate the fluxes

            
          if ( j != mesh->ystart || !has_lower_boundary ) {
            for(int k=0;k<mesh->LocalNz;k++) {
              const BoutReal common_factor =
                  0.25 * (coord->dy(i, j, k) + coord->dy(i, j + 1, k))
                  * (coord->J(i, j, k) + coord->J(i, j - 1, k));

              const BoutReal factor_lc =
                  common_factor / (coord->J(i, j, k) * coord->dy(i, j, k));
              const BoutReal factor_lm =
                  common_factor / (coord->J(i, j - 1, k) * coord->dy(i, j - 1, k));

              // Not on a domain boundary
              const BoutReal d3fdx3 = (f(i, j + 1, k) - 3. * f(i, j, k)
                                       + 3. * f(i, j - 1, k) - f(i, j - 2, k));

              result(i, j    , k) -= d3fdx3 * factor_lc; 
              result(i, j - 1, k) += d3fdx3 * factor_lm;
            }
          } else {
            // On a domain (Y) boundary
            for(int k=0;k<mesh->LocalNz;k++) {
              const BoutReal common_factor =
                  0.25 * (coord->dy(i, j, k) + coord->dy(i, j + 1, k))
                  * (coord->J(i, j, k) + coord->J(i, j - 1, k));

              const BoutReal factor_lc =
                  common_factor / (coord->J(i, j, k) * coord->dy(i, j, k));
              const BoutReal factor_lm =
                  common_factor / (coord->J(i, j - 1, k) * coord->dy(i, j - 1, k));
              const BoutReal d3fdx3 =
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
    return is_unaligned ? fromFieldAligned(result, "RGN_NOBNDRY") : result;
  }

  void communicateFluxes(Field3D &f) {
    Mesh* mesh = f.getMesh();

    // Use X=0 as temporary buffer
    if (mesh->xstart != 2)
      throw BoutException("communicateFluxes: Sorry!");

    int size = mesh->LocalNy * mesh->LocalNz;
    comm_handle xin, xout;
    // Cache results to silence spurious compiler warning about xin,
    // xout possibly being uninitialised when used
    bool not_first = !mesh->firstX();
    bool not_last = !mesh->lastX();
    if (not_first) {
      xin = mesh->irecvXIn(f(0, 0), size, 0);
    }
    if (not_last) {
      xout = mesh->irecvXOut(f(mesh->LocalNx - 1, 0), size, 1);
    }
    // Send X=1 values
    if (not_first) {
      mesh->sendXIn(f(1, 0), size, 1);
    }
    if (not_last) {
      mesh->sendXOut(f(mesh->LocalNx - 2, 0), size, 0);
    }
    // Wait
    if (not_first) {
      mesh->wait(xin);
      // Add to cells
      for (int y = mesh->ystart; y <= mesh->yend; y++)
        for (int z = 0; z < mesh->LocalNz; z++) {
          f(2, y, z) += f(0, y, z);
        }
    }
    if (not_last) {
      mesh->wait(xout);
      // Add to cells
      for (int y = mesh->ystart; y <= mesh->yend; y++)
        for (int z = 0; z < mesh->LocalNz; z++) {
          f(mesh->LocalNx - 3, y, z) += f(mesh->LocalNx - 1, y, z);
        }
    }
  }

  Field3D Div_Perp_Lap(const Field3D& a, const Field3D& f, CELL_LOC outloc) {

    Field3D result = 0.0;

    //////////////////////////////////////////
    // X-Z diffusion
    //
    //            Z
    //            |
    //
    //     o --- gU --- o
    //     |     nU     |
    //     |            |
    //    gL nL      nR gR    -> X
    //     |            |
    //     |     nD     |
    //     o --- gD --- o
    //
    Coordinates* coords = a.getCoordinates(outloc);
    Mesh* mesh = f.getMesh();

    for (int i = mesh->xstart; i <= mesh->xend; i++)
      for (int j = mesh->ystart; j <= mesh->yend; j++)
        for (int k = 0; k < mesh->LocalNz; k++) {

          // wrap k-index around as Z is (currently) periodic.
          int kp = (k + 1) % (mesh->LocalNz);
          int km = (k - 1 + mesh->LocalNz) % (mesh->LocalNz);

          // Calculate gradients on cell faces -- assumes constant grid spacing

          BoutReal gR =
              (coords->g11(i, j, k) + coords->g11(i + 1, j, k))
                  * (f(i + 1, j, k) - f(i, j, k))
                  / (coords->dx(i + 1, j, k) + coords->dx(i, j, k))
              + 0.5 * (coords->g13(i, j, k) + coords->g13(i + 1, j, k))
                    * (f(i + 1, j, kp) - f(i + 1, j, km) + f(i, j, kp) - f(i, j, km))
                    / (4. * coords->dz(i, j, k));

          BoutReal gL =
              (coords->g11(i - 1, j, k) + coords->g11(i, j, k))
                  * (f(i, j, k) - f(i - 1, j, k))
                  / (coords->dx(i - 1, j, k) + coords->dx(i, j, k))
              + 0.5 * (coords->g13(i - 1, j, k) + coords->g13(i, j, k))
                    * (f(i - 1, j, kp) - f(i - 1, j, km) + f(i, j, kp) - f(i, j, km))
                    / (4 * coords->dz(i, j, k));

          BoutReal gD =
              coords->g13(i, j, k)
                  * (f(i + 1, j, km) - f(i - 1, j, km) + f(i + 1, j, k) - f(i - 1, j, k))
                  / (4. * coords->dx(i, j, k))
              + coords->g33(i, j, k) * (f(i, j, k) - f(i, j, km)) / coords->dz(i, j, k);

          BoutReal gU =
              coords->g13(i, j, k)
                  * (f(i + 1, j, kp) - f(i - 1, j, kp) + f(i + 1, j, k) - f(i - 1, j, k))
                  / (4. * coords->dx(i, j, k))
              + coords->g33(i, j, k) * (f(i, j, kp) - f(i, j, k)) / coords->dz(i, j, k);

          // Flow right
          BoutReal flux = gR * 0.25 * (coords->J(i + 1, j, k) + coords->J(i, j, k))
                          * (a(i + 1, j, k) + a(i, j, k));
          result(i, j, k) += flux / (coords->dx(i, j, k) * coords->J(i, j, k));

          // Flow left
          flux = gL * 0.25 * (coords->J(i - 1, j, k) + coords->J(i, j, k))
                 * (a(i - 1, j, k) + a(i, j, k));
          result(i, j, k) -= flux / (coords->dx(i, j, k) * coords->J(i, j, k));

          // Flow up
          flux = gU * 0.25 * (coords->J(i, j, k) + coords->J(i, j, kp))
                 * (a(i, j, k) + a(i, j, kp));
          result(i, j, k) += flux / (coords->dz(i, j, k) * coords->J(i, j, k));

          // Flow down
          flux = gD * 0.25 * (coords->J(i, j, km) + coords->J(i, j, k))
                 * (a(i, j, km) + a(i, j, k));
          result(i, j, k) += flux / (coords->dz(i, j, k) * coords->J(i, j, k));
        }

    return result;
  }

} // Namespace FV
