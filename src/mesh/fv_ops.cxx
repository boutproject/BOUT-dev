
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

    Field3D result(mesh);
    result = 0.0;

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
	
	  BoutReal fout = (coord->J(i,j)*a(i,j,k)*coord->g11(i,j) + coord->J(i+1,j)*a(i+1,j,k)*coord->g11(i+1,j)) *
	    (f(i+1,j,k) - f(i,j,k))/(coord->dx(i,j) + coord->dx(i+1,j));
                     
	  result(i,j,k) += fout / (coord->dx(i,j)*coord->J(i,j));
	  result(i+1,j,k) -= fout / (coord->dx(i+1,j)*coord->J(i+1,j));
	}
      }
  
    // Y flux
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart-1;j<=mesh->yend;j++) {
      
	BoutReal coef = 0.5*(coord->g_23(i,j)/SQ(coord->J(i,j)*coord->Bxy(i,j)) + coord->g_23(i,j+1)/SQ(coord->J(i,j+1)*coord->Bxy(i,j+1)));
      
	for(int k=0;k<mesh->LocalNz;k++) {
	  // Calculate flux between j and j+1
	  int kp = (k + 1) % (mesh->LocalNz);
	  int km = (k - 1 + (mesh->LocalNz)) % (mesh->LocalNz);
          
	  // Calculate Z derivative at y boundary
	  BoutReal dfdz = 0.25*(f(i,j,kp) - f(i,j,km)
				+ f.yup()(i,j+1,kp) - f.yup()(i,j+1,km))/coord->dz;
	
	  // Y derivative
	  BoutReal dfdy = 2.*(f.yup()(i,j+1,k) - f(i,j,k))/(coord->dy(i,j+1) + coord->dy(i,j));
	
	  BoutReal fout = 0.5*(coord->J(i,j)*a(i,j,k)*coord->g23(i,j) + coord->J(i,j+1)*a.yup()(i,j+1,k)*coord->g23(i,j+1))
	    * ( dfdz - coef*dfdy );

	  result(i,j,k) += fout / (coord->dy(i,j)*coord->J(i,j));

          // Calculate flux between j and j-1
          dfdz = 0.25*(f(i,j,kp) - f(i,j,km)
                       + f.ydown()(i,j-1,kp) - f.ydown()(i,j-1,km))/coord->dz;
          
          dfdy = 2.*(f(i,j,k) - f.ydown()(i,j-1,k))/(coord->dy(i,j) + coord->dy(i,j-1));
          
          
          fout = 0.5*(coord->J(i,j)*a(i,j,k)*coord->g23(i,j) + coord->J(i,j-1)*a.yup()(i,j+1,k)*coord->g23(i,j+1))
	    * ( dfdz - coef*dfdy );

          result(i,j,k) -= fout / (coord->dy(i,j)*coord->J(i,j));
	}
      }
  
    // Z flux
    // Easier since all metrics constant in Z
  
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
	// Coefficient in front of df/dy term
	BoutReal coef = coord->g_23(i,j)
	  / (coord->dy(i,j+1) + 2.*coord->dy(i,j) + coord->dy(i,j-1))
	  / SQ(coord->J(i,j)*coord->Bxy(i,j));
	
	for(int k=0;k<mesh->LocalNz;k++) {
	  // Calculate flux between k and k+1
	  int kp = (k + 1) % (mesh->LocalNz);
	
	  BoutReal fout = 0.5*(a(i,j,k) + a(i,j,kp)) * coord->g33(i,j)*
	    ( 
	     // df/dz
	     (f(i,j,kp) - f(i,j,k))/coord->dz 
	   
	     // - g_yz * df/dy / SQ(J*B)
	     - coef*(f.yup()(i,j+1,k) + f.yup()(i,j+1,kp) - f.ydown()(i,j-1,k) - f.ydown()(i,j-1,kp))
	      );
	
	  result(i,j,k) += fout / coord->dz;
	  result(i,j,kp) -= fout / coord->dz;
	}
      }
  
    return result;
  }

  const Field3D Div_par_K_Grad_par(const Field3D &Kin, const Field3D &fin, bool bndry_flux) {
    TRACE("FV::Div_par_K_Grad_par");

    ASSERT2(Kin.getLocation() == fin.getLocation());

    Mesh *mesh = Kin.getMesh();
    Field3D result(0.0, mesh);

    bool use_yup_ydown = (Kin.hasYupYdown() && fin.hasYupYdown());

    const auto& K = use_yup_ydown ? Kin : mesh->toFieldAligned(Kin);
    const auto& f = use_yup_ydown ? fin : mesh->toFieldAligned(fin);

    // K and f fields in yup and ydown directions
    const auto& Kup = use_yup_ydown ? Kin.yup() : K;
    const auto& Kdown = use_yup_ydown ? Kin.ydown() : K;
    const auto& fup = use_yup_ydown ? fin.yup() : f;
    const auto& fdown = use_yup_ydown ? fin.ydown() : f;
    
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
        
        BoutReal g_22 = 0.5*(coord->g_22[i] + coord->g_22[i]);
        
        BoutReal gradient = 2.*(f[i] - fdown[iym]) / (coord->dy[i] + coord->dy[iym]);
        
        BoutReal flux = c * J * gradient / g_22;
        
        result[i] -= flux / (coord->dy[i] * coord->J[i]);
      }
    }
    
    if (!use_yup_ydown) {
      // Shifted to field aligned coordinates, so need to shift back
      result = mesh->fromFieldAligned(result);
    }
    
    return result;
  }

  const Field3D D4DY4(const Field3D &d_in, const Field3D &f_in) {
    ASSERT2(d_in.getLocation() == f_in.getLocation());

    Mesh* mesh = d_in.getMesh();
    ASSERT1(mesh = f_in.getMesh());

    Field3D result = 0.0;
    result.setLocation(f_in.getLocation());
    
    Coordinates *coord = f_in.getCoordinates();
    
    // Convert to field aligned coordinates
    Field3D d = mesh->toFieldAligned(d_in);
    Field3D f = mesh->toFieldAligned(f_in);
    
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
        BoutReal dy3 = SQ(coord->dy(i,j))*coord->dy(i,j);
        for(int k=0;k<mesh->LocalNz;k++) {
          
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
          
          if(j == mesh->ystart && (!mesh->firstY())) {
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
    return mesh->fromFieldAligned(result);
  }

  const Field3D D4DY4_Index(const Field3D &f_in, bool bndry_flux) {
    Field3D result = 0.0;
    result.setLocation(f_in.getLocation());
    
    Mesh* mesh = f_in.getMesh();

    // Convert to field aligned coordinates
    Field3D f = mesh->toFieldAligned(f_in);

    Coordinates *coord = f_in.getCoordinates();
    
    for(int i=mesh->xstart;i<=mesh->xend;i++) {
      bool yperiodic = mesh->periodicY(i);
      
      bool has_upper_boundary = !yperiodic && mesh->lastY(i);
      bool has_lower_boundary = !yperiodic && mesh->firstY(i);
      
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
        
        // Right boundary
        
        if (bndry_flux || (j != mesh->yend) || !has_upper_boundary) {
          // Calculate the fluxes
          
          // Right boundary common factors
          BoutReal common_factor = 0.25 *
            (coord->dy(i, j) + coord->dy(i, j + 1)) *
            (coord->J(i, j) + coord->J(i, j + 1));
          
          BoutReal factor_rc = common_factor / (coord->J(i,j) * coord->dy(i,j));
          BoutReal factor_rp = common_factor / (coord->J(i,j+1) * coord->dy(i,j+1));
          if ( j != mesh->yend || !has_upper_boundary ) {
            for(int k=0;k<mesh->LocalNz;k++) {
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
            
            for(int k=0;k<mesh->LocalNz;k++) {

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
        
        if (bndry_flux || (j != mesh->ystart) || !has_lower_boundary) {
          // Calculate the fluxes

          BoutReal common_factor = 0.25 *
            (coord->dy(i, j) + coord->dy(i, j + 1)) *
            (coord->J(i, j) + coord->J(i, j - 1));

          BoutReal factor_lc = common_factor / (coord->J(i, j) * coord->dy(i, j));
          BoutReal factor_lm = common_factor / (coord->J(i, j - 1) * coord->dy(i, j - 1));
            
          if ( j != mesh->ystart || !has_lower_boundary ) {
            for(int k=0;k<mesh->LocalNz;k++) {
          
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
            for(int k=0;k<mesh->LocalNz;k++) {
              BoutReal d3fdx3 = -(-(16. / 5) * 0.5 *
                                  (f(i, j - 1, k) + f(i, j, k)) // Boundary value f_b
                                  + 6. * f(i, j, k)                 // f_0
                                  - 4. * f(i, j + 1, k)             // f_1
                                  + (6. / 5) * f(i, j + 2, k)       // f_2
                                  );

              result(i, j    , k) -= d3fdx3 * factor_lc; 
              result(i, j - 1, k) += d3fdx3 * factor_lm;
            }
          }
        }
      }
    }
    
    // Convert result back to non-aligned coordinates
    return mesh->fromFieldAligned(result);
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

} // Namespace FV
