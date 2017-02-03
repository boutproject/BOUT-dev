
#include <bout/fv_ops.hxx>
#include <globals.hxx>
#include <utils.hxx>

#include <output.hxx>

namespace FV {

  // Div ( a Laplace_perp(f) )  -- Vorticity
  const Field3D Div_a_Laplace_perp(const Field3D &a, const Field3D &f) {
  
    Field3D result = 0.0;

    Coordinates *coord = mesh->coordinates();
    
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

  const Field3D Div_par_K_Grad_par(const Field3D &K, const Field3D &f, bool bndry_flux) {
    Field3D result;
    result = 0.0;
    
    Coordinates *coord = mesh->coordinates();
    
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++)
	for(int k=0;k<mesh->LocalNz;k++) {
          // Calculate flux at upper surface
        
	  if(bndry_flux || !mesh->lastY() || (j != mesh->yend)) {
	    
            BoutReal c = 0.5*(K(i,j,k) + K.yup()(i,j+1,k)); // K at the upper boundary
            BoutReal J = 0.5*(coord->J(i,j) + coord->J(i,j+1)); // Jacobian at boundary
            
            BoutReal g_22 = 0.5*(coord->g_22(i,j) + coord->g_22(i,j+1));
            
            BoutReal gradient = 2.*(f.yup()(i,j+1,k) - f(i,j,k)) / (coord->dy(i,j) + coord->dy(i,j+1));
        
            BoutReal flux = c * J * gradient / g_22;
            
            result(i,j,k) += flux / (coord->dy(i,j) * coord->J(i,j));
          }

          // Calculate flux at lower surface
          if(bndry_flux || !mesh->firstY() || (j != mesh->ystart)) {
            BoutReal c = 0.5*(K(i,j,k) + K.ydown()(i,j-1,k)); // K at the lower boundary
            BoutReal J = 0.5*(coord->J(i,j) + coord->J(i,j-1)); // Jacobian at boundary
            
            BoutReal g_22 = 0.5*(coord->g_22(i,j) + coord->g_22(i,j+1));
            
            BoutReal gradient = 2.*(f(i,j,k) - f.ydown()(i,j-1,k)) / (coord->dy(i,j) + coord->dy(i,j-1));
        
            BoutReal flux = c * J * gradient / g_22;
            
            result(i,j,k) -= flux / (coord->dy(i,j) * coord->J(i,j));
	  }
	}
    return result;
  }

} // Namespace FV
