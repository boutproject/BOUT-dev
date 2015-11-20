
#include <bout/fv_ops.hxx>
#include <globals.hxx>
#include <utils.hxx>

#include <output.hxx>

namespace FV {

  // Div ( a Laplace_perp(f) )  -- Vorticity
  const Field3D Div_a_Laplace_perp(const Field3D &a, const Field3D &f) {
  
    Field3D result = 0.0;
  
    // Flux in x

    Field3D fs = f.shiftZ(true); // Orthogonal X-Z coords
    Field3D as = a.shiftZ(true);
  
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
	for(int k=0;k<mesh->ngz-1;k++) {
	  // Calculate flux from i to i+1
	
	  BoutReal fout = (mesh->J(i,j)*as(i,j,k)*mesh->g11(i,j) + mesh->J(i+1,j)*as(i+1,j,k)*mesh->g11(i+1,j)) *
	    (fs(i+1,j,k) - fs(i,j,k))/(mesh->dx(i,j) + mesh->dx(i+1,j));
                     
	  result(i,j,k) += fout / (mesh->dx(i,j)*mesh->J(i,j));
	  result(i+1,j,k) -= fout / (mesh->dx(i+1,j)*mesh->J(i+1,j));
	}
      }

    result = result.shiftZ(false); // Field aligned
  
    // Y flux
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart-1;j<=mesh->yend;j++) {
      
	BoutReal coef = 0.5*(mesh->g_23(i,j)/SQ(mesh->J(i,j)*mesh->Bxy(i,j)) + mesh->g_23(i,j+1)/SQ(mesh->J(i,j+1)*mesh->Bxy(i,j+1)));
      
	for(int k=0;k<mesh->ngz-1;k++) {
	  // Calculate flux between j and j+1
	  int kp = (k + 1) % (mesh->ngz-1);
	  int km = (k - 1 + (mesh->ngz-1)) % (mesh->ngz-1);
	
	  // Calculate Z derivative at y boundary
	  BoutReal dfdz = 0.25*(f(i,j,kp) - f(i,j,km)
				+ f(i,j+1,kp) - f(i,j+1,km))/mesh->dz;
	
	  // Y derivative
	  BoutReal dfdy = 2.*(f(i,j+1,k) - f(i,j,k))/(mesh->dy(i,j+1) + mesh->dy(i,j));
	
	  BoutReal fout = 0.5*(mesh->J(i,j)*a(i,j,k)*mesh->g23(i,j) + mesh->J(i,j+1)*a(i,j+1,k)*mesh->g23(i,j+1))
	    * ( dfdz - coef*dfdy );

	  result(i,j,k) += fout / (mesh->dy(i,j)*mesh->J(i,j));
	  result(i,j+1,k) -= fout / (mesh->dy(i,j+1)*mesh->J(i,j+1));
	}
      }
  
    // Z flux
    // Easier since all metrics constant in Z
  
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
	// Coefficient in front of df/dy term
	BoutReal coef = mesh->g_23(i,j)
	  / (mesh->dy(i,j+1) + 2.*mesh->dy(i,j) + mesh->dy(i,j-1))
	  / SQ(mesh->J(i,j)*mesh->Bxy(i,j));
	
	for(int k=0;k<mesh->ngz-1;k++) {
	  // Calculate flux between k and k+1
	  int kp = (k + 1) % (mesh->ngz-1);
	
	  BoutReal fout = 0.5*(a(i,j,k) + a(i,j,kp)) * mesh->g33(i,j)*
	    ( 
	     // df/dz
	     (f(i,j,kp) - f(i,j,k))/mesh->dz 
	   
	     // - g_yz * df/dy / SQ(J*B)
	     - coef*(f(i,j+1,k) + f(i,j+1,kp) - f(i,j-1,k) - f(i,j-1,kp))
	      );
	
	  result(i,j,k) += fout / mesh->dz;
	  result(i,j,kp) -= fout / mesh->dz;
	}
      }
  
    return result;
  }

  const Field3D Div_par_K_Grad_par(const Field3D &K, const Field3D &f, bool bndry_flux) {
    Field3D result;
    result = 0.0;
  
    for(int i=mesh->xstart;i<=mesh->xend;i++)
      for(int j=mesh->ystart-1;j<=mesh->yend;j++)
	for(int k=0;k<mesh->ngz-1;k++) {
	  // Calculate flux at upper surface
        
	  if(!bndry_flux) {
	    if((j == mesh->yend) && mesh->lastY())
	      continue;
        
	    if((j == mesh->ystart-1) && mesh->firstY())
	      continue;
	  }
	  BoutReal c = 0.5*(K(i,j,k) + K(i,j+1,k)); // K at the upper boundary
	  BoutReal J = 0.5*(mesh->J(i,j) + mesh->J(i,j+1)); // Jacobian at boundary
        
	  BoutReal g_22 = 0.5*(mesh->g_22(i,j) + mesh->g_22(i,j+1));
        
	  BoutReal gradient = 2.*(f(i,j+1,k) - f(i,j,k)) / (mesh->dy(i,j) + mesh->dy(i,j+1));
        
	  BoutReal flux = c * J * gradient / g_22;
        
	  result(i,j,k) += flux / (mesh->dy(i,j) * mesh->J(i,j));
	  result(i,j+1,k) -= flux / (mesh->dy(i,j+1) * mesh->J(i,j+1));
	}
    return result;
  }

} // Namespace FV
