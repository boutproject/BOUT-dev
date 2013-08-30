/*******************************************************************
        2D/3D Blob simulations 
 
         NR Walkden, 20 January 2012
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <derivs.hxx>
#include <invert_laplace.hxx>

/********* Global variables***********/

// Evolving variables
Field3D n,omega;                                // Density and vorticity

// Auxilliary variables
Field3D Jpar,phi;                               // Parallel current and electrostatic potential

//Parameters
BoutReal rho_s,Omega_i,c_s,n0;                  // Bohm gyro radius, Ion cyclotron frequency, Bohm sound speed

//Collisional parameters
BoutReal sigma_par,tau_ei;                      // Parallel conduction, electron ion collision time. 

//Constants to calculate the parameters
BoutReal Te0,e,B0,D_n,D_vort,R_c,L_par,m_i,m_e;

//Boundary conditions
bool sheath_limited;                            // Apply sheath limited boudnary condition if 3D

//Model options 
bool boussinesq;                                // Use the Boussinesq approximation in vorticity
bool compressible;                              // If allow inclusion of n grad phi term in density evolution

Laplacian *phiSolver;

int physics_init(bool restarting) { 
 
  /******************Reading options *****************/

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("model");
   
  // Load system parameters
  options->get("Te0", Te0, 30);         // Temp in eV
  options->get("e", e, 1.602e-19);   
  options->get("m_i",m_i,2*1.667e-27); 
  options->get("m_e",m_e,9.11e-31);
  options->get("n0",n0,1e19);           // Background density in cubic m
  options->get("D_vort",D_vort,0);      // Viscous diffusion coefficient
  options->get("D_n",D_n,0);            // Density diffusion coefficient
  options->get("R_c",R_c,1.5);          // Radius of curvature

  // System option switches

  OPTION(options,sheath_limited,false); // Use sheath limited boundary conditions
  OPTION(options,compressible,false);   // Include compressible ExB term in density equation
  OPTION(options,boussinesq,true);      // Use Boussinesq approximation in vorticity
 
  //More options
  OPTION(options, B0, 0.35);            // Value of magnetic field strength
  
  /***************Calculate the Parameters **********/
  
  Omega_i = e*B0/m_i;           //Cyclotron Frequency
  c_s = sqrt(e * Te0/m_i);      //Bohm sound speed
  rho_s = c_s/Omega_i;          //Bohm gyro-radius
  
  output.write("\n\n\t----------Parameters: ------------ \n\tOmega_i = %e /s,\n\tc_s = %e m/s,\n\trho_s = %e m\n",
               Omega_i, c_s, rho_s);

  tau_ei = 3.44e11 * pow(Te0,1.5)/10;           //Electron-ion collision time
  output.write("\ttau_ei = %e s\n", tau_ei/n0); 
  sigma_par = 1.96 * tau_ei*e*e/m_e;            //Parallel electorn conduction (approximate)  

  /************ Create a solver for potential ********/

  phiSolver = Laplacian::create(Options::getRoot()->getSection("phiSolver")); // BOUT.inp section "phiSolver"
  phi = 0.0;
  
  /************ Tell BOUT++ what to solve ************/
  
  bout_solve(n, "Density");
  bout_solve(omega, "Vorticity");  
 
  if(!sheath_limited){Jpar.setBoundary("Jpar");}
  phi.setBoundary("phi");
  
  // Output Jpar and phi  
  SAVE_REPEAT2(Jpar,phi);                
  SAVE_ONCE3(rho_s, c_s, Omega_i);
 
  return 0;
}

// just define a macro for V_E dot Grad
// poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Can use BRACKET_STD, BRACKET_ARAKAWA or BRACKET_SIMPLE
#define vE_Grad(f, p) ( bracket(p, f, BRACKET_SIMPLE) )

void Jpar_sheath_boundary();

int physics_run(BoutReal t) {
  
  // Run communications 
  ////////////////////////////////////////////////////////////////////////////
  mesh->communicate(n,omega);

  //Invert div(n grad(phi)) = grad(n) grad(phi) + n Delp_perp^2(phi) = omega
  //////////////////////////////////////////////////////////////////////////// 

  if(!boussinesq) {
    // Including full density in vorticit inversion
    phiSolver->setCoefC(n);  // Update the 'C' coefficient. See invert_laplace.hxx
    phi = phiSolver->solve(omega / n, phi);  // Use previous solution as guess
  }else {
    // Background density only (1 in normalised units)
    phi = phiSolver->solve(omega); 
  }
  
  phi.applyBoundary();  // Apply parallel (y) boundaries
  mesh->communicate(phi);

  // Calculate Jpar from parallel ohms Law
  ////////////////////////////////////////////////////////////////////////////
  
  Jpar = Grad_par(n)/n - Grad_par(phi);
  Jpar *= sigma_par*Te0/(rho_s*n0*c_s*e);
  if(sheath_limited) {Jpar_sheath_boundary();}
  else {Jpar.applyBoundary();} 
  
  mesh->communicate(Jpar);
  
  // Density Evolution
  /////////////////////////////////////////////////////////////////////////////
  
  ddt(n) = -vE_Grad(n,phi)                                              // ExB term
         + 2*DDZ(n)*(rho_s/R_c)                                         // Curvature term
         + D_n*Delp2(n)
    ;                                                                   // Diffusion term
  if(compressible){ddt(n) -= 2*n*DDZ(phi)*(rho_s/R_c);}                 // ExB Compression term
  ddt(n) += Grad_par(Jpar);                                             // 3D parallel loss term

  // Vorticity evolution
  /////////////////////////////////////////////////////////////////////////////

  ddt(omega) = -vE_Grad(omega,phi)                                      // ExB term
             + 2*DDZ(n)*(rho_s/R_c)/n                                   // Curvature term
             + D_vort*Delp2(omega)/n                                    // Viscous diffusion term
             + Grad_par(Jpar)/n                                         // 3D parallel loss term
    ;
  
    return 0;
}



void Jpar_sheath_boundary() { //Apply Sheath limited boundary conditions to Jpar

  RangeIterator idwn = mesh->iterateBndryLowerY();
  RangeIterator iup = mesh->iterateBndryUpperY(); 

  for(iup.first();!iup.isDone();iup.next()){
    for(int jy = mesh->yend+1;jy < mesh->ngy; jy++){
      for(int jz = 0; jz < mesh->ngz;jz++){     
        
        Jpar(iup.ind, jy, jz) = n(iup.ind, jy, jz)*(1-exp(-phi(iup.ind, jy, jz)));
        
      }
    }
  }
  for(idwn.first();!idwn.isDone();idwn.next()){
    for(int jy = mesh->ystart-1;jy >= 0; jy--){
      for(int jz = 0; jz < mesh->ngz;jz++){
        
        Jpar(idwn.ind, jy, jz) = n(idwn.ind, jy, jz)*(1-exp(phi(idwn.ind, jy, jz)));
        
      }
    }
  }
}

