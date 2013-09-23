/*******************************************************************
        2D simulations 
 
         NR Walkden, B Dudson  20 January 2012
 *******************************************************************/

#include <bout.hxx>                // Commonly used BOUT++ components
#include <boutmain.hxx>            // Defines a simple main() function
#include <derivs.hxx>              // To use DDZ()
#include <invert_laplace.hxx>      // Laplacian inversion

/********* Global variables***********/

// Evolving variables
Field3D n,omega;                                // Density and vorticity

// Auxilliary variables
Field3D phi;                                    // Electrostatic potential

//Parameters
BoutReal rho_s,Omega_i,c_s,n0;                  // Bohm gyro radius, Ion cyclotron frequency, Bohm sound speed

//Constants to calculate the parameters
BoutReal Te0,e,B0,D_n,D_vort,R_c,L_par,m_i,m_e;

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
  OPTION(options, B0, 0.35);            // Value of magnetic field strength

  // System option switches
  
  OPTION(options,compressible,false);   // Include compressible ExB term in density equation
  OPTION(options,boussinesq,true);      // Use Boussinesq approximation in vorticity
  
  /***************Calculate the Parameters **********/
  
  Omega_i = e*B0/m_i;           //Cyclotron Frequency
  c_s = sqrt(e * Te0/m_i);      //Bohm sound speed
  rho_s = c_s/Omega_i;          //Bohm gyro-radius
  
  output.write("\n\n\t----------Parameters: ------------ \n\tOmega_i = %e /s,\n\tc_s = %e m/s,\n\trho_s = %e m\n",
               Omega_i, c_s, rho_s);

  /************ Create a solver for potential ********/

  if(boussinesq) {
    phiSolver = Laplacian::create(Options::getRoot()->getSection("phiBoussinesq")); // BOUT.inp section "phiBoussinesq"
  }else {
    phiSolver = Laplacian::create(Options::getRoot()->getSection("phiSolver")); // BOUT.inp section "phiSolver"
  }
  phi = 0.0; // Starting guess for first solve (if iterative)
  
  /************ Tell BOUT++ what to solve ************/
  
  SOLVE_FOR2(n, omega);
  
  // Output phi  
  SAVE_REPEAT(phi);                
  SAVE_ONCE3(rho_s, c_s, Omega_i);
 
  return 0;
}

// just define a macro for V_E dot Grad
// poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Can use BRACKET_STD, BRACKET_ARAKAWA or BRACKET_SIMPLE
#define vE_Grad(f, p) ( bracket(p, f, BRACKET_SIMPLE) )

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
  
  mesh->communicate(phi);
  
  // Density Evolution
  /////////////////////////////////////////////////////////////////////////////
  
  ddt(n) = -vE_Grad(n,phi)                                              // ExB term
         + 2*DDZ(n)*(rho_s/R_c)                                         // Curvature term
         + D_n*Delp2(n)
    ;                                                                   // Diffusion term
  if(compressible){ddt(n) -= 2*n*DDZ(phi)*(rho_s/R_c);}                 // ExB Compression term

  // Vorticity evolution
  /////////////////////////////////////////////////////////////////////////////

  ddt(omega) = -vE_Grad(omega,phi)                                      // ExB term
             + 2*DDZ(n)*(rho_s/R_c)/n                                   // Curvature term
             + D_vort*Delp2(omega)/n                                    // Viscous diffusion term
    ;
  
  return 0;
}
