
class GBS;

#ifndef __GBS_H__
#define __GBS_H__

#include <bout/physicsmodel.hxx>

#include <invert_laplace.hxx>
#include <bout/constants.hxx>

class GBS : public PhysicsModel {
protected:
  int init(bool restarting);
  int rhs(BoutReal t);
private:
  // Evolving variables
  Field3D Ne, Te;     // Electron density and temperature
  Field3D VePsi;      // Parallel electron velocity
  Field3D Vi;         // Parallel ion velocity
  Field3D Vort;       // Vorticity

  // Auxilliary variables
  Field3D phi;        // Electrostatic potential
  Field3D psi, Ve;

  // Sources of density and energy
  Field3D Sn, Sp;

  // Stress tensor
  Field3D Gi, Ge;
  
  // Collisional damping
  Field3D nu; 

  // Curvature terms
  int curv_method; // Determines which method is used
  Vector2D bxcv;    // b x kappa = (B/2)Curl(b/B)
  Field3D logB;     // Used in bracket method

  // Switches
  bool evolve_Ne, evolve_Vort, evolve_Ve, evolve_Vi, evolve_Te;
  bool ionvis, elecvis, resistivity;
  bool parallel; // Include parallel dynamics?
  bool estatic; // Electrostatic
  
  // Bracket method for advection terms // Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE      
  /*
   * Bracket method 
   *
   * BRACKET_STD      - Same as b0xGrad_dot_Grad, methods in BOUT.inp
   * BRACKET_SIMPLE   - Subset of terms, used in BOUT-06
   * BRACKET_ARAKAWA  - Arakawa central differencing (2nd order)
   * BRACKET_CTU      - 1st order upwind method
   *
   */
  BRACKET_METHOD bm_exb;

  bool mms;

  Coordinates *coords;
  
  // Stress tensor components
  BoutReal tau_e0, tau_i0;
  BoutReal Ti; // Ion temperature [eV] for stress tensor

  // Diffusion parameters
  BoutReal Dn, Dvort, Dve, Dvi, Dte;
  
  // Numerical hyper-diffusion parameters
  BoutReal Hn, Hvort, Hve, Hvi, Hte;

  // Normalisation parameters
  BoutReal Tnorm, Nnorm, Bnorm, Rnorm;
  BoutReal AA, Cs0, rho_s0, Omega_ci;
  BoutReal mi_me, beta_e;

  // Group of fields for communication
  FieldGroup evars; // Evolving variables
  
  // Initialisation
  void LoadMetric(BoutReal Lnorm, BoutReal Bnorm);

  // Operators
  const Field3D C(const Field3D &f); // Curvature operator
  const Field3D D(const Field3D &f, BoutReal d); // Diffusion operator
  const Field3D H(const Field3D &f, BoutReal h); // Hyper-diffusion
  // Powers of the mesh spacing for H operator
  Field2D dx4, dy4;
  BoutReal dz4;
  
  // Laplacian solver
  std::unique_ptr<Laplacian> phiSolver{nullptr};
  std::unique_ptr<Laplacian> aparSolver{nullptr};
};

#endif // __GBS_H__
