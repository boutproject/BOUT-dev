/*!******************************************************************
 * \file blob2d.cxx
 *
 *       2D simulations
 *
 *        NR Walkden, B Dudson  20 January 2012
 *******************************************************************/

#include <bout/physicsmodel.hxx> // Commonly used BOUT++ components
#include <derivs.hxx>            // To use DDZ()
#include <invert_laplace.hxx>    // Laplacian inversion

#include "bout/single_index_ops.hxx" // Operators at a single index

#include "RAJA/RAJA.hpp" // using RAJA lib
#include <cuda_profiler_api.h>
#include <gpu_functions.hxx>

/// 2D drift-reduced model, mainly used for blob studies
///
///
class Blob2D : public PhysicsModel {
public:
  // Evolving variables
  Field3D n, vort; ///< Density and vorticity

  // Auxilliary variables
  Field3D phi; ///< Electrostatic potential

  // Parameters
  BoutReal rho_s;   ///< Bohm gyro radius
  BoutReal Omega_i; ///< Ion cyclotron frequency
  BoutReal c_s;     ///< Bohm sound speed
  BoutReal n0;      ///< Reference density

  // Constants to calculate the parameters
  BoutReal Te0; ///< Isothermal temperature [eV]
  BoutReal B0;  ///< Constant magnetic field [T]
  BoutReal m_i; ///< Ion mass [kg]
  BoutReal m_e; ///< Electron mass [kg]
  BoutReal e;   ///< Electron charge [C]

  BoutReal D_n, D_vort; ///< Diffusion coefficients
  BoutReal R_c;         ///< Radius of curvature
  BoutReal L_par;       ///< Parallel connection length

  // Model options
  bool boussinesq;   ///< Use the Boussinesq approximation in vorticity
  bool compressible; ///< If allow inclusion of n grad phi term in density evolution
  bool sheath;       ///< Sheath connected?

  std::unique_ptr<Laplacian> phiSolver{nullptr}; ///< Performs Laplacian inversions to calculate phi

  int init(bool UNUSED(restarting)) {

    /******************Reading options *****************/

    auto& globalOptions = Options::root();
    auto& options = globalOptions["model"];

    // Load system parameters
    Te0 = options["Te0"].doc("Temperature in eV").withDefault(30.0);
    
    e = options["e"].withDefault(1.602e-19);
    m_i = options["m_i"].withDefault(2 * 1.667e-27);
    m_e = options["m_e"].withDefault(9.11e-31);

    n0 = options["n0"].doc("Background density in cubic m").withDefault(1e19);
    D_vort = options["D_vort"].doc("Viscous diffusion coefficient").withDefault(0.0);
    D_n = options["D_n"].doc("Density diffusion coefficient").withDefault(0.0);

    R_c = options["R_c"].doc("Radius of curvature").withDefault(1.5);
    L_par = options["L_par"].doc("Parallel connection length").withDefault(10.0);

    B0 = options["B0"].doc("Value of magnetic field strength").withDefault(0.35);

    // System option switches

    compressible = options["compressible"]
                       .doc("Compressible ExB term in density equation")
                       .withDefault(false);
    boussinesq = options["boussinesq"]
                     .doc("Use Boussinesq approximation in vorticity")
                     .withDefault(true);
    sheath = options["sheath"].doc("Sheath closure").withDefault(true);

    /***************Calculate the Parameters **********/

    Omega_i = e * B0 / m_i;    // Cyclotron Frequency
    c_s = sqrt(e * Te0 / m_i); // Bohm sound speed
    rho_s = c_s / Omega_i;     // Bohm gyro-radius

    output.write("\n\n\t----------Parameters: ------------ \n\tOmega_i = {:e} /s,\n\tc_s = "
                 "{:e} m/s,\n\trho_s = {:e} m\n",
                 Omega_i, c_s, rho_s);

    // Calculate delta_*, blob size scaling
    output.write("\tdelta_* = rho_s * (dn/n) * {:e} ",
                 pow(L_par * L_par / (R_c * rho_s), 1. / 5));

    /************ Create a solver for potential ********/

    if (boussinesq) {
       // BOUT.inp section "phiBoussinesq"
      phiSolver = Laplacian::create(Options::getRoot()->getSection("phiBoussinesq"));
    } else {
      // BOUT.inp section "phiSolver"
      phiSolver = Laplacian::create(Options::getRoot()->getSection("phiSolver")); 
    }
    phi = 0.0; // Starting guess for first solve (if iterative)

    /************ Tell BOUT++ what to solve ************/

    SOLVE_FOR(n, vort);

    // Output phi
    SAVE_REPEAT(phi);
    SAVE_ONCE(rho_s, c_s, Omega_i);

    return 0;
  }

  int rhs(BoutReal UNUSED(t)) {


    // RAJA::synchronize<RAJA::cuda_synchronize>();
    // Run communications
    ////////////////////////////////////////////////////////////////////////////
    mesh->communicate(n, vort);

    // Invert div(n grad(phi)) = grad(n) grad(phi) + n Delp_perp^2(phi) = omega
    ////////////////////////////////////////////////////////////////////////////

    if (!boussinesq) {
      // Including full density in vorticit inversion
      phiSolver->setCoefC(n); // Update the 'C' coefficient. See invert_laplace.hxx
      phi = phiSolver->solve(vort / n, phi); // Use previous solution as guess
    } else {
      // Background density only (1 in normalised units)
      phi = phiSolver->solve(vort);
    }

    mesh->communicate(phi);

    // Create data accessors for fast inner loop
    auto n_acc = FieldAccessor<>(n);
    auto vort_acc = FieldAccessor<>(vort);
    auto phi_acc = FieldAccessor<>(phi);

// / GPU code --------------------------start   

RAJA::synchronize<RAJA::cuda_synchronize>();
//BoutReal*   gpu_n_ddt= const_cast<BoutReal*>(n.timeDeriv()->operator()(0,0)); // copy ddt(n) to __device__
//BoutReal*   gpu_vort_ddt = const_cast<BoutReal*>(vort.timeDeriv()->operator()(0,0)); // copy ddt(vort) to __device__
 
 
  gpu_n_ddt= const_cast<BoutReal*>(n.timeDeriv()->operator()(0,0)); // copy ddt(n) to __device__
  gpu_vort_ddt = const_cast<BoutReal*>(vort.timeDeriv()->operator()(0,0)); // copy ddt(vort) to __device__

   auto region = n.getRegion("RGN_NOBNDRY"); // Region object
   auto indices = region.getIndices();   // A std::vector of Ind3D objects
   
//  Copy data to __device__
   RAJA_data_copy(n,vort,phi, phi,n_acc,phi_acc);

    RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, indices.size()), [=] RAJA_DEVICE (int i) {

                gpu_n_ddt[i] = - gpu_bracket_par(phi_acc, n_acc, i) //gpu_n_ddt[i]
                                  - 2 * gpu_DZZ_par(n_acc, i) * (rho_s / R_c)
                                  + D_n *gpu_Delp2_par(n_acc, i) 
				;

                gpu_vort_ddt[i] = - gpu_bracket_par(phi_acc, vort_acc, i)
                                     + 2 * gpu_DZZ_par(n_acc, i) * (rho_s / R_c)/ p_n[i]
                                     + D_vort *gpu_Delp2_par(vort_acc, i) / p_n[i]
				 ;

        });

//cudaDeviceSynchronize();

//RAJA::synchronize<RAJA::cuda_synchronize>();
//std::cout<<"gpu_n_ddt[99]:  "<< gpu_n_ddt[99]<< std::endl;
//std::cout<<"gpu_vort_ddt[99]:  "<< gpu_vort_ddt[99]<< std::endl;



// GPU code --------------------------end


  /*
    // Allocate arrays to store the time derivatives
    ddt(n).allocate();
    ddt(vort).allocate();
    // Iterate over the mesh except boundaries because derivatives are needed
    BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {

      // Density Evolution
      /////////////////////////////////////////////////////////////////////////////

      ddt(n)[i] = -bracket(phi_acc, n_acc, i)             // ExB term
                  + 2 * DDZ(n_acc, i) * (rho_s / R_c) // Curvature term
                  + D_n * Delp2(n_acc, i);            // Diffusion term
      
      // Vorticity evolution
      /////////////////////////////////////////////////////////////////////////////

      ddt(vort)[i] = -bracket(phi_acc, vort_acc, i)                // ExB term
                   + 2 * DDZ(n_acc, i) * (rho_s / R_c) / n_acc[i] // Curvature term
                   + D_vort * Delp2(vort_acc, i) / n_acc[i]      // Viscous diffusion term
          ;
    

	}

    if (compressible) {
      BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {
        ddt(n)[i] -= 2 * n_acc[i] * DDZ(phi_acc, i) * (rho_s / R_c); // ExB Compression term
      }
    }
    
    if (sheath) {
      // Sheath closure
      BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {
        ddt(n)[i] += n_acc[i] * phi_acc[i] * (rho_s / L_par);

        ddt(vort)[i] += phi_acc[i] * (rho_s / L_par);
      }
    }
 */ 
 
    return 0;
  }
};

// Define a standard main() function
BOUTMAIN(Blob2D);
