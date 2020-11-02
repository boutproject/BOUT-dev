/// 3D simulations of HW
/////
///// This version uses indexed operators
///// which reduce the number of loops over the domain
/////
////  GPU processing is enabled if BOUT_ENABLE_CUDA is defined
/////  Profiling markers and ranges are set if USE_NVTX is defined
/////  Baesed on Ben Duddson, Steven Glenn code, Yining Qin update 0521-2020

#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>

#include <bout/single_index_ops.hxx>

class HW3D : public PhysicsModel {
private:
  Field3D n, vort;  
  Field3D phi;      

 
  BoutReal alpha;      
  BoutReal kappa;      
  BoutReal Dvort, Dn;  
  
  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver for vort -> phi

protected:
  int init(bool UNUSED(restart)) {

    auto& options = Options::root()["hw"];
    alpha = options["alpha"].withDefault(1.0);
    kappa = options["kappa"].withDefault(0.1);
    Dvort = options["Dvort"].doc("Vorticity diffusion (normalised)").withDefault(1e-2);
    Dn = options["Dn"].doc("Density diffusion (normalised)").withDefault(1e-2);

    SOLVE_FOR(n, vort);
    SAVE_REPEAT(phi);
    
    phiSolver = Laplacian::create();
    phi = 0.; 
    
    return 0;
  }

  int rhs(BoutReal UNUSED(time)) {
   
    phi = phiSolver->solve(vort, phi);
    
    Field3D phi_minus_n = phi - n;
    
   
    mesh->communicate(n, vort, phi, phi_minus_n);

    
    auto n_acc = FieldAccessor<>(n);
    auto vort_acc = FieldAccessor<>(vort);
    auto phi_acc = FieldAccessor<>(phi);
    auto phi_minus_n_acc = FieldAccessor<>(phi_minus_n);
    
    BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {

      BoutReal div_current = alpha * Div_par_Grad_par(phi_minus_n_acc, i);

	div_current = 0.008;      
      ddt(n)[i] = -bracket(phi_acc, n_acc, i)
                  - div_current
                  - kappa * DDZ(phi_acc, i)
                  + Dn * Delp2(n_acc, i);

      ddt(vort)[i] = -bracket(phi_acc, vort_acc, i)
                     - div_current
                     + Dvort * Delp2(vort_acc, i);
    }

    return 0;
  }
};

// Define a main() function
 BOUTMAIN(HW3D);
