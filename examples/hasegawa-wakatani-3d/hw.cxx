/// 3D simulations of HW
/////
///// This version uses indexed operators
///// which reduce the number of loops over the domain
/////
////  GPU processing is enabled if BOUT_ENABLE_CUDA is defined
/////  Profiling markers and ranges are set if USE_NVTX is defined
/////  Based on Ben Dudson, Steven Glenn code, Yining Qin update 0521-2020

#include <iostream>
#include <cstdlib>
#include <stdio.h>
/// 3D simulations of HW
/////
///// This version uses indexed operators
///// which reduce the number of loops over the domain
/////
/////  Profiling markers and ranges are set if USE_NVTX is defined
/////  Baesed on Ben Duddson, Steven Glenn code, Yining Qin update 0521-2020


#include <bout/single_index_ops.hxx>
#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>

#ifdef BOUT_HAS_RAJA
#include "RAJA/RAJA.hpp" // using RAJA lib
#endif

#if defined(BOUT_USE_CUDA) && defined(__CUDACC__)
#include <cuda_profiler_api.h>
#endif



#define MAX_n_region 100000


class HW3D : public PhysicsModel {
public:
  Field3D n, vort;  // Evolving density and vorticity
  Field3D phi;      // Electrostatic potential
  // Model parameters
  BoutReal alpha;      // Adiabaticity (~conductivity)
  BoutReal kappa;      // Density gradient drive
  BoutReal Dvort, Dn;  // Diffusion
  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver for vort -> phi
   
  int* public_obj_d;
  

 int init(bool UNUSED(restart)) {

    auto& options = Options::root()["hw"];
    alpha = options["alpha"].withDefault(1.0);
    kappa = options["kappa"].withDefault(0.1);
    Dvort = options["Dvort"].doc("Vorticity diffusion (normalised)").withDefault(1e-2);
    Dn = options["Dn"].doc("Density diffusion (normalised)").withDefault(1e-2);

    SOLVE_FOR(n, vort);
    SAVE_REPEAT(phi);
    
    phiSolver = Laplacian::create();
    phi = 0.; // Starting phi
    return 0;
  }



 int rhs(BoutReal UNUSED(time)) {

	//auto start = std::chrono::steady_clock::now();   
    
	// Solve for potential
    phi = phiSolver->solve(vort, phi);    
    Field3D phi_minus_n = phi - n;
    // Communicate variables
    mesh->communicate(n, vort, phi, phi_minus_n);

    // Create accessors which enable fast access
    auto n_acc = FieldAccessor<>(n);
    auto vort_acc = FieldAccessor<>(vort);
    auto phi_acc = FieldAccessor<>(phi);
    auto phi_minus_n_acc = FieldAccessor<>(phi_minus_n); 

#if 1  //def BOUT_HAS_RAJA

    //  RAJA code ----------- start

    
    Array<int> testArray(MAX_n_region); 
    auto obj_d = testArray.begin();
    auto  indices = n.getRegion("RGN_NOBNDRY").getIndices();
    int indices_size =indices.size();
    Ind3D *ob_i = &(indices)[0];
    for(int i= 0; i< indices_size;i++){
	obj_d[i]=ob_i[i].ind;
     }
 
 	BoutReal  alpha_UM  = alpha;  
	BoutReal  kappa_UM = kappa; 
	BoutReal  Dn_UM = Dn;
	BoutReal  Dvort_UM = Dvort;

    RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0,n.getRegion("RGN_NOBNDRY").size()), [=] RAJA_DEVICE (int id) {
    int i = obj_d[id];
      
	BoutReal div_current = alpha_UM * Div_par_Grad_par_g(phi_minus_n_acc, i);
		DDT(n_acc)[i] =  - bracket_g(phi_acc, n_acc, i)
                	    - div_current
                	    - kappa_UM * DDZ_g(phi_acc, i)
                	    + Dn_UM * Delp2_g(n_acc, i)
			;

		DDT(vort_acc)[i] = - bracket_g(phi_acc, vort_acc, i)
			      - div_current
		              + Dvort_UM * Delp2_g (vort_acc, i)
			;
	  });


//  RAJA code ----------- end
#else
// -- CPU code ------------- start
   
    BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {
      
      BoutReal div_current = alpha * Div_par_Grad_par(phi_minus_n_acc, i);

      ddt(n)[i] = -bracket(phi_acc, n_acc, i)
                  - div_current
                  - kappa * DDZ(phi_acc, i)
                  + Dn * Delp2(n_acc, i)
		;

      ddt(vort)[i] = -bracket(phi_acc, vort_acc, i)
                     - div_current
                     + Dvort * Delp2(vort_acc, i)
		;
    }
   
    //--CPU code ---------------- end
#endif
//auto end = std::chrono::steady_clock::now();
//auto  time_taken = std::chrono::duration_cast<std::chrono::nanoseconds     >(end-     start);
//std::cout << "The operator since start is "<< time_taken.count()<<" nano      se     conds.\n";

    return 0;
  }  // end RHS

}; // end class HW3D
 // Define a main() function
BOUTMAIN(HW3D);

