/// 3D simulations of HW
/////
///// This version uses indexed operators
///// which reduce the number of loops over the domain
/////
////  GPU processing is enabled if BOUT_ENABLE_CUDA is defined
/////  Profiling markers and ranges are set if USE_NVTX is defined
/////  Baesed on Ben Duddson, Steven Glenn code, Yining Qin update 0521-2020

#include <iostream>
#include <cstdlib>

/// 3D simulations of HW
/////
///// This version uses indexed operators
///// which reduce the number of loops over the domain
/////
////  GPU processing is enabled if BOUT_ENABLE_CUDA is defined
/////  Profiling markers and ranges are set if USE_NVTX is defined
/////  Baesed on Ben Duddson, Steven Glenn code, Yining Qin update 0521-2020


#include <bout/single_index_ops.hxx>

//	#include <gpu_functions.hxx>

#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>
#include <bout/single_index_ops.hxx>
#include "RAJA/RAJA.hpp" // using RAJA lib
#include <cuda_profiler_api.h>

#if defined(BOUT_USE_CUDA) && defined(__CUDACC__)
#define BOUT_HOST_DEVICE __host__ __device__
#define BOUT_HOST __host__
#define BOUT_DEVICE __device__
#else
#define BOUT_HOST_DEVICE
#define BOUT_HOST
#define BOUT_DEVICE
#endif


#ifdef BOUT_USE_CUDA
const int CUDA_BLOCK_SIZE = 256;  // TODO: Make configurable
using EXEC_POL = RAJA::cuda_exec<CUDA_BLOCK_SIZE>;
#else   
using EXEC_POL = RAJA::loop_exec;
#endif  


__managed__  BoutReal* gpu_n_ddt; // copy ddt(n) to __device__


class HW3D : public PhysicsModel {
public:
  Field3D n, vort;  // Evolving density and vorticity
  Field3D phi;      // Electrostatic potential

  // Model parameters
  BoutReal alpha;      // Adiabaticity (~conductivity)
  BoutReal kappa;      // Density gradient drive
  BoutReal Dvort, Dn;  // Diffusion
  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver for vort -> phi

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
    // Solve for potential
    phi = phiSolver->solve(vort, phi);    
    Field3D phi_minus_n = phi - n;
    
    // Communicate variables
    mesh->communicate(n, vort, phi, phi_minus_n);

    // Create accessors which enable fast access
    auto n_acc = FieldAccessor<>(n);
    auto vort_acc = FieldAccessor<>(vort);
    auto vort_acc_lite = FieldAccessorLite<>(vort);
    auto phi_acc = FieldAccessor<>(phi);
    auto phi_minus_n_acc = FieldAccessor<>(phi_minus_n);

	auto indices = n.getRegion("RGN_NOBNDRY").getIndices();
	Ind3D *ob_i = &(indices[0]);

gpu_n_ddt= const_cast<BoutReal*>(n.timeDeriv()->operator()(0,0)); // copy ddt(n) to __device__


///*
//  RAJA GPU code ----------- start
 //   RAJA_data_copy(n,vort, phi,phi_minus_n,phi_acc,phi_minus_n_acc);
    
    RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, n.getRegion("RGN_NOBNDRY").getIndices().size()), [=] BOUT_HOST_DEVICE (int i) {	

		
	auto ind = ob_i[i];


	//	BoutReal test = 0;
   //   BoutReal test1 =- kappa * DDZ_g(phi_acc, ind);
	//BoutReal test2 = - Delp2_gt(vort_acc, ind) ;
	BoutReal test2 = - Delp2_gt(vort_acc_lite,ind) ;
   if(i < 16) {
      printf("test2 = %f\n",test2);
   }
	//gpu_n_ddt[i] = test2 ;

              // BoutReal test = - kappa * DDZ(phi_acc, id); 
	

	//	ddt(n)[id]= test;
	//	ddt(vort)[id]= -test;


		/*
		//auto test =  t_Div(phi_minus_n_acc, id);

	 // 	BoutReal div_current = alpha * t_Div(phi_minus_n_acc, id);
	  	
		BoutReal div_current = alpha * Div_par_Grad_par(phi_minus_n_acc, i);
		
		ddt(n)[i]= - bracket(phi_acc, n_acc, i)
			          - div_current
         	                  - kappa * DDZ(phi_acc, i)
				 // + kappa * t_DDZ(phi_acc,id)
                                  + Dn * Delp2(n_acc, i)
				 ;	

                ddt(vort)[i]= - bracket(phi_acc, vort_acc, i)
                                     - div_current
                                     + Dvort * Delp2 (vort_acc, i) 		
		;   

		*/
	});
//  RAJA GPU code ----------- end
//*/

/*
    BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {
      BoutReal div_current = alpha * Div_par_Grad_par(phi_minus_n_acc, i);
      ddt(n)[i] = -bracket(phi_acc, n_acc, i)
                  - div_current
               //   - kappa * DDZ(phi_acc, i)
                  + Dn * Delp2(n_acc, i)
		;

      ddt(vort)[i] = -bracket(phi_acc, vort_acc, i)
                     - div_current
                     + Dvort * Delp2(vort_acc, i)
		;
    }
*/


    return 0;
  }  // end RHS
}; // end class HW3D

// Define a main() function
BOUTMAIN(HW3D);

