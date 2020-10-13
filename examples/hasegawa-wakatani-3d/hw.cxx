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
#include <gpu_functions.hxx>

#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>
#include <bout/single_index_ops.hxx>
#include "RAJA/RAJA.hpp" // using RAJA lib
#include <cuda_profiler_api.h>


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
    auto phi_acc = FieldAccessor<>(phi);
    auto phi_minus_n_acc = FieldAccessor<>(phi_minus_n);

    auto region = n.getRegion("RGN_NOBNDRY"); // Region object
    auto indices = region.getIndices();   // A std::vector of Ind3D objects
//  Copy data to __device__
    RAJA_data_copy(n,vort, phi,phi_minus_n,phi_acc,phi_minus_n_acc);


//  GPU loop RAJA_DEVICE 
    {
       ArrayData<double> d2(10);
       std::iota(d2.begin(), d2.end(), 20);

       ArrayData<double> d3(10);
       std::iota(d3.begin(), d3.end(),30);

       Array<double> dd(10);
       dd[0] = 33.0;

       Array<double> dd2(10);
       dd2[0] = 34.0;

       Array<double> dd3(10);
       dd3[0] = 0.0;

       //RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, indices.size()), [=] RAJA_DEVICE (int i) {
       RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0,1), [=] RAJA_DEVICE (int i) {
#if 0          
         BoutReal div_current = alpha * gpu_Div_par_Grad_par(phi_minus_n_acc, i);
         gpu_n_ddt[i]= - gpu_bracket_par(phi_acc, n_acc, i)
                      - div_current
                                 - kappa * gpu_DZZ_par(phi_acc, i)
                                     + Dn *gpu_Delp2_par(n_acc, i) ;	

                   gpu_vort_ddt[i]= - gpu_bracket_par(phi_acc, vort_acc, i)
                                        - div_current
                                        + Dvort *gpu_Delp2_par(vort_acc, i) ;   
#endif
#if 1                                     
        int len = d2.size();                                 
        if(i < len) {
           double v2 = d2[i];
           double v3 = d3[i];
           printf("Before Assignment: d2[%d]=%f d3[%d]=%f\n",i,v2,i,v3);
        }
        d2 = d3;
        if(i < len) {
           double v2 = d2[i];
           double v3 = d3[i];
           printf("After Assignment: d2[%d]=%f d3[%d]=%f\n",i,v2,i,v3);
        }
        dd3 = dd;
        if(i==0) {
           double *b = d2.begin();
           double *e = d2.end();
           printf("begin: %p end:%p\n",b,e);
           double *dd_0 = (double*)dd.begin();
           printf("dd[0] %f %f : dd2[0] %f\n",dd_0[0],dd[0],dd2[0]);
           double *dd2Ptr = const_cast<double*>(dd2.begin());
           dd2Ptr[0] = dd[0];
           printf("dd[0] %f : dd2[0] %f dd3[0] %f\n",dd[0],dd2[0],dd3[0]);
        }

#endif
        printf("Raja Kernel[%d]\n",i);
      });
      printf("done with raja kernel dd2[0] = %f\n",dd2[0]);
    }

    return 0;
  }
};

// Define a main() function
BOUTMAIN(HW3D);
