/// 3D simulations of HW
/////
///// This version uses indexed operators
///// which reduce the number of loops over the domain
/////
////  GPU processing is enabled if BOUT_ENABLE_CUDA is defined
/////  Profiling markers and ranges are set if USE_NVTX is defined
/////  Baesed on Ben Duddson, Steven Glenn code, Yining Qin update 0521-2020

#include <iostream>
#include <gpu_functions.hxx>

#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>
#include <bout/single_index_ops.hxx>
#include "RAJA/RAJA.hpp" // using RAJA lib
#include <cuda_profiler_api.h>

#define BOUT_ENABLE_CUDA

//  CUDA settings
#define BOUT_ENABLE_CUDA
#define BOUT_DEVICE RAJA_DEVICE
#ifdef BOUT_ENABLE_CUDA
const int CUDA_BLOCK_SIZE = 256;  // TODO: Make configurable
using EXEC_POL = RAJA::cuda_exec<CUDA_BLOCK_SIZE>;
#define BOUT_DEVICE RAJA_DEVICE
#else   // BOUT_ENABLE_CUDA not defined
using EXEC_POL = RAJA::loop_exec;
#define BOUT_DEVICE
#endif  // defined(BOUT_ENABLE_CUDA)

using IdxList = RAJA::TypedListSegment<int>;


class HW3D : public PhysicsModel {
public:
  Field3D n, vort;  // Evolving density and vorticity
  Field3D phi;      // Electrostatic potential

  // Model parameters
  BoutReal alpha;      // Adiabaticity (~conductivity)
  BoutReal kappa;      // Density gradient drive
  BoutReal Dvort, Dn;  // Diffusion

  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver for vort -> phi
// variables need for RHS evaluation

  int nxMesh = 0;  // Mesh x size
  int nyMesh = 0;  // Mesh y size
  int nzMesh = 0;  // Mesh z size
  BoutReal* dxMesh = nullptr;  // pointer to Field2D data
  BoutReal* dyMesh = nullptr;  // pointer to Field2D data
  BoutReal dzMesh = 0.0;
  BoutReal* JMesh = nullptr;
  IdxList* indices = nullptr; 

  // Raw pointers for field data needed for RHS.
  // Unless noted, data are owned by Field3D objects.

  BoutReal* p_n = nullptr;
  BoutReal* p_phi = nullptr;
  BoutReal* p_phi_minus_n= nullptr;
  BoutReal* p_vort = nullptr;
  BoutReal* p_dn_dt = nullptr;
  BoutReal* p_dvort_dt = nullptr;
  BoutReal* p_G1 = nullptr;         // Field2D
  BoutReal* p_G3 = nullptr;         // Field2D
  BoutReal* p_g11 = nullptr;         // Field2D
  BoutReal* p_g13 = nullptr;         // Field2D
  BoutReal* p_g33 = nullptr;         // Field2D
  BoutReal* p_g22 = nullptr;         // Field2D

  BoutReal* phi_minus_n_acc_yup = nullptr;
  BoutReal* phi_minus_n_acc_ydown= nullptr;   


inline void BOUT_DEVICE gpu_calcTimeDerivatives(const int i) {
 
 /* ------------------- original HW3d loops
//----CPU functions starts ---------------------

BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {

      BoutReal div_current = alpha * Div_par_Grad_par(phi_minus_n_acc, i);


      ddt(n)[i] = -bracket(phi_acc, n_acc, i)
                  - div_current
                  - kappa * DDZ(phi_acc, i)
                  + Dn * Delp2(n_acc, i);


      ddt(vort)[i] = -bracket(phi_acc, vort_acc, i)
                     - div_current
                     + Dvort * Delp2(vort_acc, i);
    }
//-------- CPU ends ---------------------------------
*/

   	 
	auto dphi_dz = (p_phi[i+1] - p_phi[i-1]) / (2.0 * dzMesh);
	BoutReal div_current = alpha * gpu_Div_par_Grad_par_G(i, p_phi_minus_n,  p_g22, dxMesh,dyMesh,
                                          dzMesh,JMesh,nxMesh, nyMesh, nzMesh,phi_minus_n_acc_yup,
						phi_minus_n_acc_ydown);
	auto dn_dt = -gpu_arakawa_bracket(i, p_phi, p_n,
                                   dxMesh, dzMesh, nyMesh, nzMesh)  ;  
   
    	dn_dt +=kappa* dphi_dz;

  	dn_dt -= div_current;
	p_dn_dt[i] = dn_dt;

   	auto dvort_dt = -gpu_arakawa_bracket(i, p_phi, p_vort,
                                      dxMesh, dzMesh, nyMesh, nzMesh) ;  

   	dvort_dt += Dvort* gpu_delpsq(i, p_vort, p_G1, p_G3,
                                 p_g11, p_g13, p_g33,
                                 dxMesh, dzMesh, nyMesh, nzMesh) ;

  	dvort_dt -=div_current;
  	p_dvort_dt[i] = dvort_dt;

}


  inline void gpu_rhsImpl() {
  
 	 RAJA::forall<EXEC_POL>(*indices, [=] BOUT_DEVICE (int i) {
              
        	gpu_calcTimeDerivatives(i);
              	});
         }
                 

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
    phi = 0.; // Starting phi



    std::list<int> tmpList;
    for(auto&& regidx : n.getRegion("RGN_NOBNDRY").getIndices()) {
        tmpList.push_back(regidx.ind);
    }
    indices = new IdxList(tmpList);
 
    
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
                        

  // Copy the data members used by GPU device function gpu_rhsImpl()
    nxMesh = n.getNx();
    nyMesh = n.getNy();
    nzMesh = n.getNz();
    dxMesh = &n.getCoordinates()->dx(0,0);
    dyMesh = &n.getCoordinates()->dy(0,0);
    dzMesh = n.getCoordinates()->dz;
    JMesh = &n.getCoordinates()->J(0,0);
    p_n = const_cast<BoutReal*>(n(0,0));
    p_vort = const_cast<BoutReal*>(vort(0,0));
    p_dn_dt= const_cast<BoutReal*>(n.timeDeriv()->operator()(0,0));
    p_dvort_dt = const_cast<BoutReal*>(vort.timeDeriv()->operator()(0,0));
    p_phi = const_cast<BoutReal*>(phi(0,0));
    p_phi_minus_n = const_cast<BoutReal*>(phi_minus_n(0,0));
    p_G1 = &n.getCoordinates()->G1(0,0);
    p_G3 = &n.getCoordinates()->G3(0,0);
    p_g11 = &n.getCoordinates()->g11(0,0);
    p_g13 = &n.getCoordinates()->g13(0,0);
    p_g33 = &n.getCoordinates()->g33(0,0);
    p_g22 = &n.getCoordinates()->g22(0,0);  

	const Field3D &yup = *phi_minus_n_acc.yup;
  	const Field3D &ydown = *phi_minus_n_acc.ydown;
   	phi_minus_n_acc_yup = const_cast<BoutReal*>(yup(0,0));
   	phi_minus_n_acc_ydown = const_cast<BoutReal*>(ydown(0,0));

//----GPU function entry
   gpu_rhsImpl();

/*
//----CPU functions starts ---------------------

BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {

      BoutReal div_current = alpha * Div_par_Grad_par(phi_minus_n_acc, i);


      ddt(n)[i] = -bracket(phi_acc, n_acc, i)
                  - div_current
                  - kappa * DDZ(phi_acc, i)
                  + Dn * Delp2(n_acc, i);

     
      ddt(vort)[i] = -bracket(phi_acc, vort_acc, i)
                     - div_current
                     + Dvort * Delp2(vort_acc, i);
    }
//-------- CPU ends ---------------------------------
*/
    return 0;
  }
};

// Define a main() function
BOUTMAIN(HW3D);
