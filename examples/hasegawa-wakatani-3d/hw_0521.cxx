/// 3D simulations of HW
/////
///// This version uses indexed operators
///// which reduce the number of loops over the domain
/////
////  GPU processing is enabled if BOUT_ENABLE_CUDA is defined
/////  Profiling markers and ranges are set if USE_NVTX is defined
/////  Baesed on Ben Duddson, Steven Glenn code, Yining Qin update 0521-2020

#include <iostream>

#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>
#include <bout/single_index_ops.hxx>

#include "RAJA/RAJA.hpp" // using RAJA lib
#include <cuda_profiler_api.h>
//#include "memoryManager.hpp"


#define BOUT_ENABLE_NVTX
#ifdef BOUT_ENABLE_NVTX
#include <cuda_profiler_api.h>
#include <nvToolsExt.h>



#define PP_CAT(a, b) PP_CAT_I(a, b)
#define PP_CAT_I(a, b) PP_CAT_II(~, a ## b)
#define PP_CAT_II(p, res) res
#define UNIQUE_NVTX_RANGE_NAME PP_CAT(range_, __COUNTER__)
#define BOUT_ENABLE_CUDA
#define BOUT_ENABLE_NVTX

static const uint32_t nvtx_colors[] = { 0xff00ff00, 0xff0000ff, 0xffffff00, 0xffff00ff,
                                        0xff00ffff, 0xffff0000, 0xffff7f7f, 0xff7f007f };
static const int num_nvtx_colors = sizeof(nvtx_colors)/sizeof(uint32_t);
class RangeManager {
public:
    RangeManager(const char* name, const int color_id) {
        int cid = color_id % num_nvtx_colors;
        nvtxEventAttributes_t eventAttrib = {0};
        eventAttrib.version = NVTX_VERSION;
        eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
        eventAttrib.colorType = NVTX_COLOR_ARGB;
        eventAttrib.color = nvtx_colors[cid];
        eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
        eventAttrib.message.ascii = name;
        id = nvtxRangeStartEx(&eventAttrib);
     }

    ~RangeManager() {nvtxRangeEnd(id);}
private:
    nvtxRangeId_t id;
};
#define PROFILE_RANGE(name, color_id) RangeManager UNIQUE_NVTX_RANGE_NAME(name, color_id)
#else
#define PROFILE_RANGE(name, color_id)
#endif // BOUT_ENABLE_NVTX

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

static inline BoutReal BOUT_DEVICE arakawa_bracket(
  const int i,                // linear mesh index
  const BoutReal* f,          // pointer to array of field values
  const BoutReal* g,          // pointer to array of field values
  const BoutReal* dx,         // pointer to 2D array of x values
  const BoutReal dz,
  const int ny,
  const int nz) {

  const auto jz = i % nz;
  const auto nyz = ny*nz;
  const auto jzmax = nz - 1;

  const auto ixp = i + nyz;
  const auto ixm = i - nyz;
  const auto izp = (jz < jzmax) ? (i + 1) : (i - jzmax);  // wrap to first element
  const auto izm = (jz > 0) ? (i - 1) : (i + jzmax);      // wrap to last element

  const auto ixpzp = izp + nyz;
  const auto ixpzm = izm + nyz;

  const auto ixmzp = izp - nyz;
  const auto ixmzm = izm - nyz;

  const BoutReal fxp = f[ixp];
  const BoutReal fxm = f[ixm];
  const BoutReal fzp = f[izp];
  const BoutReal fzm = f[izm];

  const BoutReal fpp = f[ixpzp];
  const BoutReal fpm = f[ixpzm];
  const BoutReal fmp = f[ixmzp];
  const BoutReal fmm = f[ixmzm];

  const BoutReal gxp = g[ixp];
  const BoutReal gxm = g[ixm];
  const BoutReal gzp = g[izp];
  const BoutReal gzm = g[izm];

  const BoutReal gpp = g[ixpzp];
  const BoutReal gpm = g[ixpzm];
  const BoutReal gmp = g[ixmzp];
  const BoutReal gmm = g[ixmzm];

 BoutReal Jpp =
    (fzp - fzm) * (gxp - gxm)
    - (fxp - fxm) * (gzp - gzm);


  BoutReal Jpx =
    gxp * (fpp - fpm)
    - gxm * (fmp - fmm)
    - gzp * (fpp - fmp)
    + gzm * (fpm - fmm);


  BoutReal Jxp =
    gpp * (fzp - fxp)
    - gmm * (fxm - fzm)
    - gmp * (fzp - fxm)
    + gpm * (fxp - fzm);

  const int k = i / nz;  // 2D index for dx

  return (Jpp + Jpx + Jxp) / (12. * dx[k] * dz);
}


static inline BoutReal BOUT_DEVICE delpsq(
  const int i,                // linear mesh index
  const BoutReal* f,          // pointer to Field3D data
  const BoutReal* G1,         // pointer to Field2D data
  const BoutReal* G3,         // pointer to Field2D data
  const BoutReal* g11,        // pointer to Field2D data
  const BoutReal* g13,        // pointer to Field2D data
  const BoutReal* g33,        // pointer to Field2D data
  const BoutReal* dx,         // pointer to 2D array of x values
  const BoutReal dz,
  const int ny,
  const int nz) {

  const auto jz = i % nz;
  const auto nyz = ny*nz;
  const auto jzmax = nz - 1;

  const auto ixp = i + nyz;
  const auto ixm = i - nyz;
  const auto izp = (jz < jzmax) ? (i + 1) : (i - jzmax);  // wrap to first element
  const auto izm = (jz > 0) ? (i - 1) : (i + jzmax);      // wrap to last element

  const auto ixpzp = izp + nyz;
  const auto ixpzm = izm + nyz;

  const auto ixmzp = izp - nyz;
  const auto ixmzm = izm - nyz;
  const int k = i / nz;  // index for Field2D data (has no z-dependence)

  const BoutReal fxp = f[ixp];
  const BoutReal fxm = f[ixm];
  const BoutReal fzp = f[izp];
  const BoutReal fzm = f[izm];
  const BoutReal fpp = f[ixpzp];
  const BoutReal fpm = f[ixpzm];
  const BoutReal fmp = f[ixmzp];
  const BoutReal fmm = f[ixmzm];
  const BoutReal dx2 = dx[k] * dx[k];
  const BoutReal dz2 = dz * dz;
  const BoutReal dfdx = 0.5 * (fxp - fxm) / dx[k];
  const BoutReal dfdz = 0.5 * (fzp - fzm) / dz;
  const BoutReal d2f_dx2 = (fxp + fxm - 2.0*f[i]) / dx2;
  const BoutReal d2f_dz2 = (fzp + fzm - 2.0*f[i]) / dz2;
  const BoutReal d2f_dxdz = 0.25 * (fpp - fpm - fmp + fmm) / (dx[k] * dz);

  BoutReal result = (G1[k] * dfdx) + (G3[k] * dfdz)
                  + (g11[k] * d2f_dx2) + (g33[k] * d2f_dz2)
                  + (g13[k] * d2f_dxdz);
  return result;
}


static inline BoutReal BOUT_DEVICE Div_par_Grad_par_G(
  const int i,                // linear mesh index
  const BoutReal* f,          // pointer to Field3D data
  const BoutReal* g22,        // pointer to Field2D data
  const BoutReal* dx,         // pointer to 2D array of x values
  const BoutReal* dy,
  const BoutReal dz,
  const int nx,
  const int ny,
  const int nz) {

 // const auto nxz = nx*nz;

 // const auto iyp = i + nxz;
 // const auto iym = i - nxz;
 // const int k = i / nz;  // index for Field2D data (has no z-dependence)
 // const BoutReal fyp = f[iyp];
 // const BoutReal fym = f[iym];
  //const BoutReal dfdy =  (fyp - fym) / dy[k];

//  BoutReal gradient_upper = 2.*(f.yup()[iyp] - f[i]) / (metric->dy[i] + metric->dy[iyp]);
//BoutReal gradient_upper = 2.*(f[iyp] - f[i]) / dfdy;
//BoutReal gradient_upper=0;

//BoutReal flux_upper = gradient_upper * (metric->J[i] + metric->J[iyp]) / (metric->g_22[i] + metric->g_22[iyp]);
//(metric->g_22[i] + metric->g_22[iyp])
// BoutReal mg22img22iyp=(g22[iyp]-g22[iym])/dy[k]; //+ metric->g_22[iyp]
//BoutReal flux_upper = gradient_upper*0.001 ;
BoutReal flux_upper =0.0001;

// BoutReal gradient_lower = 2.*(f[i] - f.ydown()[iym]) / (metric->dy[i] + metric->dy[iyp]);
//BoutReal gradient_lower = 2.0*(f[i]-f[iym])/dfdy;
//BoutReal gradient_lower = 0;

// BoutReal flux_lower = gradient_lower * (metric->J[i] + metric->J[iym]) / (metric->g_22[i] + metric->g_22[iym]);
//BoutReal flux_lower = gradient_lower*0.0001;
BoutReal flux_lower = 0.0001;
//return (flux_upper - flux_lower) / (metric->dy[i] * metric->J[i]);
return (flux_upper - flux_lower);

}





class HW3D : public PhysicsModel {
public:
  Field3D n, vort;  // Evolving density and vorticity
  Field3D phi;      // Electrostatic potential

  // Model parameters
  BoutReal alpha;      // Adiabaticity (~conductivity)
  BoutReal kappa;      // Density gradient drive
  BoutReal Dvort, Dn;  // Diffusion

// variables need for RHS evaluation

  int nxMesh = 0;  // Mesh x size
  int nyMesh = 0;  // Mesh y size
  int nzMesh = 0;  // Mesh z size
  BoutReal* dxMesh = nullptr;  // pointer to Field2D data
  BoutReal* dyMesh = nullptr;  // pointer to Field2D data
  BoutReal dzMesh = 0.0;
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

inline void BOUT_DEVICE calcTimeDerivatives(const int i) {
    auto dn_dt = -arakawa_bracket(i, p_phi, p_n,
                                   dxMesh, dzMesh, nyMesh, nzMesh)  ;  // ExB term
    //auto dn_dz = (p_n[i+1] - p_n[i-1]) / (2.0 * dzMesh);
    auto dphi_dz = (p_phi[i+1] - p_phi[i-1]) / (2.0 * dzMesh);
    dn_dt +=kappa* dphi_dz;

   //BoutReal div_current = alpha * Div_par_Grad_par_G(i, p_phi_minus_n,  p_g22, dxMesh,dyMesh, dzMesh,nxMesh, nyMesh, nzMesh);
   BoutReal div_current = 0;
   dn_dt += Dn*delpsq(i, p_n, p_G1, p_G3,
                                  p_g11, p_g13, p_g33,
                                  dxMesh, dzMesh, nyMesh, nzMesh) ;
  dn_dt -= div_current;
  p_dn_dt[i] = dn_dt;


   auto dvort_dt = -arakawa_bracket(i, p_phi, p_vort,
                                      dxMesh, dzMesh, nyMesh, nzMesh) ;  // ExB term

   dvort_dt += Dvort* delpsq(i, p_vort, p_G1, p_G3,
                                 p_g11, p_g13, p_g33,
                                 dxMesh, dzMesh, nyMesh, nzMesh) ;

  dvort_dt -=div_current;
  p_dvort_dt[i] = dvort_dt;
}


  // This must be public since nvcc won't allow a
    // a lambda inside a protected method.
     //
        inline void rhsImpl() {
            TRACE("HW3D::rhsImpl");
                PROFILE_RANGE("rhsImpl", 2);
 
                    // Loop over the domain, avoiding boundaries, and
                       // evaluate time derivatives at each point.
                           //
            //RAJA::cuda_exec<CUDA_BLOCK_SIZE>
/*
  using BOUTCUDANestedPolicy = RAJA::KernelPolicy<
    RAJA::statement::CudaKernel<
      RAJA::statement::Tile<1, RAJA::statement::tile_fixed<32>, RAJA::cuda_block_y_loop,
          RAJA::statement::For<1, RAJA::cuda_thread_y_direct,
              RAJA::statement::Lambda<0>
        >
      >
    > >;
*/

 
//RAJA::forall<RAJA::loop_exec>(*indices, [=] BOUT_DEVICE (int i) {         
  //
//RAJA::forall<RAJA::cuda_exec<CUDA_BLOCK_SIZE>>(*indices, [=] BOUT_DEVICE (int i) {
  RAJA::forall<EXEC_POL>(*indices, [=] BOUT_DEVICE (int i) {
              
         calcTimeDerivatives(i);
              });
         }
                 

 
  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver for vort -> phi

protected:
//std::unique_ptr<Laplacian> phiSolver; // Laplacian solver for vort -> phi
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


    ddt(n).allocate();
    ddt(vort).allocate();

    // Make sure fields have Coordinates
    // This sets the Field::fast_coords member to a Coordinate*
    // Not a long-term solution, but here until a better solution is found.
    n.fast_coords = n.getCoordinates();
    vort.fast_coords = vort.getCoordinates();
    phi.fast_coords = phi.getCoordinates();
    phi_minus_n.fast_coords = phi_minus_n.getCoordinates();



 // Set the data members used by rhsImpl() and calcTimeDerivatives()
    nxMesh = n.getNx();
    nyMesh = n.getNy();
    nzMesh = n.getNz();
    dxMesh = &n.getCoordinates()->dx(0,0);
    dyMesh = &n.getCoordinates()->dy(0,0);
    dzMesh = n.getCoordinates()->dz;
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



rhsImpl();

 /*
   BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {
//auto i = n.getRegion("RGN_NOBNDRY").begin();	
//	 for (auto &i : n.getRegion("RGN_NOBNDRY")) {
      
//	BoutReal div_current = alpha * Div_par_Grad_par(phi_minus_n, i);
	BoutReal div_current = 0;
      // Density equation
     ddt(n)[i] =-bracket(phi, n, i) - div_current - kappa * DDZ(phi, i) + Dn * Delp2(n, i);

//    ddt(n)[i] =-bracket(phi, n, i) + Dn * Delp2(n, i)- kappa * DDZ(phi, i);
      // Vorticity equation
 ddt(vort)[i] =-bracket(phi, vort, i)- div_current  + Dvort * Delp2(vort, i);
 //    ddt(vort)[i] =-bracket(phi, vort, i)+ Dvort * Delp2(vort, i);



    }
*/


    return 0;
  }
};

// Define a main() function
BOUTMAIN(HW3D);
