
 ////  GPU processing is enabled if BOUT_ENABLE_CUDA is defined
 /////  GPU funcition for BOUT++ 
 /////  Yining Qin update 0730-2020

#include <iostream>
#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>
#include "RAJA/RAJA.hpp" // using RAJA lib
#include <cuda_profiler_api.h>
#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"


#define BOUT_ENABLE_CUDA
#define UMPIRE_ENABLE_CUDA

#include "bout/field_accessor.hxx"

 
namespace {
template <typename G>
using itr = G*;
template <typename G>
using const_itr = const G*;
}
template <typename G>
  struct UMData {
    int len;    ///< Size of the array
    G *data;    ///< Array of data

    UMData(int size) : len(size) {
#ifdef BOUT_ENABLE_CUDA
      const std::string MEM_RESOURCE_NAME{"UM"};
#else
      const std::string MEM_RESOURCE_NAME{"HOST"};
#endif
      auto allocator = umpire::ResourceManager::getInstance().getAllocator(MEM_RESOURCE_NAME);
      data = static_cast<G*>(allocator.allocate(size*sizeof(G)));
    }
    ~UMData() {
      auto allocator = umpire::ResourceManager::getInstance().getAllocator(data);
      allocator.deallocate(data);
    }
  };


void allocate_and_deallocate(const std::string& resource)
{
  constexpr std::size_t SIZE = 1024;

  auto& rm = umpire::ResourceManager::getInstance();

  umpire::Allocator allocator = rm.getAllocator(resource);

  double* data =
      static_cast<double*>(allocator.allocate(SIZE * sizeof(double)));

  std::cout << "Allocated " << (SIZE * sizeof(double)) << " bytes using the "
            << allocator.getName() << " allocator...";

  allocator.deallocate(data);

  std::cout << " deallocated." << std::endl;
}


void copy_data(double* source_data, std::size_t size,
               const std::string& destination)
{
  auto& rm = umpire::ResourceManager::getInstance();
  auto dest_allocator = rm.getAllocator(destination);

  double* dest_data =
      static_cast<double*>(dest_allocator.allocate(size * sizeof(double)));

  rm.copy(dest_data, source_data);

  std::cout << "Copied source data (" << source_data << ") to destination "
            << destination << " (" << dest_data << ")" << std::endl;

 // dest_allocator.deallocate(dest_data);
}





//----prepare varaibles for GPU----------------------------------------------- Start
  __managed__  int nxMesh = 0;  // Mesh x size
  __managed__  int nyMesh = 0;  // Mesh y size
  __managed__ int nzMesh = 0;  // Mesh z size
   
  __managed__ double test1 = 0.999;

  __managed__ BoutReal* dxMesh = nullptr;  // pointer to Field2D data
  __managed__ BoutReal* dyMesh = nullptr;  // pointer to Field2D data
  __managed__ BoutReal dzMesh = 0.0;
  __managed__ BoutReal* JMesh = nullptr;

  constexpr std::size_t SIZE = 1;
    auto& rm = umpire::ResourceManager::getInstance();
    auto dest_allocator = rm.getAllocator("UM");
 
    double* dest_data = static_cast<double*>(dest_allocator.allocate(SIZE * sizeof(double)));
  
// Raw pointers to copy data to GPU.
  // Unless noted, data are owned by Field3D objects.
  
 __managed__ BoutReal* p_n = nullptr;
 __managed__ BoutReal* p_phi = nullptr;
 __managed__ BoutReal* p_phi_minus_n= nullptr;
 __managed__ BoutReal* p_vort = nullptr;
 __managed__ BoutReal* p_dn_dt = nullptr;
 __managed__ BoutReal* p_dvort_dt = nullptr;
 __managed__ BoutReal* p_G1 = nullptr;         // Field2D
 __managed__ BoutReal* p_G3 = nullptr;         // Field2D
 __managed__ BoutReal* p_g11 = nullptr;         // Field2D
 __managed__ BoutReal* p_g13 = nullptr;         // Field2D
 __managed__ BoutReal* p_g33 = nullptr;         // Field2D
 __managed__ BoutReal* p_g22 = nullptr;         // Field2D
 __managed__ BoutReal* phi_minus_n_acc_yup = nullptr;
 __managed__ BoutReal* phi_minus_n_acc_ydown= nullptr;   


 __managed__ BoutReal r_gpu_DZZ_par = 1.987654321432143214321;


static void RAJA_data_copy(const Field3D n, const Field3D vort, const Field3D phi,const Field3D phi_minus_n, const FieldAccessor<> &phi_acc, const FieldAccessor<> &phi_minus_n_acc, int size) {
  
    nxMesh = n.getNx();
    nyMesh = n.getNy();
    nzMesh = n.getNz();
    dxMesh = &n.getCoordinates()->dx(0,0);
    dyMesh = &n.getCoordinates()->dy(0,0);
    dzMesh = n.getCoordinates()->dz;
    JMesh = &n.getCoordinates()->J(0,0);
    p_n = const_cast<BoutReal*>(n(0,0));
    p_vort = const_cast<BoutReal*>(vort(0,0));
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

 // allocate_and_deallocate("UM");
 
 constexpr std::size_t SIZE = 1;

  auto& rm = umpire::ResourceManager::getInstance();

  auto allocator = rm.getAllocator("UM");

  double* data =
      static_cast<double*>(allocator.allocate(SIZE * sizeof(double)));
      data[0] = 0.00;

   rm.copy(dest_data, data);
   test1 = dest_data[0];
//#if defined(UMPIRE_ENABLE_CUDA)
 // copy_data(data, SIZE, "DEVICE");
 // copy_data(data, SIZE, "UM");
 // copy_data(data, SIZE, "PINNED");
//#endif 

}


static  BoutReal RAJA_DEVICE gpu_delpsq(
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

static BoutReal RAJA_DEVICE gpu_arakawa_bracket(
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


static  BoutReal RAJA_DEVICE gpu_Div_par_Grad_par_G(
  const int i,                // linear mesh index
  const BoutReal* f,          // pointer to Field3D data
  const BoutReal* g22,        // pointer to Field2D data
  const BoutReal* dx,         // pointer to 2D array of x values
  const BoutReal* dy,
  const BoutReal dz,
  const BoutReal* J,
  const int nx,
  const int ny,
  const int nz,
  const BoutReal* yup,
  const BoutReal* ydown) {

  const auto iyp = i + nz; // auto iyp = i.yp();
  const auto iym = i - nz; // auto iym = i.ym():
  const int  k = i / nz;  // index for Field2D data (has no z-dependence)
  const int  kyp = iyp/ nz;
  const int  kym = iym/nz;


BoutReal gradient_upper = 2.*(yup[kyp] - f[k]) / (dy[k]+dy[kyp]);
//BoutReal gradient_upper = 2.*(yup[kyp] - f[i]) / (dy[k]+dy[kyp]);
//BoutReal gradient_upper = 2.*(yup[iyp] - f[i]) / (dy[k]+dy[kyp]);

BoutReal flux_upper = gradient_upper*(J[k]+J[kyp])/ (g22[k]+g22[kyp]) ;
//BoutReal flux_upper = gradient_upper*(J[k]+J[kyp])/ (g22[k]+g22[kyp]) ;
//BoutReal flux_upper = gradient_upper*(J[i]+J[kyp])/ (g22[k]+g22[iyp]) ;

BoutReal gradient_lower = 2.0*(f[k]-ydown[kym])/(dy[k]+ dy[kyp]);
//BoutReal gradient_lower = 2.0*(f[i]-ydown[kym])/(dy[k]+ dy[kyp]);
//BoutReal gradient_lower = 2.0*(f[i]-ydown[kym])/(dy[k]+ dy[kyp]);


BoutReal flux_lower = gradient_lower*(J[k]+J[kym])/(g22[k]+ g22[kym]);
//BoutReal flux_lower = gradient_lower*(J[k]+J[kym])/(g22[k]+ g22[kym]);
 //BoutReal flux_lower = gradient_lower*(J[k]+J[kym])/(g22[k]+ g22[iym]);

BoutReal output =(flux_upper - flux_lower) /(dy[k]*J[k]);
//BoutReal output =(flux_upper - flux_lower) /(dy[k] * J[k]) ;
return output;
}



static  BoutReal RAJA_DEVICE gpu_DZZ(
const int i,
BoutReal* f,
const BoutReal dz
){
auto dphi_dz = (f[i+1] - f[i-1]) / (2.0 * dz);
return dphi_dz;

}




template<CELL_LOC location>
static BoutReal RAJA_DEVICE gpu_Delp2_par(const FieldAccessor<location> &f, const int i) {

auto dp2=gpu_delpsq(i, p_n, p_G1, p_G3, p_g11, p_g13, p_g33, dxMesh, dzMesh, nyMesh, nzMesh) ;

return dp2;
}

template<CELL_LOC location>
static BoutReal RAJA_DEVICE gpu_bracket_par(const FieldAccessor<location> &f, const FieldAccessor<location> &g, const int i) {
	//  Coordinates *metric = g.coords;
 	
	auto gbp = gpu_arakawa_bracket(i, p_phi, p_n, dxMesh, dzMesh, nyMesh, nzMesh);
	return gbp;

}

template< CELL_LOC location>
static BoutReal RAJA_DEVICE  gpu_DZZ_par(const FieldAccessor<location> &f, const int i) {
  	// Index offsets
//    	auto izm = i.zm();
//     	auto izp = i.zp();
//	return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
//	return (f[izp] - f[izm] / (2. * f.coords->dz);
	auto  ddz = gpu_DZZ(i,p_phi,dzMesh);
//	return ddz;
//	r_gpu_DZZ_par = (p_phi[i+1] - p_phi[i-1]) / (2.0 * dzMesh);	
     //  printf("GPU ddz = %f\n", ddz);
//	cudaDeviceSynchronize();
 //      printf("GPU nxMesh = %d\n", nxMesh);
//       r_gpu_DZZ_par = 0.01 * nxMesh;
       return ddz;
//	return dphi_dz;
//	return 0.0005;
}


template< CELL_LOC location>
static BoutReal RAJA_DEVICE  gpu_Div_par_Grad_par(const FieldAccessor<location> &f, const int i) {  
//	Coordinates *metric = f.coords;

  // Index offsets
  // auto iyp = i.yp();
  // auto iym = i.ym();
	  //Use the raw pointers to yup/ydown fields. These must have been set before calling
  //   const Field3D &yup = *f.yup;
  //   const Field3D &ydown = *f.ydown; 
//auto dx = dxMesh;
//  BoutReal dy = metric->dy[i];
//BoutReal dz = metric->dz;
  //BoutReal g_22 = metric->g_22(0,0);

auto div_p =  gpu_Div_par_Grad_par_G(i, p_phi_minus_n,  p_g22, dxMesh,dyMesh,
                                           dzMesh,JMesh,nxMesh, nyMesh, nzMesh,phi_minus_n_acc_yup, phi_minus_n_acc_ydown);


  return div_p;
//return 0.001;
     }




