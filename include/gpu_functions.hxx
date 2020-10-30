
 ////  GPU processing is enabled if BOUT_ENABLE_CUDA is defined
 /////  GPU funcition for BOUT++ 
 /////  Yining Qin and Holger Jones update from 0730-2020

#include <iostream>
#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>
#include "RAJA/RAJA.hpp" // using RAJA lib
#include <cuda_profiler_api.h>
#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"

#include "bout/field_accessor.hxx"

#define BOUT_ENABLE_CUDA
#define UMPIRE_ENABLE_CUDA

#define BOUT_HAS_UMPIRE
#define BOUT_USE_CUDA

#define ddt(name) gpu_ddt_helper(#name, name )


//--  RAJA CUDA settings--------------------------------------------------------start
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
////-----------CUDA settings------------------------------------------------------end

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
//----prepare varaibles for GPU----------------------------------------------- Start
// int region_size =0;

UMData<int> test(3);

constexpr std::size_t SIZE = 10;
auto& rm = umpire::ResourceManager::getInstance();
umpire::Allocator allocator = rm.getAllocator("UM");

double* g_data = static_cast<double*>(allocator.allocate(SIZE * sizeof(double)));




__managed__  int nxMesh = 0;  // Mesh x size
__managed__  int nyMesh = 0;  // Mesh y size
__managed__ int nzMesh = 0;  // Mesh z size
__managed__ BoutReal* dxMesh = nullptr;  // pointer to Field2D data
__managed__ BoutReal* dyMesh = nullptr;  // pointer to Field2D data
__managed__ BoutReal dzMesh = 0.0;
__managed__ BoutReal* JMesh = nullptr;


// Raw pointers to copy data to GPU. Unless noted, data are owned by Field3D objects.

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

__managed__  BoutReal* gpu_n_ddt; // copy ddt(n) to __device__
__managed__  BoutReal* gpu_vort_ddt; // copy ddt(vort) to __device__  


//static BoutReal* gpu_ddt( char* name, fiedl3d &f){
//	gpu_n_ddt= const_cast<BoutReal*>(f.timeDeriv()->operator()(0,0));
//	return gpu_n_ddt;

//}
static BoutReal* RAJA_DEVICE gpu_ddt_helper(const char* p_name, Field3D &f ){
auto ddt = gpu_n_ddt;
//std::cout<<"Name is  "<< p_name <<std::endl;
//if (p_name == "n"){std::cout<<"I am n:     "<< p_name <<std::endl;}
//if (p_name == "vort"){std::cout<<"I am vort:     "<< p_name <<std::endl;}
if (p_name == "n"){ddt =gpu_n_ddt;}
if (p_name == "vort"){ddt = gpu_vort_ddt;}
return ddt;

}

// copy data from host to device
static void RAJA_data_copy( Field3D &n, Field3D &vort,  Field3D &phi, Field3D &phi_minus_n,  const FieldAccessor<> &phi_acc, const FieldAccessor<> &phi_minus_n_acc) {
  
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
 

 //   RAJA::synchronize<RAJA::cuda_synchronize>();	
    gpu_n_ddt= const_cast<BoutReal*>(n.timeDeriv()->operator()(0,0)); // copy ddt(n) to __device__
    gpu_vort_ddt = const_cast<BoutReal*>(vort.timeDeriv()->operator()(0,0)); // copy ddt(vort) to __device__  

    //auto region = n.getRegion("RGN_NOBNDRY"); // Region object
    //auto indices = region.getIndices();   // A std::vector of Ind3D objects 
   // region_size = indices.size();


//g_data[0]= 0.0099;

  for (std::size_t i = 0; i < SIZE; i++) {
    g_data[i] = 0.0;
  }

g_data[0]= 0.9988776;

}

static  BoutReal RAJA_HOST_DEVICE g_test( const BoutReal t){
	 return t;
}


static  BoutReal RAJA_HOST_DEVICE gpu_delpsq(
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

static BoutReal RAJA_HOST_DEVICE gpu_arakawa_bracket(
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


static  BoutReal RAJA_HOST_DEVICE gpu_Div_par_Grad_par_G(
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
  const int  kyp = iyp / nz;
  const int  kym = iym / nz;


	BoutReal gradient_upper = 2.*(yup[iyp] - f[i]) / (dy[k] + dy[kyp]);
	BoutReal flux_upper = gradient_upper*(J[k]+J[kyp]) / (g22[k] + g22[kyp]) ;
	BoutReal gradient_lower = 2.0*(f[i] - ydown[iym])/(dy[k] + dy[kyp]);
	BoutReal flux_lower = gradient_lower*(J[k] + J[kym])/(g22[k]+ g22[kym]);
	BoutReal output =(flux_upper - flux_lower) /(dy[k]*J[k]);
	return output;

}



static  BoutReal RAJA_HOST_DEVICE gpu_DZZ(
const int i,
BoutReal* f,
const BoutReal dz
){
	auto dphi_dz = (f[i+1] - f[i-1]) / (2.0 * dz);
	return dphi_dz;

}


template<CELL_LOC location>
static BoutReal RAJA_HOST_DEVICE Delp2(const FieldAccessor<location> &f, const int i) {
        Coordinates *metric = f.coords;
        BoutReal dz = metric->dz;
	auto dp2=gpu_delpsq(i, p_n, p_G1, p_G3, p_g11, p_g13, p_g33, dxMesh, dz, nyMesh, nzMesh) ;
	return dp2;
}

template<CELL_LOC location>
static BoutReal RAJA_HOST_DEVICE bracket(const FieldAccessor<location> &f, const FieldAccessor<location> &g, const int i) {
	Coordinates *metric = g.coords;
	BoutReal dz = metric->dz;
	auto gbp = gpu_arakawa_bracket(i, p_phi, p_n, dxMesh, dz, nyMesh, nzMesh);
	return gbp;

}

/*
template<IND_TYPE N, CELL_LOC location>
BoutReal DDZ(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
}
*/

template< CELL_LOC location>
static BoutReal RAJA_HOST_DEVICE  DDZ(const FieldAccessor<location> &f, const int i) {

	Coordinates *metric = f.coords;
	auto dzz = metric->dz;
       auto ddz_o = 0.5* (p_phi[i+1] - p_phi[i-1]) / dzz;	
       return ddz_o;
}

template<IND_TYPE N, CELL_LOC location>
static BoutReal RAJA_HOST_DEVICE t_DDZ(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
	Coordinates *metric = f.coords;
	auto dzz = metric->dz;
	int  dz = 1;
	auto nz = nzMesh;
	dz = dz <= nz ? dz :dz % nz;
        // ind.zp();
	int izp_ind = ((ind.ind + dz) % nz < dz ? ind.ind - nz + dz : ind.ind + dz);
	SpecificInd<IND_TYPE::IND_3D> izp = ind;
	izp.ind = izp_ind;
	// ind.zm()
	int izm_ind = ind.ind % nz < dz ? ind.ind + nz - dz : ind.ind - dz;  // izm = ind.zm();
	SpecificInd<IND_TYPE::IND_3D> izm = ind;
	izm.ind = izm_ind;
	auto ddz_o = 0.5*(f[izp]-f[izm]) / dzz; // return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
	return ddz_o;
}



template< CELL_LOC location>
static BoutReal RAJA_HOST_DEVICE  Div_par_Grad_par(const FieldAccessor<location> &f, const int i) {  
	Coordinates *metric = f.coords;
//	BoutReal t =g_data[0];
	BoutReal dz = metric->dz;
	auto div_p =  gpu_Div_par_Grad_par_G(i, p_phi_minus_n,  p_g22, dxMesh,dyMesh,
                                           dz,JMesh,nxMesh, nyMesh, nzMesh,phi_minus_n_acc_yup, phi_minus_n_acc_ydown);
  return div_p;
     }


template<IND_TYPE N, CELL_LOC location>
static BoutReal RAJA_HOST_DEVICE t_Div(const FieldAccessor<location> &f, const SpecificInd<N> &i) {

	Coordinates *metric = f.coords;
	auto dzz = metric->dz;
	int  ind_dy = 1;
	auto nz = nzMesh;
	auto ny = nyMesh;
  	// Index offsets
  	// auto iyp = i.yp();
	int ind_iyp = i.ind + (ind_dy * nz);
	int ind_iym = i.ind + (-ind_dy * nz);
	SpecificInd<IND_TYPE::IND_3D> iyp = i;
	iyp.ind = ind_iyp;
	SpecificInd<IND_TYPE::IND_3D> iym = i;
	iym.ind = ind_iym;

	// Use the raw pointers to yup/ydown fields. These must have been set before calling
	const Field3D &yup = *f.yup;
  	const Field3D &ydown = *f.ydown;
 

	  // Fetch values used more than once
	SpecificInd<IND_TYPE::IND_2D> i_2D ;
	i_2D.ind = i.ind/nz;
	i_2D.ny= ny;
	i_2D.nz = 1;
	//BoutReal tempArr= metric->dy.data.ptr.data[i_2D.ind];
//	ArrayData<BoutReal>* arrPtrTemp = &(metric->dy.data.ptr);
//	BoutReal dy = tempArr.data[i_2D.ind];	
	//BoutReal dy = metric->dy[i_2D];
	  //     BoutReal J = metric->J[i];
	  //       BoutReal g_22 = metric->g_22[i];	

//   phi_minus_n_acc_yup = const_cast<BoutReal*>(yup(0,0));
//   phi_minus_n_acc_ydown = const_cast<BoutReal*>(ydown(0,0));
	
	


	//BoutReal gradient_upper = 2.*(yup[iyp] - f[i]) / (dy + metric->dy[iyp]);
//	BoutReal gradient_upper = 2.*(phi_minus_n_acc_yup[iyp.ind] - f[i]) / (dyMesh[k] + dyMesh[k]);
	
//	BoutReal flux_upper = gradient_upper * (J + metric->J[iyp]) / (g_22 + metric->g_22[iyp]);
	
//	BoutReal flux_upper = gradient_upper * (JMesh[k] + JMesh[kyp]) / (p_g22[k] + p_g22[kyp]);


//  	BoutReal gradient_lower = 2.*(f[i] - ydown[iym]) / (dy + metric->dy[iyp]);

  //	BoutReal gradient_lower = 2.*(f[i] - phi_minus_n_acc_ydown[iym.ind]) / (dyMesh[k] + dyMesh[kyp]);
  	

//	BoutReal flux_lower = gradient_lower * (J + metric->J[iym]) / (g_22 + metric->g_22[iym]);

//	BoutReal flux_lower = gradient_lower * (JMesh[k] + JMesh[kym]) / (p_g22[k] + p_g22[kym]);
  
//	BoutReal r =  (flux_upper - flux_lower) / (dyMesh[k] * JMesh[k]);
	BoutReal r = 0.001;
	return r;

//	return 0.001;
}


