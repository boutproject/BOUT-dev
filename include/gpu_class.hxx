
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


//----prepare varaibles for GPU----------------------------------------------- Start
  int nxMesh = 0;  // Mesh x size
  int nyMesh = 0;  // Mesh y size
  int nzMesh = 0;  // Mesh z size
  BoutReal* dxMesh = nullptr;  // pointer to Field2D data
  BoutReal* dyMesh = nullptr;  // pointer to Field2D data
  BoutReal dzMesh = 0.0;
  BoutReal* JMesh = nullptr;


  
// Raw pointers to copy data to GPU.
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




void data_copy(const Field3D n, const Field3D vort, const FieldAccessor<> &f,const FieldAccessor<> &g,int size) {
    nxMesh = n.getNx();
    nyMesh = n.getNy();
    nzMesh = n.getNz();
    dxMesh = &n.getCoordinates()->dx(0,0);
    dyMesh = &n.getCoordinates()->dy(0,0);
    dzMesh = n.getCoordinates()->dz;
    JMesh = &n.getCoordinates()->J(0,0);
    p_n = const_cast<BoutReal*>(n(0,0));
    p_vort = const_cast<BoutReal*>(vort(0,0));
    p_G1 = &n.getCoordinates()->G1(0,0);
    p_G3 = &n.getCoordinates()->G3(0,0);
    p_g11 = &n.getCoordinates()->g11(0,0);
    p_g13 = &n.getCoordinates()->g13(0,0);
    p_g33 = &n.getCoordinates()->g33(0,0);
    p_g22 = &n.getCoordinates()->g22(0,0);

}



template<IND_TYPE N, CELL_LOC location>
BoutReal RAJA_DEVICE g_gpu_Delp2_par(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
return 0.0001;

}

template<IND_TYPE N, CELL_LOC location>
static BoutReal RAJA_DEVICE gpu_bracket_par(const FieldAccessor<location> &f, const FieldAccessor<location> &g, const SpecificInd<N> &ind) {
        //  Coordinates *metric = g.coords;
         return 0.001;
 }

template<IND_TYPE N, CELL_LOC location>
static BoutReal RAJA_DEVICE  gpu_DZZ_par(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
        // Index offsets
        // //      auto izm = i.zm();
        // //      auto izp = i.zp();
        // //      return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
        // //      return (f[izp] - f[izm] / (2. * f.coords->dz);
        //         return 0.001 / (2. * f.coords->dz);
              return 0.0005;
  }
        
tatic BoutReal RAJA_DEVICE  gpu_Div_par_Grad_par(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
        Coordinates *metric = f.coords;

  // Index offsets
  //  auto iyp = i.yp();
  //  auto iym = i.ym();
  //   Use the raw pointers to yup/ydown fields. These must have been set before calling
    const Field3D &yup = *f.yup;
    const Field3D &ydown = *f.ydown;
    //auto dx = dxMesh;
    //  BoutReal dy = metric->dy[i];
    BoutReal dz = metric->dz;
  //                           //BoutReal g_22 = metric->g_22(0,0);
                                  return 0.0009;
  
                                            }
  








