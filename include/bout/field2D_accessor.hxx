// GPU version Field-Accessor, updated by Dr. Yining Qin, Oct.27, 2020

#pragma once
#ifndef FIELD2D_ACCESSOR_H__
#define FIELD2D_ACCESSOR_H__

#include "../bout_types.hxx"
#include "../field3d.hxx"
#include "../field.hxx"
#include "../field2d.hxx"


#if defined(BOUT_USE_CUDA) && defined(__CUDACC__)
#define BOUT_HOST_DEVICE __host__ __device__
#define BOUT_HOST __host__
#define BOUT_DEVICE __device__
#else
#define BOUT_HOST_DEVICE
#define BOUT_HOST
#define BOUT_DEVICE
#endif

template<CELL_LOC location = CELL_CENTRE>
struct Field2DAccessor {
  explicit Field2DAccessor(Field2D &f) {
    ASSERT0(f.getLocation() == location);
    coords = f.getCoordinates();
    
    ASSERT0(f.isAllocated());
  
    data = &f(0,0);
    //----- Field 2d data -> array for GPU    
    f_data = &f(0,0);

    // Field size for GPU 
    f_nx = f.getNx();
    f_ny = f.getNy();
    f_nz = f.getNz();
    
   //------ Field2D data to array for GPU
    f2d_dx= &f.getCoordinates()->dx(0,0);
    f2d_dy= &f.getCoordinates()->dy(0,0);
    f2d_dz= f.getCoordinates()->dz;
    f2d_J= &f.getCoordinates()->J(0,0);
    
    f2d_G1= &f.getCoordinates()->G1(0,0);
    f2d_G3= &f.getCoordinates()->G3(0,0);
    f2d_g11= &f.getCoordinates()->g11(0,0);
    f2d_g12= &f.getCoordinates()->g12(0,0);
    f2d_g13= &f.getCoordinates()->g13(0,0);
    f2d_g22= &f.getCoordinates()->g22(0,0);
    f2d_g23= &f.getCoordinates()->g23(0,0);
    f2d_g33= &f.getCoordinates()->g33(0,0);
 //---------------------------------------

  if (f.hasParallelSlices()) {
      yup = &f.yup();
      ydown = &f.ydown();
    
     // set field region index for GPU
     // auto indices = f.getRegion("RGN_NOBNDRY").getIndices(); //set index of region
     // Ind3D *ob_i = &(indices[0]);
     //---------------------


	}
  }

    BoutReal& BOUT_HOST_DEVICE operator[](const Ind2D &d) {
    return data[d.ind];
  }
  const BoutReal& BOUT_HOST_DEVICE operator[](const Ind2D &d) const {
    return data[d.ind];
  }
  

 Field2D   *yup {nullptr};

 Field2D   *ydown {nullptr};

 Coordinates*  coords {nullptr};

 BoutReal* data;

// --------------  define for GPU
 BoutReal* f2d_dx = nullptr;
 BoutReal* f2d_dy = nullptr;  // pointer to Field2D data
 BoutReal f2d_dz = 0.0;       // dz of Field2D

 BoutReal* f2d_J = nullptr;  // pointer to Field2D data
 
 BoutReal* f2d_G1 = nullptr;
 BoutReal* f2d_G3 = nullptr;
 BoutReal* f2d_g11 = nullptr;
 BoutReal* f2d_g12 = nullptr;
 BoutReal* f2d_g13 = nullptr;
 BoutReal* f2d_g22 = nullptr;
 BoutReal* f2d_g23 = nullptr;
 BoutReal* f2d_g33 = nullptr;

 BoutReal* f_yup = nullptr;
 BoutReal* f_ydown = nullptr;
 BoutReal* f_data = nullptr;

 BoutReal* f_ddt = nullptr;

 int* ind_2D_index = nullptr;
  
 int f_nx =0;
 int f_ny =0;
 int f_nz =0;
//--------------------------------

};

#endif 
