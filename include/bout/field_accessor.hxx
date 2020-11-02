<<<<<<< HEAD
=======
// GPU version Field-Accessor, updated by Dr. Yining Qin, Oct.27, 2020

>>>>>>> next-outerloop-GPU-umpire
#pragma once
#ifndef FIELD_ACCESSOR_H__
#define FIELD_ACCESSOR_H__

#include "../bout_types.hxx"
#include "../field3d.hxx"
<<<<<<< HEAD
=======
#include "../field.hxx"
#include "../field2d.hxx"

>>>>>>> next-outerloop-GPU-umpire

template<CELL_LOC location = CELL_CENTRE>
struct FieldAccessor {
  explicit FieldAccessor(Field3D &f) {
    ASSERT0(f.getLocation() == location);

<<<<<<< HEAD
    // Get raw pointers
    coords = f.getCoordinates();

    ASSERT0(f.isAllocated());
    
    // Underlying data array
    data = &f(0,0,0);
    
    if (f.hasParallelSlices()) {
      yup = &f.yup();
      ydown = &f.ydown();
    }
  }

  BoutReal& operator[](const Ind3D &d) {
    return data[d.ind];
  }
  const BoutReal& operator[](const Ind3D &d) const {
    return data[d.ind];
  }
  
  /// These are for fast (inner loop) access to yup/ydown
  /// Used in single_index_ops.hxx functions
  Field3D *yup {nullptr}, *ydown {nullptr};

  Coordinates* coords {nullptr};

  /// Internal data array
  BoutReal* data;
};

template<CELL_LOC location = CELL_CENTRE>
struct FieldAccessorLite {
  explicit FieldAccessorLite(Field3D &f) {
    ASSERT0(f.getLocation() == location);

    // Get raw pointers
    Coordinates* coords = f.getCoordinates();
    G1=&coords->G1(0,0,0);
    g_11=&coords->g_11(0,0,0);

    ASSERT0(f.isAllocated());
    
    // Underlying data array
    data = &f(0,0,0);
    
    if (f.hasParallelSlices()) {
      auto test = f.yup()(0,0);
      printf("type test = %s\n",typeid(test).name());
      yup = f.yup()(0,0);
      yup[3] = 3.33;
      ydown = f.ydown()(0,0);
    }
  }

  BoutReal& operator[](const Ind3D &d) {
    return data[d.ind];
  }
  const BoutReal& operator[](const Ind3D &d) const {
    return data[d.ind];
  }
  
  /// These are for fast (inner loop) access to yup/ydown
  /// Used in single_index_ops.hxx functions
  //Field3D *yup {nullptr}, *ydown {nullptr};
  BoutReal *yup {nullptr}, *ydown {nullptr};

  //Coordinates* coords {nullptr};

  /// Internal data array
  BoutReal* data;
  BoutReal* G1;
  BoutReal* g_11;
};

#endif // FIELD_ACCESSOR_H__
=======

    coords = f.getCoordinates();
    
    ASSERT0(f.isAllocated());
  
    data = &f(0,0,0);

    //----- Field 3d data -> array for GPU    
    f_data = f(0,0);

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
    f2d_g13= &f.getCoordinates()->g13(0,0);
    f2d_g22= &f.getCoordinates()->g22(0,0);
    f2d_g33= &f.getCoordinates()->g33(0,0);
 //---------------------------------------

  if (f.hasParallelSlices()) {
      yup = &f.yup();
      ydown = &f.ydown();
     
    //----- Field3D data -> array for GPU
     f_yup =  static_cast<BoutReal*>(f.yup()(0,0));
     f_ydown =  static_cast<BoutReal*>(f.ydown()(0,0));
     // ddt() array data for GPU
     f_ddt =  static_cast<BoutReal*>(f.timeDeriv()->operator()(0,0));

    // set field region index for GPU
     // auto indices = f.getRegion("RGN_NOBNDRY").getIndices(); //set index of region
     // Ind3D *ob_i = &(indices[0]);
     //---------------------


	}
  }

    BoutReal& __host__ __device__ operator[](const Ind3D &d) {
    return data[d.ind];
  }
  const BoutReal& __host__ __device__ operator[](const Ind3D &d) const {
    return data[d.ind];
  }
  

 Field3D   *yup {nullptr};

 Field3D   *ydown {nullptr};

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
 BoutReal* f2d_g13 = nullptr;
 BoutReal* f2d_g22 = nullptr;
 BoutReal* f2d_g33 = nullptr;

 BoutReal* f_yup = nullptr;
 BoutReal* f_ydown = nullptr;
 BoutReal* f_data = nullptr;

 BoutReal* f_ddt = nullptr;

 int* ind_3D_index = nullptr;
  
 int f_nx =0;
 int f_ny =0;
 int f_nz =0;
//--------------------------------

};

#endif 
>>>>>>> next-outerloop-GPU-umpire
