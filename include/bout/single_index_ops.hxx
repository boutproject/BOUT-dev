
//--  GPU version Update by Dr. Yining Qin, Oct27, 2020

#pragma once
#ifndef SINGLE_INDEX_OPS_H
#define SINGLE_INDEX_OPS_H

#include "field_accessor.hxx"
<<<<<<< HEAD

#if defined(BOUT_USE_CUDA) && defined(__CUDACC__)
#define BOUT_HOST_DEVICE __host__ __device__
#define BOUT_HOST __host__
#define BOUT_DEVICE __device__
#else
#define BOUT_HOST_DEVICE
#define BOUT_HOST
#define BOUT_DEVICE
#endif


=======

#define BOUT_HOST_DEVICE __host__ __device__
#define RAJA_HOST_DEVICE __host__ __device__

#define BOUT_HAS_UMPIRE
#define BOUT_USE_CUDA

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

template<CELL_LOC location>
RAJA_DEVICE BoutReal* DDT( const FieldAccessor<location> &f){

return f.f_ddt;

} 

// Ind3D: i.zp():
BOUT_DEVICE int i_zp(const int id, const int nz){
 int jz = id % nz;
 int jzmax = nz -1;
 return (jz < jzmax) ? (id + 1) : (id - jzmax);
}
// Ind3D: i.zm():
BOUT_DEVICE int i_zm(const int id, const int nz){
 int jz = id % nz;
 int jzmax = nz -1;
 return (jz > 0) ? (id - 1) : (id + jzmax);
}

// Ind3D: i.ym(): 
BOUT_DEVICE int i_yp(const int id, const int nz){
return id + nz;
}

// Ind3D: i.yp();
BOUT_DEVICE int i_ym(const int id, const int nz){
return id - nz;
}


// Ind3D: i.xp();
BOUT_DEVICE int i_xp(const int id, const int ny, const int nz){
return id + ny*nz;
}


// Ind3D: i.xm();
BOUT_DEVICE int i_xm(const int id, const int ny, const int nz){
return id - ny*nz;
}

template<CELL_LOC location>
BOUT_DEVICE BoutReal  bracket_g(const FieldAccessor<location> &f, const FieldAccessor<location> &g, const int i) {

  BoutReal* dx = g.f2d_dx;
  BoutReal dz = g.f2d_dz;
  BoutReal* f_a =f.f_data;
  BoutReal* g_a =g.f_data;

  int nx = g.f_nx;
  int ny = g.f_ny;
  int nz = g.f_nz;
  int ind_2d = i / nz;  // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)
 

  // Offset indices
//  auto ixp = ind.xp();
//  auto ixm = ind.xm();
//  auto izp = ind.zp();
//  auto izm = ind.zm();

//  auto izpxp = izp.xp();
//  auto izpxm = izp.xm();
//  auto izmxp = izm.xp();
//  auto izmxm = izm.xm();

   int ixp = i_xp(i,ny,nz);
   int ixm = i_xm(i,ny,nz);
   int izp = i_zp(i,nz);
   int izm = i_zm(i,nz);

   int izpxp = i_xp(izp,ny,nz);
   int izpxm = i_xm(izp,ny,nz);
   int izmxp = i_xp(izm,ny,nz);
   int izmxm = i_xm(izm,ny,nz);

    
  // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
  BoutReal Jpp = ((f_a[izp] - f_a[izm]) * (g_a[ixp] - g_a[ixm])
                  - (f_a[ixp] - f_a[ixm]) * (g_a[izp] - g_a[izm]));
  
  // J+x
  BoutReal Jpx = (g_a[ixp] * (f_a[izpxp] - f_a[izmxp]) - g_a[ixm] * (f_a[izpxm] - f_a[izmxm])
                  - g_a[izp] * (f_a[izpxp] - f_a[izpxm]) + g_a[izm] * (f_a[izmxp] - f_a[izmxm]));

  // Jx+
  BoutReal Jxp = (g_a[izpxp] * (f_a[izp] - f_a[ixp]) - g_a[izmxm] * (f_a[ixm] - f_a[izm])
                  - g_a[izpxm] * (f_a[izp] - f_a[ixm]) + g_a[izmxp] * (f_a[ixp] - f_a[izm]));
  
  return (Jpp + Jpx + Jxp) / (12 * dx[ind_2d] * dz);
}


>>>>>>> next-outerloop-GPU-umpire
template<IND_TYPE N, CELL_LOC location>
BoutReal bracket(const FieldAccessor<location> &f, const FieldAccessor<location> &g, const SpecificInd<N> &ind) {
  Coordinates *metric = g.coords;

  // Offset indices
  auto ixp = ind.xp();
  auto ixm = ind.xm();
  auto izp = ind.zp();
  auto izm = ind.zm();

  auto izpxp = izp.xp();
  auto izpxm = izp.xm();
  auto izmxp = izm.xp();
  auto izmxm = izm.xm();
  
  // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
  BoutReal Jpp = ((f[izp] - f[izm]) * (g[ixp] - g[ixm])
                  - (f[ixp] - f[ixm]) * (g[izp] - g[izm]));
  
  // J+x
  BoutReal Jpx = (g[ixp] * (f[izpxp] - f[izmxp]) - g[ixm] * (f[izpxm] - f[izmxm])
                  - g[izp] * (f[izpxp] - f[izpxm]) + g[izm] * (f[izmxp] - f[izmxm]));

  // Jx+
  BoutReal Jxp = (g[izpxp] * (f[izp] - f[ixp]) - g[izmxm] * (f[ixm] - f[izm])
                  - g[izpxm] * (f[izp] - f[ixm]) + g[izmxp] * (f[ixp] - f[izm]));
  
  return (Jpp + Jpx + Jxp) / (12 * metric->dx[ind] * metric->dz);
}

template<IND_TYPE N, CELL_LOC location>
<<<<<<< HEAD
BoutReal DDX(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.xp()] - f[ind.xm()]) / (2. * f.coords->dx[ind]);
}

template<IND_TYPE N, CELL_LOC location>
BoutReal DDY(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return ( (*f.yup)[ind.yp()] - (*f.ydown)[ind.ym()]) / (2. * f.coords->dy[ind]);
}

template<IND_TYPE N, CELL_LOC location>
BoutReal DDZ(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
}

=======
BoutReal DDX_g(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
   
  int nz = f.f_nz;
  auto dx = f.f2d_dx;

  int ind_2d = ind.ind / nz;  // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)

  return (f[ind.xp()] - f[ind.xm()]) / (2. * dx[ind_2d]);
}

template<IND_TYPE N, CELL_LOC location>
BoutReal DDX(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.xp()] - f[ind.xm()]) / (2. * f.coords->dx[ind]);
}

template<IND_TYPE N, CELL_LOC location>
BoutReal DDY(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return ( (*f.yup)[ind.yp()] - (*f.ydown)[ind.ym()]) / (2. * f.coords->dy[ind]);
}

template<CELL_LOC location>
BOUT_DEVICE BoutReal DDZ_g(const FieldAccessor<location> &f, const int i) {
   BoutReal dz = f.f2d_dz;
   int nz = f.f_nz;
   BoutReal* f_a = f.f_data;
   int izp = i_zp(i,nz);
   int izm = i_zm(i,nz);
  
   return (f_a[izp] - f_a[izm]) / (2. * dz);
}


template<IND_TYPE N, CELL_LOC location>
BoutReal  DDZ(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
}


template<CELL_LOC location>
RAJA_DEVICE BoutReal Delp2_g(const FieldAccessor<location> &f, const int i) {

  BoutReal* dx = f.f2d_dx;
  BoutReal dz = f.f2d_dz;
  BoutReal* f_a = f.f_data;

  BoutReal* G1 = f.f2d_G1;
  BoutReal* G3 = f.f2d_G3;
  BoutReal* g11 = f.f2d_g11;
  BoutReal* g13 = f.f2d_g13;
  BoutReal* g33 = f.f2d_g33;

  int nx = f.f_nx;
  int ny = f.f_ny;
  int nz = f.f_nz;
  int ind_2d = i / nz;  // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)
 

//  auto izm = i.zm();
//  auto izp = i.zp();
//  auto ixm = i.xm();
//  auto ixp = i.xp();

   int ixp = i_xp(i,ny,nz);
   int ixm = i_xm(i,ny,nz);
   int izp = i_zp(i,nz);
   int izm = i_zm(i,nz);

   int izpxp = i_xp(izp,ny,nz);
   int izpxm = i_xm(izp,ny,nz);
   int izmxp = i_xp(izm,ny,nz);
   int izmxm = i_xm(izm,ny,nz);

  BoutReal DDX = G1[ind_2d] * (f_a[ixp] - f_a[ixm]) / (2.0 * dx[ind_2d]); //DDX
  BoutReal DDZ =  G3[ind_2d] * (f_a[izp] - f_a[izm]) / (2.0 * dz);  // DDZ
  BoutReal D2DX2 =  g11[ind_2d] * (f_a[ixp] - 2.*f_a[i] + f_a[ixm]) / (dx[ind_2d]*dx[ind_2d]);  // D2DX2
  BoutReal D2DZ2 =  g33[ind_2d] * (f_a[izp] - 2.0*f_a[i] + f_a[izm]) /(dz*dz); // D2DZ2 
  BoutReal D2DXDZ = 2 * g13[ind_2d] * ((f_a[izpxp] - f_a[izpxm]) - (f_a[izmxp] - f_a[izmxm])) / (4. * dz * dx[ind_2d]) ; // D2DXDZ
  
  BoutReal res = DDX + DDZ + D2DX2 + D2DZ2 + D2DXDZ;
  return res;
}


>>>>>>> next-outerloop-GPU-umpire
template<IND_TYPE N, CELL_LOC location>
BoutReal Delp2(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
  Coordinates *metric = f.coords;

  // Index offsets
  auto izm = i.zm();
  auto izp = i.zp();
  auto ixm = i.xm();
  auto ixp = i.xp();
 
  return metric->G1[i] * (f[ixp] - f[ixm]) / (2.*metric->dx[i])  // DDX
    + metric->G3[i] * (f[izp] - f[izm]) / (2.*metric->dz)  // DDZ
    + metric->g11[i] * (f[ixp] - 2.*f[i] + f[ixm]) / SQ(metric->dx[i])  // D2DX2
    + metric->g33[i] * (f[izp] - 2.*f[i] + f[izm]) / SQ(metric->dz) // D2DZ2
    + 2 * metric->g13[i] * ((f[izp.xp()] - f[izp.xm()]) - (f[izm.xp()] - f[izm.xm()])) / (4. * metric->dz * metric->dx[i]) // D2DXDZ
    ;
}

<<<<<<< HEAD
#if 0
template<IND_TYPE N, CELL_LOC location>
BOUT_HOST_DEVICE BoutReal Delp2_gt(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
  Coordinates *metric = f.coords;

  // Index offsets
  auto izm = i.zm();
  auto izp = i.zp();
  auto ixm = i.xm();
  auto ixp = i.xp();

  auto g1= metric->G1[i];
  //return metric->G1[i] * (f[ixp] - f[ixm]) / (2.*metric->dx[i])  // DDX
  return 0.5*(f[ixm] - f[ixp]) ;

}
#else

template<IND_TYPE N, CELL_LOC location>
BOUT_HOST_DEVICE BoutReal Delp2_gt(const FieldAccessorLite<location> &f,const SpecificInd<N> &i) {

  // Index offsets
  auto izm = i.zm();

  auto yup = f.yup;

  //metrics
  auto G1 = f.G1;
  auto g_11 = f.g_11;

  //return f[izm] * G1[i.ind] * g_11[i.ind] * yup[i.ind] ;
  return yup[i.ind];

}
#endif
=======
template< CELL_LOC location>
RAJA_DEVICE BoutReal Div_par_Grad_par_g(const FieldAccessor<location> &f, const int i) {

// use raw pointer to field data 
  BoutReal* f_a = f.f_data;
  //auto dx = f.f2d_dx;
    BoutReal* dy = f.f2d_dy;
  //auto dz = f.f2d_dz;
  BoutReal* J = f.f2d_J;
  BoutReal* g22 = f.f2d_g22;
  int  nz = f.f_nz;
 
  // Use the raw pointers to yup/ydown fields. 
  BoutReal* yup = f.f_yup;
  BoutReal* ydown = f.f_ydown;
  
  // Index offsets
  //auto iyp = i.yp();
  //auto iym = i.ym();
  //const int  ind_2d = i.ind / nz;  // index for Field2D data (has no z-dependence)
  //int  iyp_2d = iyp.ind / nz; // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)
  //int  iym_2d = iym.ind / nz; // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)
 
  int iyp = i_yp(i,nz); // auto iyp = i.yp();
  int iym = i_ym(i,nz); // auto iym = i.ym():
  int  ind_2d = i / nz;  // index for Field2D data (has no z-dependence)
  int  iyp_2d = iyp/ nz;
  int  iym_2d = iym/nz;

  BoutReal gradient_upper = 2.*(yup[iyp] - f_a[i]) / (dy[ind_2d] + dy[iyp_2d]); // metric->dy[i] is dy
  BoutReal flux_upper = gradient_upper * (J[ind_2d] + J[iyp_2d]) / (g22[ind_2d] + g22[iyp_2d]); //metric->g22[i]
  BoutReal gradient_lower = 2.*(f_a[i] - ydown[iym]) / (dy[ind_2d] + dy[iyp_2d]);
  BoutReal flux_lower = gradient_lower * (J[ind_2d] + J[iym_2d]) / (g22[ind_2d] + g22[iym_2d]);

  BoutReal output =  (flux_upper - flux_lower) / (dy[ind_2d] * J[ind_2d]);
  return output;

}

>>>>>>> next-outerloop-GPU-umpire

template<IND_TYPE N, CELL_LOC location>
BoutReal Div_par_Grad_par(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
  Coordinates *metric = f.coords;
  // Index offsets
  auto iyp = i.yp();
  auto iym = i.ym();

  // Use the raw pointers to yup/ydown fields. These must have been set before calling
  const Field3D &yup = *f.yup;
  const Field3D &ydown = *f.ydown;

  // Fetch values used more than once
  BoutReal dy = metric->dy[i];
  BoutReal J = metric->J[i];
  BoutReal g_22 = metric->g_22[i];
<<<<<<< HEAD
  
  BoutReal gradient_upper = 2.*(yup[iyp] - f[i]) / (dy + metric->dy[iyp]);
  BoutReal flux_upper = gradient_upper * (J + metric->J[iyp]) / (g_22 + metric->g_22[iyp]);
  
=======
  
  BoutReal gradient_upper = 2.*(yup[iyp] - f[i]) / (dy + metric->dy[iyp]);
  BoutReal flux_upper = gradient_upper * (J + metric->J[iyp]) / (g_22 + metric->g_22[iyp]);
  
>>>>>>> next-outerloop-GPU-umpire
  BoutReal gradient_lower = 2.*(f[i] - ydown[iym]) / (dy + metric->dy[iyp]);
  BoutReal flux_lower = gradient_lower * (J + metric->J[iym]) / (g_22 + metric->g_22[iym]);

  return (flux_upper - flux_lower) / (dy * J);
}




#endif // SINGLE_INDEX_OPS_H
