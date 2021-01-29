
//--  GPU version Update by Dr. Yining Qin, Oct27, 2020

#pragma once
#ifndef SINGLE_INDEX_OPS_H
#define SINGLE_INDEX_OPS_H

#include "field2D_accessor.hxx"
#include "field_accessor.hxx"

#if defined(BOUT_USE_CUDA) && defined(__CUDACC__)
#define BOUT_HOST_DEVICE __host__ __device__
#define BOUT_HOST __host__
#define BOUT_DEVICE __device__
#else
#define BOUT_HOST_DEVICE
#define BOUT_HOST
#define BOUT_DEVICE
#endif


//--  RAJA CUDA settings--------------------------------------------------------start
#ifdef BOUT_USE_CUDA
const int CUDA_BLOCK_SIZE = 256;  // TODO: Make configurable
using EXEC_POL = RAJA::cuda_exec<CUDA_BLOCK_SIZE>;
#else   // BOUT_ENABLE_CUDA not defined
using EXEC_POL = RAJA::loop_exec;
#endif  // defined(BOUT_ENABLE_CUDA)
////-----------CUDA settings------------------------------------------------------end

template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal* DDT( const FieldAccessor<location> &f){

return f.f_ddt;

} 
template< CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal* FIELD_DATA(const FieldAccessor<location> &f) {

// use raw pointer to field3D data
   return f.f_data;
}

template< CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal* FIELD2D_DATA(const Field2DAccessor<location> &f) {

// use raw pointer to fieldi2D data
   return f.f_data;
}

template< CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal FIELD2D_3DINDEX_DATA(const FieldAccessor<location> &f3d, const Field2DAccessor<location> &f2d, const int i ) {

int  nz = f3d.f_nz;
int  ind_2d = i / nz;  // index for Field2D data (has no z-dependence)
    
return f2d.f_data[ind_2d];
   
 }




// Ind3D: i.zp():
BOUT_HOST_DEVICE inline int i_zp(const int id, const int nz){
 int jz = id % nz;
 int jzmax = nz -1;
 return (jz < jzmax) ? (id + 1) : (id - jzmax);
}
// Ind3D: i.zm():
 BOUT_HOST_DEVICE inline int i_zm(const int id, const int nz){
 int jz = id % nz;
 int jzmax = nz -1;
 return (jz > 0) ? (id - 1) : (id + jzmax);
}

// Ind3D: i.ym(): 
BOUT_HOST_DEVICE inline int i_yp(const int id, const int nz){
return id + nz;
}

// Ind3D: i.yp();
BOUT_HOST_DEVICE inline int i_ym(const int id, const int nz){
return id - nz;
}


// Ind3D: i.xp();
BOUT_HOST_DEVICE inline int i_xp(const int id, const int ny, const int nz){
return id + ny*nz;
}


// Ind3D: i.xm();
BOUT_HOST_DEVICE inline int i_xm(const int id, const int ny, const int nz){
return id - ny*nz;
}

template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal  bracket_2d3D_g(const Field2DAccessor<location> &f2d, const FieldAccessor<location> &f3d, const int i) {

BoutReal result =  b0xGrad_dot_Grad_g(f2d,f3d,i);
return result;
}
template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal  bracket_g(const FieldAccessor<location> &f, const FieldAccessor<location> &g, const int i) {

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


template<IND_TYPE N, CELL_LOC location>
inline BoutReal bracket(const FieldAccessor<location> &f, const FieldAccessor<location> &g, const SpecificInd<N> &ind) {
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

template<CELL_LOC location>

BOUT_HOST_DEVICE inline BoutReal DDX_g(const Field2DAccessor<location> &f, const int i) {
   
   BoutReal* dx = f.f2d_dx; 
   BoutReal* f_a = f.f_data;
   //int nx = f.f_nx;
   int ny = f.f_ny;
   int nz = f.f_nz;
   int ind_2d = i / nz;  // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)

   int ixp = i_xp(i,ny,nz);
   int ixm = i_xm(i,ny,nz);

  
    return 0.5 * (f_a[ixp] - f_a[ixm]) / dx[ind_2d];

}


template<CELL_LOC location>

BOUT_HOST_DEVICE inline BoutReal DDX_g(const FieldAccessor<location> &f, const int i) {
   
   BoutReal* dx = f.f2d_dx; 
   BoutReal* f_a = f.f_data;
   int nx = f.f_nx;
   int ny = f.f_ny;
   int nz = f.f_nz;
   int ind_2d = i / nz;  // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)

   int ixp = i_xp(i,ny,nz);
   int ixm = i_xm(i,ny,nz);

  
    return 0.5 * (f_a[ixp] - f_a[ixm]) / dx[ind_2d];

}

template<IND_TYPE N, CELL_LOC location>
BoutReal DDX(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.xp()] - f[ind.xm()]) / (2. * f.coords->dx[ind]);
}

template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal DDY_g(const Field2DAccessor<location> &f, const int i) {
   
 
// use raw pointer to field data 
 // BoutReal* f_a = f.f_data;
  //auto dx = f.f2d_dx;
    BoutReal* dy = f.f2d_dy;
  //auto dz = f.f2d_dz;
  //BoutReal* J = f.f2d_J;
  //BoutReal* g22 = f.f2d_g22;
  int  nz = f.f_nz;
 
  // Use the raw pointers to yup/ydown fields. 
  BoutReal* yup = f.f_yup;
  BoutReal* ydown = f.f_ydown;
  
  int iyp = i_yp(i,nz); // auto iyp = i.yp();
  int iym = i_ym(i,nz); // auto iym = i.ym():
  int ind_2d = i / nz;  // index for Field2D data (has no z-dependence)

  BoutReal result = 0.5*(yup[iyp] - ydown[iym]) / dy[ind_2d]; // metric->dy[i] is dy
  
    //return ( (*f.yup)[ind.yp()] - (*f.ydown)[ind.ym()]) / (2. * f.coords->dy[ind]); // original operator
   return result;
}


template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal DDY_g(const FieldAccessor<location> &f, const int i) {
   
 
// use raw pointer to field data 
 // BoutReal* f_a = f.f_data;
  //auto dx = f.f2d_dx;
    BoutReal* dy = f.f2d_dy;
  //auto dz = f.f2d_dz;
  //BoutReal* J = f.f2d_J;
  //BoutReal* g22 = f.f2d_g22;
  int  nz = f.f_nz;
 
  // Use the raw pointers to yup/ydown fields. 
  BoutReal* yup = f.f_yup;
  BoutReal* ydown = f.f_ydown;
  
  int iyp = i_yp(i,nz); // auto iyp = i.yp();
  int iym = i_ym(i,nz); // auto iym = i.ym():
  int ind_2d = i / nz;  // index for Field2D data (has no z-dependence)

  BoutReal result = 0.5*(yup[iyp] - ydown[iym]) / dy[ind_2d]; // metric->dy[i] is dy
  
    //return ( (*f.yup)[ind.yp()] - (*f.ydown)[ind.ym()]) / (2. * f.coords->dy[ind]); // original operator
   return result;
}


template<IND_TYPE N, CELL_LOC location>
inline BoutReal DDY(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return ( (*f.yup)[ind.yp()] - (*f.ydown)[ind.ym()]) / (2. * f.coords->dy[ind]);
}

template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal DDZ_g(const FieldAccessor<location> &f, const int i) {
   BoutReal dz = f.f2d_dz;
   int nz = f.f_nz;
   BoutReal* f_a = f.f_data;
   int izp = i_zp(i,nz);
   int izm = i_zm(i,nz);
  
   return (f_a[izp] - f_a[izm]) / (2. * dz);
}


template<IND_TYPE N, CELL_LOC location>
inline BoutReal  DDZ(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
}


template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal Delp2_g(const FieldAccessor<location> &f, const int i) {

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


template<IND_TYPE N, CELL_LOC location>
inline BoutReal Delp2(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
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


template< CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal Div_par_Grad_par_g(const FieldAccessor<location> &f, const int i) {

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


template<IND_TYPE N, CELL_LOC location>
inline BoutReal Div_par_Grad_par(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
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
  
  BoutReal gradient_upper = 2.*(yup[iyp] - f[i]) / (dy + metric->dy[iyp]);
  BoutReal flux_upper = gradient_upper * (J + metric->J[iyp]) / (g_22 + metric->g_22[iyp]);
  
  BoutReal gradient_lower = 2.*(f[i] - ydown[iym]) / (dy + metric->dy[iyp]);
  BoutReal flux_lower = gradient_lower * (J + metric->J[iym]) / (g_22 + metric->g_22[iym]);

  return (flux_upper - flux_lower) / (dy * J);
}

template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal b0xGrad_dot_Grad_g(const FieldAccessor<location> &f3d, const Field2DAccessor<location> &f2d, const int i) {


return 0.0001;

}


template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal VDDX_g(const Field2DAccessor<location> &f2d, const FieldAccessor<location> &f3d, const int i) {
 
 return 0.0001;
}


template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal VDDY_g(const Field2DAccessor<location> &f2d, const FieldAccessor<location> &f3d, const int i) {
  

return 0.0001;

}

template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal VDDZ_g(const Field2DAccessor<location> &f2d, const FieldAccessor<location> &f3d, const int i) {
 

return 0.0001;
 
}


template<CELL_LOC location>
BOUT_HOST_DEVICE  BoutReal SQ_g(const FieldAccessor<location> &f3d,const Field2DAccessor<location> &f2d,const int i) {

 int  nz = f3d.f_nz;
 int  ind_2d = i / nz;  // index for Field2D data (has no z-dependence)

 //BoutReal* f3d_a = f3d.f_data;
 BoutReal* f2d_a = f2d.f_data;
 return f2d_a[ind_2d] * f2d_a[ind_2d];

}


template<CELL_LOC location>
BOUT_HOST_DEVICE  BoutReal b0xGrad_dot_Grad_3D2D_g(const FieldAccessor<location> &f3d,const Field2DAccessor<location> &f2d,const int i) {

  
  int  nz = f3d.f_nz;
  int  ind_2d = i / nz;  // index for Field2D data (has no z-dependence)

  BoutReal* J = f3d.f2d_J;   
   BoutReal* g_23 = f2d.f2d_g23;
   BoutReal* g_12 = f2d.f2d_g12;
   BoutReal* g_22 = f2d.f2d_g22;

  // Calculate phi derivatives
   BoutReal dpdx = DDX_g(f3d,  i);
   BoutReal dpdy = DDY_g(f3d,  i);
   BoutReal dpdz = DDZ_g(f3d,  i);
   
 
  // Calculate advection velocity
   BoutReal vx = g_22[ind_2d] * dpdz -  g_23[ind_2d] * dpdy;
   BoutReal vy = g_23[ind_2d] * dpdx -  g_12[ind_2d] * dpdz;
  
  // VDDX,VDDY
   BoutReal* f3d_a = f3d.f_data;
   BoutReal* dx = f2d.f2d_dx;
   BoutReal* dy = f2d.f2d_dy;
   
   BoutReal vddx = f3d_a[i] / dx[ind_2d];  // VDDX()
   BoutReal vddy = f3d_a[i] / dy[ind_2d];  // VDDY()

 // result 
   BoutReal p1 = vddx + vddy;
   BoutReal p2 = J[ind_2d] * sqrt(g_22[ind_2d]);

   BoutReal result = p1;
   result /= p2;

   return result;



}


template<CELL_LOC location>
BOUT_HOST_DEVICE  BoutReal b0xGrad_dot_Grad_g(const Field2DAccessor<location> &f2d, const FieldAccessor<location> &f3d, const int i) {
  
int  nz = f3d.f_nz;
int  ind_2d = i / nz;  // index for Field2D data (has no z-dependence)


  // Calculate phi derivatives
   BoutReal dpdx = DDX_g(f2d,  ind_2d);
   BoutReal dpdy = DDY_g(f2d,  ind_2d);
   
   BoutReal* J = f2d.f2d_J;   
   BoutReal* g_23 = f2d.f2d_g23;
   BoutReal* g_12 = f2d.f2d_g12;
   BoutReal* g_22 = f2d.f2d_g22;

  // Calculate advection velocity
   BoutReal vx = g_23[ind_2d] * dpdy;
   BoutReal vy = g_23[ind_2d] * dpdx;
   BoutReal vz = g_12[ind_2d] * dpdy - g_22[ind_2d] * dpdx;
// VDDX,VDDY,VDDZ
   BoutReal* f2d_a = f2d.f_data;
   BoutReal* dx = f3d.f2d_dx;
   BoutReal* dy = f3d.f2d_dy;
   BoutReal dz = f3d.f2d_dz;
   
   BoutReal vddx = f2d_a[ind_2d] / dx[i];  // VDDX()
   BoutReal vddy = f2d_a[ind_2d] / dy[i];  // VDDY()
   BoutReal vddz = f2d_a[ind_2d] / dz;          // VDDZ()

 // result 
   BoutReal p1 = vddx + vddy + vddz;
   BoutReal p2 = J[ind_2d] * sqrt(g_22[ind_2d]);

   BoutReal result = p1;
   result /= p2;
   if ((result < -1) || (result >1)){result = 0.0;}
   (isinf(result) == 1) ? (result = 0.0): result = result;
   (isnan(result) == 1) ? (result = 0.0): result = result;
   
   return result*0.001;



}

template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal D2DY2_g(const FieldAccessor<location> &f3d, const int i) {

int  nz = f3d.f_nz;
int  ind_2d = i / nz;  // index for Field2D data (has no z-dependence)

BoutReal* f3d_a = f3d.f_data;
BoutReal* dy = f3d.f2d_dy;
BoutReal result = f3d_a[i] / (dy[ind_2d]* dy[ind_2d]);
return result;

}

template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal Grad_parP_g(const FieldAccessor<location> &f3d, const int i) {

BoutReal result =  Grad_par_g(f3d,i);

   //if (nonlinear) {
   //   result -= bracket(interp_to(Psi, loc), f, bm_mag) * B0;

     // if (include_rmp) {
      //  result -= bracket(interp_to(rmp_Psi, loc), f, bm_mag) * B0;
     // }
   // }


return result;
}



template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal Grad_par_g(const FieldAccessor<location> &f3d, const int i) {

BoutReal* g22 = f3d.f2d_g22;
int  nz = f3d.f_nz;
int  ind_2d = i / nz;  // index for Field2D data (has no z-dependence)

BoutReal ddy = DDY_g(f3d,i);
BoutReal result = ddy / sqrt(g22[ind_2d]); 
return result;
//return ::DDY(var, outloc, method) / sqrt(g_22); // original operator
}


template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal* Grad_g(const FieldAccessor<location> &f3d, const int i) {


BoutReal ddx = DDX_g(f3d,  i);
BoutReal ddy = DDY_g(f3d,  i);
BoutReal ddz = DDZ_g(f3d,  i);

BoutReal result[2] ;

result[0] = ddx;
result[1] = ddy;
result[2] = ddz;

return result;

}


template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal filter_g(const FieldAccessor<location> &f3d, const int N0,const int i) {

int ncx = f3d.f_nx;
int ncy = f3d.f_ny;
int ncz = f3d.f_nz;

 /// Convenience functions for converting to (x, y, z)
 //  int x() const { return (ind / nz) / ny; }
 //    int y() const { return (ind / nz) % ny; }
 //      int z() const { return (ind % nz); }
 //
 //int z()
int  ind_z = i % ncz;  // index for Z data

BoutReal* f3d_a = f3d.f_data;
//BoutReal re = f3d_a[i];

      //for (int jz = 0; jz <= ncz / 2; jz++) {
        //if (jz != N0) {
          // Zero this component
          //           f[jz] = 0.0;
          //                   }
          //                         }
     if  (( ind_z >= 0) && (ind_z <= ncz/2) && (ind_z != N0)  )
         {f3d_a[i] = 0.0;}

return f3d_a[i];


}


template<CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal lowPass_g(const FieldAccessor<location> &f3d, const int zmax, const bool keep_zonal, const int i) {

int ncx = f3d.f_nx;
int ncy = f3d.f_ny;
int ncz = f3d.f_nz;
int  ind_z = i % ncz;  // index for Z data

BoutReal* f3d_a = f3d.f_data;
BoutReal re =  f3d_a[i];

  if (((zmax >= ncz / 2) || (zmax < 0)) && keep_zonal) {
 
         re =  f3d_a[i];;
   
  } else {
     if  ((ind_z >= 0) && (ind_z <= ncz/2) )      
	re = 0.0;
	if (!keep_zonal){f3d_a[0]=0.0;}
	
    }
return re;

}





#endif // SINGLE_INDEX_OPS_H
