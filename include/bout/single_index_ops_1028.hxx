
//--  GPU version Update by Dr. Yining Qin, Oct27, 2020

#pragma once
#ifndef SINGLE_INDEX_OPS_H
#define SINGLE_INDEX_OPS_H

#include "field_accessor.hxx"

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


RAJA_HOST_DEVICE BoutReal* DDT( Field3D &f){

return static_cast<BoutReal*>(f.timeDeriv()->operator()(0,0));

} 


template<IND_TYPE N, CELL_LOC location>
BOUT_HOST_DEVICE BoutReal  bracket_g(const FieldAccessor<location> &f, const FieldAccessor<location> &g, const SpecificInd<N> &ind) {


  // Offset indices
  auto ixp = ind.xp();
  auto ixm = ind.xm();
  auto izp = ind.zp();
  auto izm = ind.zm();

  auto izpxp = izp.xp();
  auto izpxm = izp.xm();
  auto izmxp = izm.xp();
  auto izmxm = izm.xm();
 
  int nz = g.f_nz;
  int ind_2d = ind.ind / nz;  // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)

  auto dx = g.f2d_dx;
  auto dz = g.f2d_dz;
   
  // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
  BoutReal Jpp = ((f[izp] - f[izm]) * (g[ixp] - g[ixm])
                  - (f[ixp] - f[ixm]) * (g[izp] - g[izm]));
  
  // J+x
  BoutReal Jpx = (g[ixp] * (f[izpxp] - f[izmxp]) - g[ixm] * (f[izpxm] - f[izmxm])
                  - g[izp] * (f[izpxp] - f[izpxm]) + g[izm] * (f[izmxp] - f[izmxm]));

  // Jx+
  BoutReal Jxp = (g[izpxp] * (f[izp] - f[ixp]) - g[izmxm] * (f[ixm] - f[izm])
                  - g[izpxm] * (f[izp] - f[ixm]) + g[izmxp] * (f[ixp] - f[izm]));
  
  return (Jpp + Jpx + Jxp) / (12 * dx[ind_2d] * dz);
}


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
BoutReal DDX(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.xp()] - f[ind.xm()]) / (2. * f.coords->dx[ind]);
}

template<IND_TYPE N, CELL_LOC location>
BoutReal DDY(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return ( (*f.yup)[ind.yp()] - (*f.ydown)[ind.ym()]) / (2. * f.coords->dy[ind]);
}

template<IND_TYPE N, CELL_LOC location>
BOUT_HOST_DEVICE BoutReal DDZ_g(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
   auto dz = f.f2d_dz;
   return (f[ind.zp()] - f[ind.zm()]) / (2. * dz);
}

template<IND_TYPE N, CELL_LOC location>
BoutReal  DDZ(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
}


template<IND_TYPE N, CELL_LOC location>
RAJA_HOST_DEVICE BoutReal Delp2_g(const FieldAccessor<location> &f, const SpecificInd<N> &i) {

  auto izm = i.zm();
  auto izp = i.zp();
  auto ixm = i.xm();
  auto ixp = i.xp();
  
  int nz = f.f_nz;
  int ind_2d = i.ind / nz;  // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)

  auto dx = f.f2d_dx;
  auto dz = f.f2d_dz;
  auto G1 = f.f2d_G1;
  auto G3 = f.f2d_G3;
  auto g11 = f.f2d_g11;
  auto g13 = f.f2d_g13;
  auto g33 = f.f2d_g33;

  BoutReal DDX = G1[ind_2d] * (f[ixp] - f[ixm]) / (2.0 * dx[ind_2d]); //DDX
  BoutReal DDZ =  G3[ind_2d] * (f[izp] - f[izm]) / (2.0 * dz);  // DDZ
  BoutReal D2DX2 =  g11[ind_2d] * (f[ixp] - 2.*f[i] + f[ixm]) / (dx[ind_2d]*dx[ind_2d]);  // D2DX2
  BoutReal D2DZ2 =  g33[ind_2d] * (f[izp] - 2.0*f[i] + f[izm]) /(dz*dz); // D2DZ2 
  BoutReal D2DXDZ = 2 * g13[ind_2d] * ((f[izp.xp()] - f[izp.xm()]) - (f[izm.xp()] - f[izm.xm()])) / (4. * dz * dx[ind_2d]) ; // D2DXDZ
  
  BoutReal res = DDX + DDZ + D2DX2 + D2DZ2 + D2DXDZ;
  return res;
}


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


template<IND_TYPE N, CELL_LOC location>
RAJA_HOST_DEVICE BoutReal Div_par_Grad_par_g(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
  
  //auto f_a = f.f_data;
  //auto dx = f.f2d_dx;
  auto dy = f.f2d_dy;
  //auto dz = f.f2d_dz;
  auto J = f.f2d_J;
  
  auto g22 = f.f2d_g22;
  

  int  nz = f.f_nz;
  int  ind_2d = i.ind / nz; // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)
  
  // Index offsets
  auto iyp = i.yp();
  auto iym = i.ym();
  int  iyp_2d = iyp.ind / nz; // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)
  int  iym_2d = iym.ind / nz; // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)

  // Use the raw pointers to yup/ydown fields. 
  auto yup = f.f_yup;
  auto ydown = f.f_ydown;
 
  BoutReal gradient_upper = 2.*(yup[iyp.ind] - f[i]) / (dy[ind_2d] + dy[iyp_2d]); // metric->dy[i] is dy
  BoutReal flux_upper = gradient_upper * (J[ind_2d] + J[iyp_2d]) / (g22[ind_2d] + g22[iyp_2d]); //metric->g22[i]
  BoutReal gradient_lower = 2.*(f[i] - ydown[iym.ind]) / (dy[ind_2d] + dy[iyp_2d]);
  BoutReal flux_lower = gradient_lower * (J[ind_2d] + J[iym_2d]) / (g22[ind_2d] + g22[iym_2d]);

  return (flux_upper - flux_lower) / (dy[ind_2d] * J[ind_2d]);
}



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
  
  BoutReal gradient_upper = 2.*(yup[iyp] - f[i]) / (dy + metric->dy[iyp]);
  BoutReal flux_upper = gradient_upper * (J + metric->J[iyp]) / (g_22 + metric->g_22[iyp]);
  
  BoutReal gradient_lower = 2.*(f[i] - ydown[iym]) / (dy + metric->dy[iyp]);
  BoutReal flux_lower = gradient_lower * (J + metric->J[iym]) / (g_22 + metric->g_22[iym]);

  return (flux_upper - flux_lower) / (dy * J);
}

#endif // SINGLE_INDEX_OPS_H
