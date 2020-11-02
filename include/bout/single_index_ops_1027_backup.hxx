
#pragma once
#ifndef SINGLE_INDEX_OPS_H
#define SINGLE_INDEX_OPS_H

#include "field_accessor.hxx"

#define BOUT_HOST_DEVICE __host__ __device__
#define RAJA_HOST_DEVICE __host__ __device__

#define BOUT_HAS_UMPIRE
#define BOUT_USE_CUDA
#define gpu_input(name) gpu_input_helper(#name, name )

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


// whole setting up
__managed__  BoutReal* gpu_n_ddt;     //  ddt(n) to __device__
__managed__  BoutReal* gpu_vort_ddt;  // ddt(vort) to __device__

__managed__  int t_nx = 0;  // Mesh x size
__managed__  int t_ny = 0;  // Mesh y size
__managed__  int t_nz = 0;   // Mesh z size


// t_acc
__managed__ BoutReal* t_dx = nullptr;  // pointer to Field2D data
__managed__ BoutReal* t_dy = nullptr;  // pointer to Field2D data
__managed__ BoutReal t_dz = 0.0;       // dz of Field2D
__managed__ BoutReal* t_J = nullptr;   // pointer to Field2D data

__managed__ BoutReal* t_G1 = nullptr;         // Field2D
__managed__ BoutReal* t_G3 = nullptr;         // Field2D
__managed__ BoutReal* t_g11 = nullptr;        // Field2D
__managed__ BoutReal* t_g13 = nullptr;        // Field2D
__managed__ BoutReal* t_g33 = nullptr;        // Field2D
__managed__ BoutReal* t_g22 = nullptr;        // Field2D


//n_acc
__managed__ BoutReal* n_dx = nullptr;  // pointer to Field2D data
__managed__ BoutReal* n_dy = nullptr;  // pointer to Field2D data
__managed__ BoutReal n_dz = 0.0;       // dz of Field2D
__managed__ BoutReal* n_J = nullptr;   // pointer to Field2D data

__managed__ BoutReal* n_G1 = nullptr;         // Field2D
__managed__ BoutReal* n_G3 = nullptr;         // Field2D
__managed__ BoutReal* n_g11 = nullptr;        // Field2D
__managed__ BoutReal* n_g13 = nullptr;        // Field2D
__managed__ BoutReal* n_g33 = nullptr;        // Field2D
__managed__ BoutReal* n_g22 = nullptr;        // Field2D

//vort_acc

__managed__ BoutReal* v_dx = nullptr;  // pointer to Field2D data
__managed__ BoutReal* v_dy = nullptr;  // pointer to Field2D data
__managed__ BoutReal  v_dz = 0.0;       // dz of Field2D
__managed__ BoutReal* v_J = nullptr;   // pointer to Field2D data

__managed__ BoutReal* v_G1 = nullptr;         // Field2D
__managed__ BoutReal* v_G3 = nullptr;         // Field2D
__managed__ BoutReal* v_g11 = nullptr;        // Field2D
__managed__ BoutReal* v_g13 = nullptr;        // Field2D
__managed__ BoutReal* v_g33 = nullptr;        // Field2D
__managed__ BoutReal* v_g22 = nullptr;        // Field2D




//phi_minus_n_acc
 
__managed__ BoutReal* pmn_yup = nullptr;
__managed__ BoutReal* pmn_ydown = nullptr;


__managed__ BoutReal* pmn_dx = nullptr;  // pointer to Field2D data
__managed__ BoutReal* pmn_dy = nullptr;  // pointer to Field2D data
__managed__ BoutReal  pmn_dz = 0.0;       // dz of Field2D
__managed__ BoutReal* pmn_J = nullptr;   // pointer to Field2D data

__managed__ BoutReal* pmn_G1 = nullptr;         // Field2D
__managed__ BoutReal* pmn_G3 = nullptr;         // Field2D
__managed__ BoutReal* pmn_g11 = nullptr;        // Field2D
__managed__ BoutReal* pmn_g13 = nullptr;        // Field2D
__managed__ BoutReal* pmn_g33 = nullptr;        // Field2D
__managed__ BoutReal* pmn_g22 = nullptr;        // Field2D



template<CELL_LOC location>
void metrics (Field3D &n, Field3D &vort,  Field3D &phi, Field3D &phi_minus_n, 
			const FieldAccessor<location> &n_acc, const FieldAccessor<location> &vort_acc,
			const FieldAccessor<> &phi_acc, const FieldAccessor<> &phi_minus_n_acc )
 {

    gpu_n_ddt= const_cast<BoutReal*>(n.timeDeriv()->operator()(0,0)); // copy ddt(n) to __device__
    gpu_vort_ddt = const_cast<BoutReal*>(vort.timeDeriv()->operator()(0,0)); // copy ddt(vort) to __device__

    t_nx = n.getNx();
    t_ny = n.getNy();
    t_nz = n.getNz();

  // n_acc 
    Coordinates* metric = n_acc.coords;
    n_G1= &metric->G1(0,0);
    n_G3 = &metric->G3(0,0);
 
    n_g11 = &metric->g11(0,0);
    n_g13 = &metric->g13(0,0);
    n_g33 = &metric->g33(0,0);
    n_g22 = &metric->g22(0,0);

    n_dx = &metric->dx(0,0);
    n_dy = &metric->dy(0,0);
    n_dz = metric->dz;
    n_J  = &metric->J(0,0);
// vort_acc

    metric = vort_acc.coords;
    v_G1= &metric->G1(0,0);
    v_G3 = &metric->G3(0,0);
 
    v_g11 = &metric->g11(0,0);
    v_g13 = &metric->g13(0,0);
    v_g33 = &metric->g33(0,0);
    v_g22 = &metric->g22(0,0);

    v_dx = &metric->dx(0,0);
    v_dy = &metric->dy(0,0);
    v_dz = metric->dz;
    v_J  = &metric->J(0,0);
// t_acc
    
    t_G1 = v_G1;
    t_G3 = v_G3;
    t_g11 = v_g11;
    t_g13 = v_g13;
    t_g33 = v_g33;
    t_g22 = v_g22;

    t_dx = v_dx;
    t_dy = v_dy;
    t_dz = v_dz;
    t_J = v_J;


// phi_minus_n_acc 
   metric = phi_minus_n_acc.coords;
    const Field3D &yup = *phi_minus_n_acc.yup;
    const Field3D &ydown = *phi_minus_n_acc.ydown;
    pmn_yup = const_cast<BoutReal*>(yup(0,0));
    pmn_ydown = const_cast<BoutReal*>(ydown(0,0));

    pmn_dy = &metric->dy(0,0);
    pmn_J = &metric->J(0,0);
    pmn_g22 = &metric->g22(0,0);

}

template<CELL_LOC location>
BOUT_DEVICE static void gpu_input_helper(const char* p_name, const FieldAccessor<location> &f ){

if (p_name == "n_acc"){
    
    t_G1 = n_G1;
    t_G3 = n_G3;
    t_g11 = n_g11;
    t_g13 = n_g13;
    t_g33 = n_g33;
    t_g22 = n_g22;

    t_dx = n_dx;
    t_dy = n_dy;
    t_dz = n_dz;
    t_J = n_J;
}


if (p_name == "vort_acc"){

    
    t_G1 = v_G1;
    t_G3 = v_G3;
    t_g11 = v_g11;
    t_g13 = v_g13;
    t_g33 = v_g33;
    t_g22 = v_g22;

    t_dx = v_dx;
    t_dy = v_dy;
    t_dz = v_dz;
    t_J = v_J;


}

 

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
 
  int  ind_2d = ind.ind / t_nz;  // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)
 
  // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
  BoutReal Jpp = ((f[izp] - f[izm]) * (g[ixp] - g[ixm])
                  - (f[ixp] - f[ixm]) * (g[izp] - g[izm]));
  
  // J+x
  BoutReal Jpx = (g[ixp] * (f[izpxp] - f[izmxp]) - g[ixm] * (f[izpxm] - f[izmxm])
                  - g[izp] * (f[izpxp] - f[izpxm]) + g[izm] * (f[izmxp] - f[izmxm]));

  // Jx+
  BoutReal Jxp = (g[izpxp] * (f[izp] - f[ixp]) - g[izmxm] * (f[ixm] - f[izm])
                  - g[izpxm] * (f[izp] - f[ixm]) + g[izmxp] * (f[ixp] - f[izm]));
  
  return (Jpp + Jpx + Jxp) / (12 * t_dx[ind_2d] * t_dz);
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
  return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
}

template<IND_TYPE N, CELL_LOC location>
BoutReal  DDZ(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
}


template<IND_TYPE N, CELL_LOC location>
RAJA_HOST_DEVICE BoutReal Delp2_g(const FieldAccessor<location> &f, const SpecificInd<N> &i) {

gpu_input(f);

  auto izm = i.zm();
  auto izp = i.zp();
  auto ixm = i.xm();
  auto ixp = i.xp();
  int ind_2d = i.ind / t_nz;  // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)

  BoutReal DDX_g = t_G1[ind_2d] * (f[ixp] - f[ixm]) / (2.0 * t_dx[ind_2d]); //DDX
  BoutReal DDZ_g =  t_G3[ind_2d] * (f[izp] - f[izm]) / (2.0 * t_dz);  // DDZ
  BoutReal D2DX2_g =  t_g11[ind_2d] * (f[ixp] - 2.*f[i] + f[ixm]) / (t_dx[ind_2d]*t_dx[ind_2d]);  // D2DX2
  BoutReal D2DZ2_g = t_g33[ind_2d] * (f[izp] - 2.0*f[i] + f[izm]) /(t_dz*t_dz); // D2DZ2 
  BoutReal D2DXDZ_g = 2 * t_g13[ind_2d] * ((f[izp.xp()] - f[izp.xm()]) - (f[izm.xp()] - f[izm.xm()])) / (4. * t_dz * t_dx[ind_2d]) ; // D2DXDZ
  
  BoutReal res = DDX_g + DDZ_g + D2DX2_g + D2DZ2_g + D2DXDZ_g;
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
  
  // Index offsets
  auto iyp = i.yp();
  auto iym = i.ym();
int  ind_2d = i.ind / t_nz; // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)
  int  iyp_2d = iyp.ind / t_nz; // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)
  int  iym_2d = iym.ind / t_nz; // Map3Dto2Da: index for metrics(Field2D) data (has no z-dependence)

 // Use the raw pointers to yup/ydown fields by t_yup and t_down. They both use 3D index;
  // Fetch values used more than once
  BoutReal dy = pmn_dy[ind_2d];
  BoutReal J  = pmn_J[ind_2d];
  BoutReal g_22 = pmn_g22[ind_2d];
 
//  printf("f[i] %f ,pmn_yup[iyp.ind] %f \n", f[i], pmn_yup[iyp.ind]);
 
  BoutReal gradient_upper = 2.*(pmn_yup[iyp.ind] - f[i]) / (dy + pmn_dy[iyp_2d]); // metric->dy[i] is t_dy
  BoutReal flux_upper = gradient_upper * (J + pmn_J[iyp_2d]) / (g_22 + pmn_g22[iyp_2d]); //metric->g22[i]
  BoutReal gradient_lower = 2.*(f[i] - pmn_ydown[iym.ind]) / (dy + pmn_dy[iyp_2d]);
  BoutReal flux_lower = gradient_lower * (J + pmn_J[iym_2d]) / (g_22 + pmn_g22[iym_2d]);

  return (flux_upper - flux_lower) / (dy * J);
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
