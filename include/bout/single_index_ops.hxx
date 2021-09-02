
//--  GPU version Update by Dr. Yining Qin, Oct27, 2020

#pragma once
#ifndef SINGLE_INDEX_OPS_H
#define SINGLE_INDEX_OPS_H

#include "field_accessor.hxx"

#ifdef BOUT_HAS_RAJA
//--  RAJA CUDA settings--------------------------------------------------------start
#ifdef BOUT_USE_CUDA
const int CUDA_BLOCK_SIZE = 256; // TODO: Make configurable
using EXEC_POL = RAJA::cuda_exec<CUDA_BLOCK_SIZE>;
#else  // BOUT_ENABLE_CUDA not defined
using EXEC_POL = RAJA::loop_exec;
#endif // defined(BOUT_ENABLE_CUDA)
////-----------CUDA settings------------------------------------------------------end
#endif

// Ind3D: i.zp():
BOUT_HOST_DEVICE inline int i_zp(const int id, const int nz) {
  int jz = id % nz;
  int jzmax = nz - 1;
  return (jz < jzmax) ? (id + 1) : (id - jzmax);
}
// Ind3D: i.zm():
BOUT_HOST_DEVICE inline int i_zm(const int id, const int nz) {
  int jz = id % nz;
  int jzmax = nz - 1;
  return (jz > 0) ? (id - 1) : (id + jzmax);
}

// Ind3D: i.ym():
BOUT_HOST_DEVICE inline int i_yp(const int id, const int nz) { return id + nz; }

// Ind3D: i.yp();
BOUT_HOST_DEVICE inline int i_ym(const int id, const int nz) { return id - nz; }

// Ind3D: i.xp();
BOUT_HOST_DEVICE inline int i_xp(const int id, const int ny, const int nz) {
  return id + ny * nz;
}

// Ind3D: i.xm();
BOUT_HOST_DEVICE inline int i_xm(const int id, const int ny, const int nz) {
  return id - ny * nz;
}

//////////////////////////////////////////////////////////////////////////////
// bracket

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal bracket(const Field2DAccessor<location>& f,
                                         const FieldAccessor<location>& g, const int i) {
  return b0xGrad_dot_Grad(f, g, i) * f.J[i]
         / sqrt(f.g_22[i]); // Dividing by B = sqrt(g_22) / J;
}

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal bracket(const FieldAccessor<location>& f,
                                         const FieldAccessor<location>& g, const int i) {
  const BoutReal* f_a = f.data;
  const BoutReal* g_a = g.data;

  int ny = g.ny;
  int nz = g.nz;

  // Offset indices

  int ixp = i_xp(i, ny, nz);
  int ixm = i_xm(i, ny, nz);
  int izp = i_zp(i, nz);
  int izm = i_zm(i, nz);

  int izpxp = i_xp(izp, ny, nz);
  int izpxm = i_xm(izp, ny, nz);
  int izmxp = i_xp(izm, ny, nz);
  int izmxm = i_xm(izm, ny, nz);

  // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
  BoutReal Jpp = ((f_a[izp] - f_a[izm]) * (g_a[ixp] - g_a[ixm])
                  - (f_a[ixp] - f_a[ixm]) * (g_a[izp] - g_a[izm]));

  // J+x
  BoutReal Jpx =
      (g_a[ixp] * (f_a[izpxp] - f_a[izmxp]) - g_a[ixm] * (f_a[izpxm] - f_a[izmxm])
       - g_a[izp] * (f_a[izpxp] - f_a[izpxm]) + g_a[izm] * (f_a[izmxp] - f_a[izmxm]));

  // Jx+
  BoutReal Jxp =
      (g_a[izpxp] * (f_a[izp] - f_a[ixp]) - g_a[izmxm] * (f_a[ixm] - f_a[izm])
       - g_a[izpxm] * (f_a[izp] - f_a[ixm]) + g_a[izmxp] * (f_a[ixp] - f_a[izm]));

  return (Jpp + Jpx + Jxp) / (12 * f.dx[i] * f.dz[i]);
}

template<CELL_LOC location, class FieldType_f, class FieldType_g>
BOUT_HOST_DEVICE inline BoutReal bracket(const FieldAccessor<location, FieldType_f>& f,
                                         const FieldAccessor<location, FieldType_g>& g, const Ind3D &ind) {
  return bracket(f, g, ind.ind);
}

//////////////////////////////////////////////////////////////////////////////
// DDX

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal DDX(const Field2DAccessor<location>& f, const int i) {
  const int ny = f.ny;

  int i2d = i / f.mesh_nz; // Field2D index
  int ixp = i_xp(i2d, ny, 1);
  int ixm = i_xm(i2d, ny, 1);

  return 0.5 * (f.data[ixp] - f.data[ixm]) / f.dx[i];
}

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal DDX(const FieldAccessor<location>& f, const int i) {

  int ny = f.ny;
  int nz = f.nz;

  int ixp = i_xp(i, ny, nz);
  int ixm = i_xm(i, ny, nz);

  return 0.5 * (f.data[ixp] - f.data[ixm]) / f.dx[i];
}

template<CELL_LOC location, class FieldType>
BOUT_HOST_DEVICE inline BoutReal DDX(const FieldAccessor<location, FieldType>& f, const Ind3D &ind) {
  return DDX(f, ind.ind);
}

//////////////////////////////////////////////////////////////////////////////
// DDY

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal DDY(const Field2DAccessor<location>& f, const int i) {
  const int i2d = i / f.mesh_nz; // Field2D index

  int iyp = i_yp(i2d, 1);
  int iym = i_ym(i2d, 1);

  return 0.5 * (f.yup[iyp] - f.ydown[iym]) / f.dy[i];
}

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal DDY(const FieldAccessor<location>& f, const int i) {
  const int nz = f.nz;

  int iyp = i_yp(i, nz);
  int iym = i_ym(i, nz);

  return 0.5 * (f.yup[iyp] - f.ydown[iym]) / f.dy[i];
}

template<CELL_LOC location, class FieldType>
BOUT_HOST_DEVICE inline BoutReal DDY(const FieldAccessor<location, FieldType>& f, const Ind3D &ind) {
  return DDY(f, ind.ind);
}

//////////////////////////////////////////////////////////////////////////////
// DDZ

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal DDZ(const FieldAccessor<location>& f, const int i) {
  int izp = i_zp(i, f.nz);
  int izm = i_zm(i, f.nz);

  return 0.5 * (f[izp] - f[izm]) / f.dz[i];
}

template<CELL_LOC location, class FieldType>
BOUT_HOST_DEVICE inline BoutReal DDZ(const FieldAccessor<location, FieldType>& f, const Ind3D &ind) {
  return DDZ(f, ind.ind);
}

//////////////////////////////////////////////////////////////////////////////
// Delp2

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal Delp2(const FieldAccessor<location>& f, const int i) {
  const int ny = f.ny;
  const int nz = f.nz;

  const int ixp = i_xp(i, ny, nz);
  const int ixm = i_xm(i, ny, nz);
  const int izp = i_zp(i, nz);
  const int izm = i_zm(i, nz);

  const int izpxp = i_xp(izp, ny, nz);
  const int izpxm = i_xm(izp, ny, nz);
  const int izmxp = i_xp(izm, ny, nz);
  const int izmxm = i_xm(izm, ny, nz);

  const BoutReal dx = f.dx[i];
  const BoutReal dz = f.dz[i];

  return f.G1[i] * (f[ixp] - f[ixm]) / (2.0 * dx)             // DDX
         + f.G3[i] * (f[izp] - f[izm]) / (2.0 * dz)           // DDZ
         + f.g11[i] * (f[ixp] - 2.0 * f[i] + f[ixm]) / SQ(dx) // D2DX2
         + f.g33[i] * (f[izp] - 2.0 * f[i] + f[izm]) / SQ(dz) // D2DZ2
         + 2 * f.g13[i] * ((f[izpxp] - f[izpxm]) - (f[izmxp] - f[izmxm]))
               / (4. * dz * dx); // D2DXDZ
}

template<CELL_LOC location, class FieldType>
BOUT_HOST_DEVICE inline BoutReal Delp2(const FieldAccessor<location, FieldType>& f, const Ind3D &ind) {
  return Delp2(f, ind.ind);
}

//////////////////////////////////////////////////////////////////////////////
// Div_par_Grad_par

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal Div_par_Grad_par(const FieldAccessor<location>& f,
                                                  const int i) {

  const int nz = f.nz;

  const int iyp = i_yp(i, nz);
  const int iym = i_ym(i, nz);

  BoutReal gradient_upper = 2. * (f.yup[iyp] - f[i]) / (f.dy[i] + f.dy[iyp]);
  BoutReal flux_upper = gradient_upper * (f.J[i] + f.J[iyp]) / (f.g_22[i] + f.g_22[iyp]);

  BoutReal gradient_lower = 2. * (f[i] - f.ydown[iym]) / (f.dy[i] + f.dy[iym]);
  BoutReal flux_lower =
      gradient_lower * (f.J[i] + f.J[iym]) / (f.g_22[i] + f.g_22[iym]);

  BoutReal output = (flux_upper - flux_lower) / (f.dy[i] * f.J[i]);
  return output;
}

template<CELL_LOC location, class FieldType>
BOUT_HOST_DEVICE inline BoutReal Div_par_Grad_par(const FieldAccessor<location, FieldType>& f, const Ind3D &ind) {
  return Div_par_Grad_par(f, ind.ind);
}

//////////////////////////////////////////////////////////////////////////////
// b0xGrad_dot_Grad

template <CELL_LOC location>
BOUT_HOST_DEVICE BoutReal b0xGrad_dot_Grad(const FieldAccessor<location>& phi,
                                           const Field2DAccessor<location>& f,
                                           const int i) {
  const int ny = f.ny;

  // Calculate phi derivatives
  BoutReal dpdx = DDX(phi, i);
  BoutReal dpdy = DDY(phi, i);
  BoutReal dpdz = DDZ(phi, i);

  // Calculate advection velocity
  BoutReal vx = phi.g_22[i] * dpdz - phi.g_23[i] * dpdy;
  BoutReal vy = phi.g_23[i] * dpdx - phi.g_12[i] * dpdz;

  // Advect using gradients of f
  // Note: Need to convert indices to 2D

  const int i2d = i / f.mesh_nz;
  const int ixp = i_xp(i2d, ny, 1);
  const int ixm = i_xm(i2d, ny, 1);
  const int iyp = i_yp(i2d, 1);
  const int iym = i_ym(i2d, 1);

  BoutReal vddx = vx * (f[ixp] - f[ixm]) / (2. * f.dx[i]);
  BoutReal vddy = vy * (f.yup[iyp] - f.ydown[iym]) / (2. * f.dy[i]);

  return (vddx + vddy) / (f.J[i] * sqrt(f.g_22[i]));
}

template <CELL_LOC location>
BOUT_HOST_DEVICE BoutReal b0xGrad_dot_Grad(const Field2DAccessor<location>& phi,
                                           const FieldAccessor<location>& f,
                                           const int i) {
  int ny = f.ny;
  int nz = f.nz;

  // Calculate phi derivatives
  BoutReal dpdx = DDX(phi, i);
  BoutReal dpdy = DDY(phi, i);

  // Calculate advection velocity
  BoutReal vx = -f.g_23[i] * dpdy;
  BoutReal vy = f.g_23[i] * dpdx;
  BoutReal vz = f.g_12[i] * dpdy - f.g_22[i] * dpdx;

  int ixp = i_xp(i, ny, nz);
  int ixm = i_xm(i, ny, nz);
  int iyp = i_yp(i, nz);
  int iym = i_ym(i, nz);
  int izp = i_zp(i, nz);
  int izm = i_zm(i, nz);

  BoutReal vddx = vx * (f[ixp] - f[ixm]) / (2. * f.dx[i]);
  BoutReal vddy = vy * (f.yup[iyp] - f.ydown[iym]) / (2. * f.dy[i]);
  BoutReal vddz = vz * (f[izp] - f[izm]) / (2. * f.dz[i]);

  return (vddx + vddy + vddz) / (f.J[i] * sqrt(f.g_22[i]));
}

template<CELL_LOC location, class FieldType_f, class FieldType_g>
BOUT_HOST_DEVICE inline BoutReal b0xGrad_dot_Grad(const FieldAccessor<location, FieldType_f>& f,
                                                  const FieldAccessor<location, FieldType_g>& g, const Ind3D &ind) {
  return b0xGrad_dot_Grad(f, g, ind.ind);
}


//////////////////////////////////////////////////////////////////////////////
// D2DY2

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal D2DY2(const FieldAccessor<location>& f, const int i) {
  const int nz = f.mesh_nz;

  int iyp = i_yp(i, nz);
  int iym = i_ym(i, nz);

  BoutReal result = (f.yup[iyp] - 2. * f[i] + f.ydown[iym]) / SQ(f.dy[i]);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
// Grad_par

template <CELL_LOC location, class FieldType>
BOUT_HOST_DEVICE inline BoutReal Grad_parP(const FieldAccessor<location, FieldType>& f,
                                           const int i) {
  BoutReal result = Grad_par(f, i);

  // if (nonlinear) {
  //   result -= bracket(interp_to(Psi, loc), f, bm_mag) * B0;

  // if (include_rmp) {
  //  result -= bracket(interp_to(rmp_Psi, loc), f, bm_mag) * B0;
  // }
  // }

  return result;
}

template <CELL_LOC location>
BOUT_HOST_DEVICE inline BoutReal Grad_par(const FieldAccessor<location>& f, const int i) {
  return DDY(f, i) / sqrt(f.g_22[i]);
}

template<CELL_LOC location, class FieldType>
BOUT_HOST_DEVICE inline BoutReal Grad_par(const FieldAccessor<location, FieldType>& f, const Ind3D &ind) {
  return Grad_par(f, ind.ind);
}

#endif // SINGLE_INDEX_OPS_H
