
#pragma once

#ifndef __OPERATORS_DI_H__
#define __OPERATORS_DI_H__

#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/mesh.hxx>

/// \brief 2nd order central differencing in X
///
/// Performs calculation at a single point. The input field
/// must be defined at i.xp() and i.xm()
///
/// @param f  The field to be differentiated
/// @param i  The point where the result is calculated
///
inline BoutReal DDX_C2(const Field3D &f, const DataIterator &i) {
  return (f[i.xp()] - f[i.xm()])/(2.*f.getCoordinates()->dx[i]);
}

/// \brief 2nd order central differencing in Y using yup/down fields
///
/// @param f  The field to be differentiated
/// @param i  The point where the result is calculated
/// 
inline BoutReal DDY_C2(const Field3D &f, const DataIterator &i) {
  return (f.yup()[i.yp()] - f.ydown()[i.ym()])/(2.*f.getCoordinates()->dy[i]);
}

/// \brief 2nd order central differencing in Y assuming field aligned
///
/// @param f  The field to be differentiated (field aligned)
/// @param i  The point where the result is calculated
/// 
inline BoutReal DDY_C2_FA(const Field3D &f, const DataIterator &i) {
  return (f[i.yp()] - f[i.ym()])/(2.*f.getCoordinates()->dy[i]);
}

/// \brief 2nd order central differencing in Z 
///
/// @param f  The field to be differentiated
/// @param i  The point where the result is calculated
/// 
inline BoutReal DDZ_C2(const Field3D &f, const DataIterator &i) {
  return (f[i.zp()] - f[i.zm()])/(2.*f.getCoordinates()->dz);
}


///
/// d/dx( g^xx df/dx + g^xz df/dz) + d/dz( g^zz df/dz + g^xz df/dx)
/// 
BoutReal Delp2_C2(const Field3D &f, const DataIterator &i) {
  Coordinates *metric = f.getCoordinates();

  //  o -- UZ -- o
  //  |          |
  //  LX        RX
  //  |          |
  //  o -- DZ -- o
  
  // Upper Z boundary (UZ)
  
  
}

/// \brief Arakawa bracket in X-Z [f,g]
///
/// @param[in] f  Field to be differentiated
/// @param[in] g  Field to be differentiated
/// @param[in] i  The point where the result is calculated
BoutReal bracket_arakawa(const Field3D &f, const Field3D &g, const DataIterator &i) {
  Coordinates *metric = f.getCoordinates();

  // Indexing
  const auto xp = i.xp();
  const auto xm = i.xm();
  const auto zp = i.zp();
  const auto zm = i.zm();

  const auto xpzp = i.offset(1,0,1);
  const auto xpzm = i.offset(1,0,-1);
  const auto xmzp = i.offset(-1,0,1);
  const auto xmzm = i.offset(-1,0,-1);
  
  const BoutReal fxp = f[xp];
  const BoutReal fxm = f[xm];
  const BoutReal fzp = f[zp];
  const BoutReal fzm = f[zm];

  const BoutReal fpp = f[xpzp];
  const BoutReal fpm = f[xpzm];
  const BoutReal fmp = f[xmzp];
  const BoutReal fmm = f[xmzm];
  
  const BoutReal gxp = g[xp];
  const BoutReal gxm = g[xm];
  const BoutReal gzp = g[zp];
  const BoutReal gzm = g[zm];

  // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
  BoutReal Jpp =
    (fzp - fzm) * (gxp - gxm)
    - (fxp - fxm) * (gzp - gzm);

  // J+x
  BoutReal Jpx =
    gxp * (fpp - fpm)
    - gxm * (fmp - fmm)
    - gzp * (fpp - fmp)
    + gzm * (fpm - fmm);

  // Jx+
  BoutReal Jxp =
    g[xpzp] * (fzp - fxp)
    - g[xmzm] * (fxm - fzm)
    - g[xmzp] * (fzp - fxm)
    + g[xpzm] * (fxp - fzm);

  return (Jpp + Jpx + Jxp) / (12. * metric->dx[i] * metric->dz);
}

#endif // __OPERATORS_DI_H__
