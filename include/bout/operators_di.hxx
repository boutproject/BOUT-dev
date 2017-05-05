
#pragma once

#ifndef __OPERATORS_DI_H__
#define __OPERATORS_DI_H__

/// \brief 2nd order central differencing in X
///
/// Performs calculation at a single point. The input field
/// must be defined at i.xp() and i.xm()
///
/// @param f  The field to be differentiated
/// @param i  The point where the result is calculated
///
inline BoutReal DDX_C2(const Field3D &f, const DataIterator &i) {
  return (f[i.xp()] - f[i.xm()])/(2.*mesh->coordinates()->dx[i]);
}

/// \brief 2nd order central differencing in Y using yup/down fields
///
/// @param f  The field to be differentiated
/// @param i  The point where the result is calculated
/// 
inline BoutReal DDY_C2(const Field3D &f, const DataIterator &i) {
  return (f.yup()[i.yp()] - f.ydown()[i.ym()])/(2.*mesh->coordinates()->dy[i]);
}

/// \brief 2nd order central differencing in Y assuming field aligned
///
/// @param f  The field to be differentiated (field aligned)
/// @param i  The point where the result is calculated
/// 
inline BoutReal DDY_C2_FA(const Field3D &f, const DataIterator &i) {
  return (f[i.yp()] - f[i.ym()])/(2.*mesh->coordinates()->dy[i]);
}

/// \brief 2nd order central differencing in Z 
///
/// @param f  The field to be differentiated
/// @param i  The point where the result is calculated
/// 
inline BoutReal DDZ_C2(const Field3D &f, const DataIterator &i) {
  return (f[i.zp()] - f[i.zm()])/(2.*mesh->coordinates()->dz);
}


///
/// d/dx( g^xx df/dx + g^xz df/dz) + d/dz( g^zz df/dz + g^xz df/dx)
/// 
BoutReal Delp2_C2(const Field3D &f, const DataIterator &i) {
  Coordinates *metric = mesh->coordinates();

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
/// @param[in[ i  The point where the result is calculated
///
BoutReal bracket_arakawa(const Field3D &f, const Field3D &g, const DataIterator &i) {
  Coordinates *metric = mesh->coordinates();
  
  BoutReal fxp = f[i.xp()];
  BoutReal fxm = f[i.xm()];
  BoutReal fzp = f[i.zp()];
  BoutReal fzm = f[i.zm()];

  BoutReal fpp = f[i.offset(1,0,1)];
  BoutReal fpm = f[i.offset(1,0,-1)];
  BoutReal fmp = f[i.offset(-1,0,1)];
  BoutReal fmm = f[i.offset(-1,0,-1)];
  
  BoutReal gxp = g[i.xp()];
  BoutReal gxm = g[i.xm()];
  BoutReal gzp = g[i.zp()];
  BoutReal gzm = g[i.zm()];

  // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
  BoutReal Jpp = 0.25*( (fzp - fzm)*
                        (gxp - gxm) -
                        (fxp - fxm)*
                        (gzp - gzm) );
  
  // J+x
  BoutReal Jpx = 0.25*( gxp*(fpp-fpm) -
                        gxm*(fmp-fmm) -
                        gzp*(fpp-fmp) +
                        gzm*(fpm-fmm));
  // Jx+
  BoutReal Jxp = 0.25*( g[i.offset(1,0,1)]*(fzp-fxp) -
                        g[i.offset(-1,0,-1)]*(fxm-fzm) -
                        g[i.offset(-1,0,1)]*(fzp-fxm) +
                        g[i.offset(1,0,-1)]*(fxp-fzm) );
  
  return (Jpp + Jpx + Jxp) / (3. * metric->dx[i] * metric->dz);
}

#endif // __OPERATORS_DI_H__
