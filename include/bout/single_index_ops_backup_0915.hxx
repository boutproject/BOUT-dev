#pragma once
#ifndef SINGLE_INDEX_OPS_H
#define SINGLE_INDEX_OPS_H

#include "field_accessor.hxx"

template<IND_TYPE N, CELL_LOC location>
BoutReal bracket(const FieldAccessor<location> &f, const FieldAccessor<location> &g, const SpecificInd<N> &ind) {
  Coordinates *metric = g.coords;

 
  auto ixp = ind.xp();
  auto ixm = ind.xm();
  auto izp = ind.zp();
  auto izm = ind.zm();

  auto izpxp = izp.xp();
  auto izpxm = izp.xm();
  auto izmxp = izm.xp();
  auto izmxm = izm.xm();
  

  BoutReal Jpp = ((f[izp] - f[izm]) * (g[ixp] - g[ixm])
                  - (f[ixp] - f[ixm]) * (g[izp] - g[izm]));
  
 
  BoutReal Jpx = (g[ixp] * (f[izpxp] - f[izmxp]) - g[ixm] * (f[izpxm] - f[izmxm])
                  - g[izp] * (f[izpxp] - f[izpxm]) + g[izm] * (f[izmxp] - f[izmxm]));

 
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
BoutReal DDZ(const FieldAccessor<location> &f, const SpecificInd<N> &ind) {
  return (f[ind.zp()] - f[ind.zm()]) / (2. * f.coords->dz);
}

template<IND_TYPE N, CELL_LOC location>
BoutReal Delp2(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
  Coordinates *metric = f.coords;

 
  auto izm = i.zm();
  auto izp = i.zp();
  auto ixm = i.xm();
  auto ixp = i.xp();
 
  return metric->G1[i] * (f[ixp] - f[ixm]) / (2.*metric->dx[i])  
    + metric->G3[i] * (f[izp] - f[izm]) / (2.*metric->dz)  
    + metric->g11[i] * (f[ixp] - 2.*f[i] + f[ixm]) / SQ(metric->dx[i])  
    + metric->g33[i] * (f[izp] - 2.*f[i] + f[izm]) / SQ(metric->dz) 
    + 2 * metric->g13[i] * ((f[izp.xp()] - f[izp.xm()]) -
                            (f[izm.xp()] - f[izm.xm()])) / (4. * metric->dz * metric->dx[i]) 
;
}

template<IND_TYPE N, CELL_LOC location>
BoutReal Div_par_Grad_par(const FieldAccessor<location> &f, const SpecificInd<N> &i) {
  Coordinates *metric = f.coords;
  
  auto iyp = i.yp();
  auto iym = i.ym();

 
  const Field3D &yup = *f.yup;
  const Field3D &ydown = *f.ydown;

 
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
