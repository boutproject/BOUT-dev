#pragma once
#ifndef BOUT_STENCIL_EXPR_HXX
#define BOUT_STENCIL_EXPR_HXX

#include "bout/coordinates_accessor.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldops.hxx"
#include "bout/mesh.hxx"
#include "bout/single_index_ops.hxx"

namespace bout::stencil {

struct DDZ_C2_Op {
  CoordinatesAccessor coords;
  int nz{0};

  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& lhs,
                                                        const RView&) const {
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);

    return 0.5 * (lhs(izp) - lhs(izm)) / coords.dz(idx);
  }
};

struct DDZ_C4_Op {
  CoordinatesAccessor coords;
  int nz{0};

  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& lhs,
                                                        const RView&) const {
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);
    const int izp2 = i_zp(izp, nz);
    const int izm2 = i_zm(izm, nz);

    return (-lhs(izp2) + 8.0 * lhs(izp) - 8.0 * lhs(izm) + lhs(izm2))
           / (12.0 * coords.dz(idx));
  }
};

struct BracketArakawaOp {
  CoordinatesAccessor coords;
  int ny{0};
  int nz{0};

  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& f,
                                                        const RView& g) const {
    const int ixp = i_xp(idx, ny, nz);
    const int ixm = i_xm(idx, ny, nz);
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);

    const int izpxp = i_xp(izp, ny, nz);
    const int izpxm = i_xm(izp, ny, nz);
    const int izmxp = i_xp(izm, ny, nz);
    const int izmxm = i_xm(izm, ny, nz);

    // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
    const BoutReal Jpp =
        ((f(izp) - f(izm)) * (g(ixp) - g(ixm)) - (f(ixp) - f(ixm)) * (g(izp) - g(izm)));

    // J+x
    const BoutReal Jpx =
        (g(ixp) * (f(izpxp) - f(izmxp)) - g(ixm) * (f(izpxm) - f(izmxm))
         - g(izp) * (f(izpxp) - f(izpxm)) + g(izm) * (f(izmxp) - f(izmxm)));

    // Jx+
    const BoutReal Jxp = (g(izpxp) * (f(izp) - f(ixp)) - g(izmxm) * (f(ixm) - f(izm))
                          - g(izpxm) * (f(izp) - f(ixm)) + g(izmxp) * (f(ixp) - f(izm)));

    return (Jpp + Jpx + Jxp) / (12.0 * coords.dx(idx) * coords.dz(idx));
  }
};

using DDZExprC2 = BinaryExpr<Field3D, Field3D, Field3D, DDZ_C2_Op>;
using DDZExprC4 = BinaryExpr<Field3D, Field3D, Field3D, DDZ_C4_Op>;
using BracketArakawaExpr = BinaryExpr<Field3D, Field3D, Field3D, BracketArakawaOp>;

} // namespace bout::stencil

inline bout::stencil::DDZExprC2 DDZ_C2(const Field3D& f) {
  checkData(f);

  const auto region_id = f.getMesh()->getRegionID("RGN_NOBNDRY");

  return bout::stencil::DDZExprC2{
      static_cast<Field3D::View>(f),
      static_cast<Field3D::View>(f),
      bout::stencil::DDZ_C2_Op{CoordinatesAccessor{f.getCoordinates()}, f.getNz()},
      f.getMesh(),
      f.getLocation(),
      f.getDirections(),
      region_id,
      f.getMesh()->getRegion("RGN_NOBNDRY")};
}

inline bout::stencil::DDZExprC4 DDZ_C4(const Field3D& f) {
  checkData(f);

  const auto region_id = f.getMesh()->getRegionID("RGN_NOBNDRY");

  return bout::stencil::DDZExprC4{
      static_cast<Field3D::View>(f),
      static_cast<Field3D::View>(f),
      bout::stencil::DDZ_C4_Op{CoordinatesAccessor{f.getCoordinates()}, f.getNz()},
      f.getMesh(),
      f.getLocation(),
      f.getDirections(),
      region_id,
      f.getMesh()->getRegion("RGN_NOBNDRY")};
}

inline bout::stencil::BracketArakawaExpr bracket_arakawa(const Field3D& f,
                                                         const Field3D& g) {
  checkData(f);
  checkData(g);
  ASSERT1_FIELDS_COMPATIBLE(f, g);

  const auto region_id = f.getMesh()->getRegionID("RGN_NOBNDRY");

  return bout::stencil::BracketArakawaExpr{
      static_cast<Field3D::View>(f),
      static_cast<Field3D::View>(g),
      bout::stencil::BracketArakawaOp{CoordinatesAccessor{f.getCoordinates()}, f.getNy(),
                                      f.getNz()},
      f.getMesh(),
      f.getLocation(),
      f.getDirections(),
      region_id,
      f.getMesh()->getRegion("RGN_NOBNDRY")};
}

#endif
