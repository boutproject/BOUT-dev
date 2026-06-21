#pragma once
#ifndef BOUT_STENCIL_EXPR_HXX
#define BOUT_STENCIL_EXPR_HXX

#include "bout/coordinates_accessor.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldops.hxx"
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

using DDZExprC2 = BinaryExpr<Field3D, Field3D, Field3D, DDZ_C2_Op>;
using DDZExprC4 = BinaryExpr<Field3D, Field3D, Field3D, DDZ_C4_Op>;

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

#endif
