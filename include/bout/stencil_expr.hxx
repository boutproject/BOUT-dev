#pragma once
#ifndef BOUT_STENCIL_EXPR_HXX
#define BOUT_STENCIL_EXPR_HXX

#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/build_config.hxx"
#include "bout/coordinates_accessor.hxx"
#include "bout/field.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldops.hxx"
#include "bout/mesh.hxx"
#include "bout/single_index_ops.hxx"
#include "bout/utils.hxx"

#include <string>

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

struct DDZ_Dispatch_Op {
  CoordinatesAccessor coords;
  int nz{0};
  STAGGER stagger{STAGGER::None};
  DIFF_METHOD method{DIFF_DEFAULT};

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c2_none(int idx,
                                                           const LView& lhs) const {
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);

    return 0.5 * (lhs(izp) - lhs(izm)) / coords.dz(idx);
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c4_none(int idx,
                                                           const LView& lhs) const {
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);
    const int izp2 = i_zp(izp, nz);
    const int izm2 = i_zm(izm, nz);

    return (-lhs(izp2) + 8.0 * lhs(izp) - 8.0 * lhs(izm) + lhs(izm2))
           / (12.0 * coords.dz(idx));
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c2_c2l(int idx,
                                                          const LView& lhs) const {
    const int izm = i_zm(idx, nz);

    return (lhs(idx) - lhs(izm)) / coords.dz(idx);
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c4_c2l(int idx,
                                                          const LView& lhs) const {
    const int izm = i_zm(idx, nz);
    const int izp = i_zp(idx, nz);
    const int izm2 = i_zm(izm, nz);

    return (27.0 * (lhs(idx) - lhs(izm)) - (lhs(izp) - lhs(izm2)))
           / (24.0 * coords.dz(idx));
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c2_l2c(int idx,
                                                          const LView& lhs) const {
    const int izp = i_zp(idx, nz);

    return (lhs(izp) - lhs(idx)) / coords.dz(idx);
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c4_l2c(int idx,
                                                          const LView& lhs) const {
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);
    const int izp2 = i_zp(izp, nz);

    return (27.0 * (lhs(izp) - lhs(idx)) - (lhs(izp2) - lhs(izm)))
           / (24.0 * coords.dz(idx));
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c2(int idx, const LView& lhs) const {
    switch (stagger) {
    case STAGGER::None:
      return apply_c2_none(idx, lhs);
    case STAGGER::C2L:
      return apply_c2_c2l(idx, lhs);
    case STAGGER::L2C:
      return apply_c2_l2c(idx, lhs);
    }
    return 0.0;
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c4(int idx, const LView& lhs) const {
    switch (stagger) {
    case STAGGER::None:
      return apply_c4_none(idx, lhs);
    case STAGGER::C2L:
      return apply_c4_c2l(idx, lhs);
    case STAGGER::L2C:
      return apply_c4_l2c(idx, lhs);
    }
    return 0.0;
  }

  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& lhs,
                                                        const RView&) const {
    switch (method) {
    case DIFF_C2:
      return apply_c2(idx, lhs);
    case DIFF_C4:
      return apply_c4(idx, lhs);
    default:
      return 0.0;
    }
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
using DDZDispatchExpr = BinaryExpr<Field3D, Field3D, Field3D, DDZ_Dispatch_Op>;
using BracketArakawaExpr = BinaryExpr<Field3D, Field3D, Field3D, BracketArakawaOp>;

} // namespace bout::stencil

namespace {

inline DIFF_METHOD parseDDZMethodString(const std::string& method) {
  auto normalized = uppercase(method);
  if (normalized.rfind("DIFF_", 0) == 0) {
    normalized = normalized.substr(5);
  }

  if (normalized == "DEFAULT") {
    return DIFF_DEFAULT;
  }
  if (normalized == "C2") {
    return DIFF_C2;
  }
  if (normalized == "C4") {
    return DIFF_C4;
  }
  if (normalized == "U1") {
    return DIFF_U1;
  }
  if (normalized == "U2") {
    return DIFF_U2;
  }
  if (normalized == "U3") {
    return DIFF_U3;
  }
  if (normalized == "W2") {
    return DIFF_W2;
  }
  if (normalized == "W3") {
    return DIFF_W3;
  }
  if (normalized == "S2") {
    return DIFF_S2;
  }
  if (normalized == "FFT") {
    return DIFF_FFT;
  }
  if (normalized == "SPLIT") {
    return DIFF_SPLIT;
  }

  throw BoutException("Unknown DDZ method '{:s}'", method);
}

} // namespace

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

inline bout::stencil::DDZDispatchExpr DDZ(const Field3D& f,
                                          CELL_LOC outloc = CELL_DEFAULT,
                                          DIFF_METHOD method = DIFF_DEFAULT,
                                          const std::string& region = "RGN_NOBNDRY") {
  checkData(f);

  const auto resolved_outloc = (outloc == CELL_DEFAULT) ? f.getLocation() : outloc;
  const auto stagger =
      f.getMesh()->getStagger(f.getLocation(), resolved_outloc, CELL_ZLOW);

  if (method == DIFF_DEFAULT) {
    method = f.getMesh()->getDefaultMethod(DIRECTION::Z, DERIV::Standard, stagger);
  }

  if ((method != DIFF_C2) && (method != DIFF_C4)) {
    throw BoutException("DDZ only supports DIFF_C2 and DIFF_C4, got {:s}",
                        toString(method));
  }

  const auto region_id = f.getMesh()->getRegionID(region);

  return bout::stencil::DDZDispatchExpr{
      static_cast<Field3D::View>(f),
      static_cast<Field3D::View>(f),
      bout::stencil::DDZ_Dispatch_Op{
          CoordinatesAccessor{f.getCoordinates(resolved_outloc)}, f.getNz(), stagger,
          method},
      f.getMesh(),
      resolved_outloc,
      f.getDirections(),
      region_id,
      f.getMesh()->getRegion(region)};
}

inline bout::stencil::DDZDispatchExpr DDZ(const Field3D& f, CELL_LOC outloc,
                                          const std::string& method,
                                          const std::string& region = "RGN_NOBNDRY") {
  return DDZ(f, outloc, parseDDZMethodString(method), region);
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
