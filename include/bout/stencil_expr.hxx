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
#include "bout/index_derivs.hxx"
#include "bout/mesh.hxx"
#include "bout/single_index_ops.hxx"
#include "bout/utils.hxx"

#include <cmath>
#include <string>

namespace bout::stencil {

struct DDX_C2_Op {
  CoordinatesAccessor coords;
  int ny{0};
  int nz{0};
  bool inc_int_shear{false};

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_ddz_c2(int idx,
                                                          const LView& lhs) const {
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);

    return 0.5 * (lhs(izp) - lhs(izm)) / coords.dz(idx);
  }

  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& lhs,
                                                        const RView&) const {
    const int ixp = i_xp(idx, ny, nz);
    const int ixm = i_xm(idx, ny, nz);

    BoutReal result = 0.5 * (lhs(ixp) - lhs(ixm)) / coords.dx(idx);
    if (inc_int_shear) {
      result += coords.IntShiftTorsion(idx) * apply_ddz_c2(idx, lhs);
    }
    return result;
  }
};

struct DDX_C4_Op {
  CoordinatesAccessor coords;
  int ny{0};
  int nz{0};
  bool inc_int_shear{false};

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_ddz_c4(int idx,
                                                          const LView& lhs) const {
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);
    const int izp2 = i_zp(izp, nz);
    const int izm2 = i_zm(izm, nz);

    return (-lhs(izp2) + 8.0 * lhs(izp) - 8.0 * lhs(izm) + lhs(izm2))
           / (12.0 * coords.dz(idx));
  }

  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& lhs,
                                                        const RView&) const {
    const int ixp = i_xp(idx, ny, nz);
    const int ixm = i_xm(idx, ny, nz);
    const int ixp2 = i_xp(ixp, ny, nz);
    const int ixm2 = i_xm(ixm, ny, nz);

    BoutReal result = (-lhs(ixp2) + 8.0 * lhs(ixp) - 8.0 * lhs(ixm) + lhs(ixm2))
                      / (12.0 * coords.dx(idx));
    if (inc_int_shear) {
      result += coords.IntShiftTorsion(idx) * apply_ddz_c4(idx, lhs);
    }
    return result;
  }
};

struct DDX_Dispatch_Op {
  CoordinatesAccessor coords;
  int ny{0};
  int nz{0};
  STAGGER stagger{STAGGER::None};
  DIFF_METHOD method{DIFF_DEFAULT};
  bool inc_int_shear{false};

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_ddz_c2(int idx,
                                                          const LView& lhs) const {
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);

    return 0.5 * (lhs(izp) - lhs(izm)) / coords.dz(idx);
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_ddz_c4(int idx,
                                                          const LView& lhs) const {
    const int izp = i_zp(idx, nz);
    const int izm = i_zm(idx, nz);
    const int izp2 = i_zp(izp, nz);
    const int izm2 = i_zm(izm, nz);

    return (-lhs(izp2) + 8.0 * lhs(izp) - 8.0 * lhs(izm) + lhs(izm2))
           / (12.0 * coords.dz(idx));
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c2_none(int idx,
                                                           const LView& lhs) const {
    const int ixp = i_xp(idx, ny, nz);
    const int ixm = i_xm(idx, ny, nz);

    return 0.5 * (lhs(ixp) - lhs(ixm)) / coords.dx(idx);
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c4_none(int idx,
                                                           const LView& lhs) const {
    const int ixp = i_xp(idx, ny, nz);
    const int ixm = i_xm(idx, ny, nz);
    const int ixp2 = i_xp(ixp, ny, nz);
    const int ixm2 = i_xm(ixm, ny, nz);

    return (-lhs(ixp2) + 8.0 * lhs(ixp) - 8.0 * lhs(ixm) + lhs(ixm2))
           / (12.0 * coords.dx(idx));
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_w2_none(int idx,
                                                           const LView& lhs) const {
    const int ixp = i_xp(idx, ny, nz);
    const int ixm = i_xm(idx, ny, nz);

    const BoutReal dc = 0.5 * (lhs(ixp) - lhs(ixm));
    const BoutReal dl = lhs(idx) - lhs(ixm);
    const BoutReal dr = lhs(ixp) - lhs(idx);

    const BoutReal isl = SQ(dl);
    const BoutReal isr = SQ(dr);
    const BoutReal isc = (13.0 / 3.0) * SQ(lhs(ixp) - 2.0 * lhs(idx) + lhs(ixm))
                         + 0.25 * SQ(lhs(ixp) - lhs(ixm));

    const BoutReal al = 0.25 / SQ(WENO_SMALL + isl);
    const BoutReal ar = 0.25 / SQ(WENO_SMALL + isr);
    const BoutReal ac = 0.5 / SQ(WENO_SMALL + isc);
    const BoutReal sa = al + ar + ac;

    return (al * dl + ar * dr + ac * dc) / (sa * coords.dx(idx));
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_w3_none(int idx,
                                                           const LView& lhs) const {
    const int ixp = i_xp(idx, ny, nz);
    const int ixm = i_xm(idx, ny, nz);
    const int ixp2 = i_xp(ixp, ny, nz);
    const int ixm2 = i_xm(ixm, ny, nz);

    BoutReal ma = fabs(lhs(idx));
    ma = BOUTMAX(ma, fabs(lhs(ixm)));
    ma = BOUTMAX(ma, fabs(lhs(ixp)));
    ma = BOUTMAX(ma, fabs(lhs(ixm2)));
    ma = BOUTMAX(ma, fabs(lhs(ixp2)));

    const BoutReal sp_mm = lhs(ixm2) + ma;
    const BoutReal sp_m = lhs(ixm) + ma;
    const BoutReal sp_c = lhs(idx) + ma;
    const BoutReal sp_p = lhs(ixp) + ma;
    const BoutReal sm_m = ma - lhs(ixm);
    const BoutReal sm_c = ma - lhs(idx);
    const BoutReal sm_p = ma - lhs(ixp);
    const BoutReal sm_pp = ma - lhs(ixp2);

    BoutReal r = (WENO_SMALL + SQ(sp_c - 2.0 * sp_m + sp_mm))
                 / (WENO_SMALL + SQ(sp_p - 2.0 * sp_c + sp_m));
    BoutReal deriv = -sp_mm + 3.0 * sp_m - 3.0 * sp_c + sp_p;
    const BoutReal w_pos = 1.0 / (1.0 + 2.0 * r * r);
    const BoutReal pos = 0.25 * ((sp_p - sp_m) - w_pos * deriv);

    r = (WENO_SMALL + SQ(sm_pp - 2.0 * sm_p + sm_c))
        / (WENO_SMALL + SQ(sm_p - 2.0 * sm_c + sm_m));
    deriv = -sm_m + 3.0 * sm_c - 3.0 * sm_p + sm_pp;
    const BoutReal w_neg = 1.0 / (1.0 + 2.0 * r * r);
    const BoutReal neg = -0.25 * ((sm_p - sm_m) - w_neg * deriv);

    return (pos + neg) / coords.dx(idx);
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_s2_none(int idx,
                                                           const LView& lhs) const {
    const int ixp = i_xp(idx, ny, nz);
    const int ixm = i_xm(idx, ny, nz);
    const int ixp2 = i_xp(ixp, ny, nz);
    const int ixm2 = i_xm(ixm, ny, nz);

    BoutReal result = (-lhs(ixp2) + 8.0 * lhs(ixp) - 8.0 * lhs(ixm) + lhs(ixm2)) / 12.0;
    result += SIGN(lhs(idx))
              * (lhs(ixp2) - 4.0 * lhs(ixp) + 6.0 * lhs(idx) - 4.0 * lhs(ixm) + lhs(ixm2))
              / 12.0;

    return result / coords.dx(idx);
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c2_c2l(int idx,
                                                          const LView& lhs) const {
    const int ixm = i_xm(idx, ny, nz);

    return (lhs(idx) - lhs(ixm)) / coords.dx(idx);
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c4_c2l(int idx,
                                                          const LView& lhs) const {
    const int ixm = i_xm(idx, ny, nz);
    const int ixp = i_xp(idx, ny, nz);
    const int ixm2 = i_xm(ixm, ny, nz);

    return (27.0 * (lhs(idx) - lhs(ixm)) - (lhs(ixp) - lhs(ixm2)))
           / (24.0 * coords.dx(idx));
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c2_l2c(int idx,
                                                          const LView& lhs) const {
    const int ixp = i_xp(idx, ny, nz);

    return (lhs(ixp) - lhs(idx)) / coords.dx(idx);
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_c4_l2c(int idx,
                                                          const LView& lhs) const {
    const int ixp = i_xp(idx, ny, nz);
    const int ixm = i_xm(idx, ny, nz);
    const int ixp2 = i_xp(ixp, ny, nz);

    return (27.0 * (lhs(ixp) - lhs(idx)) - (lhs(ixp2) - lhs(ixm)))
           / (24.0 * coords.dx(idx));
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

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_ddx(int idx, const LView& lhs) const {
    switch (method) {
    case DIFF_C2:
      return apply_c2(idx, lhs);
    case DIFF_C4:
      return apply_c4(idx, lhs);
    case DIFF_W2:
      return apply_w2_none(idx, lhs);
    case DIFF_W3:
      return apply_w3_none(idx, lhs);
    case DIFF_S2:
      return apply_s2_none(idx, lhs);
    default:
      return 0.0;
    }
  }

  template <typename LView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal apply_torsion(int idx,
                                                           const LView& lhs) const {
    if (!inc_int_shear) {
      return 0.0;
    }

    switch (method) {
    case DIFF_C2:
      return coords.IntShiftTorsion(idx) * apply_ddz_c2(idx, lhs);
    case DIFF_C4:
      return coords.IntShiftTorsion(idx) * apply_ddz_c4(idx, lhs);
    default:
      return 0.0;
    }
  }

  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& lhs,
                                                        const RView&) const {
    return apply_ddx(idx, lhs) + apply_torsion(idx, lhs);
  }
};

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
using DDXExprC2 = BinaryExpr<Field3D, Field3D, Field3D, DDX_C2_Op>;
using DDXExprC4 = BinaryExpr<Field3D, Field3D, Field3D, DDX_C4_Op>;
using DDXDispatchExpr = BinaryExpr<Field3D, Field3D, Field3D, DDX_Dispatch_Op>;
using BracketArakawaExpr = BinaryExpr<Field3D, Field3D, Field3D, BracketArakawaOp>;

} // namespace bout::stencil

namespace {

inline int requiredDDXGuards(DIFF_METHOD method, STAGGER stagger) {
  if (stagger != STAGGER::None) {
    switch (method) {
    case DIFF_C2:
      return 1;
    case DIFF_C4:
      return 2;
    default:
      return -1;
    }
  }

  switch (method) {
  case DIFF_C2:
  case DIFF_W2:
    return 1;
  case DIFF_C4:
  case DIFF_W3:
  case DIFF_S2:
    return 2;
  default:
    return -1;
  }
}

inline DIFF_METHOD parseField3DMethodString(const std::string& method) {
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

  throw BoutException("Unknown field derivative method '{:s}'", method);
}

} // namespace

inline bout::stencil::DDXExprC2 DDX_C2(const Field3D& f) {
  checkData(f);

  const auto region_id = f.getMesh()->getRegionID("RGN_NOBNDRY");

  return bout::stencil::DDXExprC2{
      static_cast<Field3D::View>(f),
      static_cast<Field3D::View>(f),
      bout::stencil::DDX_C2_Op{CoordinatesAccessor{f.getCoordinates()}, f.getNy(),
                               f.getNz(), f.getMesh()->IncIntShear},
      f.getMesh(),
      f.getLocation(),
      f.getDirections(),
      region_id,
      f.getMesh()->getRegion("RGN_NOBNDRY")};
}

inline bout::stencil::DDXExprC4 DDX_C4(const Field3D& f) {
  checkData(f);

  const auto region_id = f.getMesh()->getRegionID("RGN_NOBNDRY");

  return bout::stencil::DDXExprC4{
      static_cast<Field3D::View>(f),
      static_cast<Field3D::View>(f),
      bout::stencil::DDX_C4_Op{CoordinatesAccessor{f.getCoordinates()}, f.getNy(),
                               f.getNz(), f.getMesh()->IncIntShear},
      f.getMesh(),
      f.getLocation(),
      f.getDirections(),
      region_id,
      f.getMesh()->getRegion("RGN_NOBNDRY")};
}

inline bout::stencil::DDXDispatchExpr DDX(const Field3D& f,
                                          CELL_LOC outloc = CELL_DEFAULT,
                                          DIFF_METHOD method = DIFF_DEFAULT,
                                          const std::string& region = "RGN_NOBNDRY") {
  checkData(f);

  const auto resolved_outloc = (outloc == CELL_DEFAULT) ? f.getLocation() : outloc;
  const auto stagger =
      f.getMesh()->getStagger(f.getLocation(), resolved_outloc, CELL_XLOW);

  if (method == DIFF_DEFAULT) {
    method = f.getMesh()->getDefaultMethod(DIRECTION::X, DERIV::Standard, stagger);
  }

  const bool supported_none = (method == DIFF_C2) || (method == DIFF_C4)
                              || (method == DIFF_W2) || (method == DIFF_W3)
                              || (method == DIFF_S2);
  const bool supported_staggered = (method == DIFF_C2) || (method == DIFF_C4);

  if ((stagger == STAGGER::None) ? !supported_none : !supported_staggered) {
    throw BoutException("DDX only supports DIFF_C2, DIFF_C4, DIFF_W2, DIFF_W3, "
                        "and DIFF_S2 for unstaggered grids, and DIFF_C2/DIFF_C4 "
                        "for staggered grids; got {:s}",
                        toString(method));
  }

  if (f.getMesh()->IncIntShear && (method != DIFF_C2) && (method != DIFF_C4)) {
    throw BoutException("DDX with integrated shear only supports DIFF_C2 and "
                        "DIFF_C4, got {:s}",
                        toString(method));
  }

  const auto required_guards = requiredDDXGuards(method, stagger);
  ASSERT2(required_guards >= 0);
  ASSERT2(f.getMesh()->getNguard(DIRECTION::X) >= required_guards);

  const auto region_id = f.getMesh()->getRegionID(region);

  return bout::stencil::DDXDispatchExpr{
      static_cast<Field3D::View>(f),
      static_cast<Field3D::View>(f),
      bout::stencil::DDX_Dispatch_Op{
          CoordinatesAccessor{f.getCoordinates(resolved_outloc)}, f.getNy(), f.getNz(),
          stagger, method, f.getMesh()->IncIntShear},
      f.getMesh(),
      resolved_outloc,
      f.getDirections(),
      region_id,
      f.getMesh()->getRegion(region)};
}

inline bout::stencil::DDXDispatchExpr DDX(const Field3D& f, CELL_LOC outloc,
                                          const std::string& method,
                                          const std::string& region = "RGN_NOBNDRY") {
  return DDX(f, outloc, parseField3DMethodString(method), region);
}

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

inline bout::stencil::DDZDispatchExpr
DDZ_stencil(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
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
    throw BoutException("DDZ_stencil only supports DIFF_C2 and DIFF_C4, got {:s}",
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

inline bout::stencil::DDZDispatchExpr
DDZ_stencil(const Field3D& f, CELL_LOC outloc, const std::string& method,
            const std::string& region = "RGN_NOBNDRY") {
  return DDZ_stencil(f, outloc, parseField3DMethodString(method), region);
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
