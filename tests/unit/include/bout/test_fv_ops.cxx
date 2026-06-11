#include "gtest/gtest.h"

#include <bout/fv_ops.hxx>

TEST(FVOpsLimiterTest, VanAlbadaConstant) {
  FV::Stencil1D s{};
  s.m = 1.0;
  s.c = 1.0;
  s.p = 1.0;

  FV::VanAlbada limiter;
  limiter(s);

  EXPECT_DOUBLE_EQ(s.L, 1.0);
  EXPECT_DOUBLE_EQ(s.R, 1.0);
}

TEST(FVOpsLimiterTest, VanAlbadaLinear) {
  FV::Stencil1D s{};
  s.m = 0.0;
  s.c = 1.0;
  s.p = 2.0;

  FV::VanAlbada limiter;
  limiter(s);

  EXPECT_NEAR(s.L, 0.5, 1e-14);
  EXPECT_NEAR(s.R, 1.5, 1e-14);
}

TEST(FVOpsLimiterTest, VanAlbadaExtremumGivesZeroSlope) {
  FV::Stencil1D s{};
  s.m = 0.0;
  s.c = 1.0;
  s.p = 0.0;

  FV::VanAlbada limiter;
  limiter(s);

  EXPECT_NEAR(s.L, 1.0, 1e-14);
  EXPECT_NEAR(s.R, 1.0, 1e-14);
}

TEST(FVOpsLimiterTest, VanAlbadaLimitsToSmallerGradient) {
  FV::Stencil1D s{};
  s.m = 0.0;
  s.c = 1.0;
  s.p = 1.5;

  FV::VanAlbada limiter;
  limiter(s);

  // dl = 1, dr = 0.5 => slope = (dl*dr*(dl+dr)) / (dl^2+dr^2) = 0.6
  EXPECT_NEAR(s.L, 0.7, 1e-12);
  EXPECT_NEAR(s.R, 1.3, 1e-12);
}

TEST(FVOpsLimiterTest, WENO3Constant) {
  FV::Stencil1D s{};
  s.m = 1.0;
  s.c = 1.0;
  s.p = 1.0;

  FV::WENO3 recon;
  recon(s);

  EXPECT_DOUBLE_EQ(s.L, 1.0);
  EXPECT_DOUBLE_EQ(s.R, 1.0);
}

TEST(FVOpsLimiterTest, WENO3Linear) {
  FV::Stencil1D s{};
  s.m = 0.0;
  s.c = 1.0;
  s.p = 2.0;

  FV::WENO3 recon;
  recon(s);

  EXPECT_NEAR(s.L, 0.5, 1e-14);
  EXPECT_NEAR(s.R, 1.5, 1e-14);
}
