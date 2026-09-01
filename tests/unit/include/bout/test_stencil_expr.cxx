#include "gtest/gtest.h"

#include "test_extras.hxx"

#include "bout/derivs.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/stencil_expr.hxx"

#include <string>

#include "fake_mesh.hxx"
#include "fake_mesh_fixture.hxx"

namespace {

Field3D makeTestField(Mesh* mesh, CELL_LOC location) {
  Field3D result(mesh, location);
  result.allocate();

  for (auto i : result.getRegion("RGN_ALL")) {
    result[i] = 0.1 * i.x() + 0.2 * i.y() + 0.3 * i.z() + 0.05 * i.x() * i.z();
  }

  return result;
}

void fillTestField(Field3D& result) {
  for (auto i : result.getRegion("RGN_ALL")) {
    result[i] = 0.1 * i.x() + 0.2 * i.y() + 0.3 * i.z() + 0.05 * i.x() * i.z();
  }
}

Field3DParallel makeParallelTestField(Mesh* mesh, CELL_LOC location) {
  Field3DParallel result(mesh, location);
  result.allocate();
  fillTestField(result);
  result.splitParallelSlices();

  for (size_t slice = 0; slice < result.numberParallelSlices(); ++slice) {
    result.yup(slice).allocate();
    result.ydown(slice).allocate();
    fillTestField(result.yup(slice));
    fillTestField(result.ydown(slice));
  }

  return result;
}

Field3D expectedDDX(const Field3D& input, CELL_LOC outloc, const std::string& method,
                    const std::string& region = "RGN_NOBNDRY") {
  const auto resolved_outloc = (outloc == CELL_DEFAULT) ? input.getLocation() : outloc;
  const auto* coords = input.getCoordinates(resolved_outloc);
  auto dx = Field3D{coords->dx()}.setLocation(resolved_outloc);
  Field3D result = bout::derivatives::index::DDX(input, resolved_outloc, method, region);
  result /= dx;

  if (input.getMesh()->IncIntShear) {
    auto dz = Field3D{coords->dz()}.setLocation(resolved_outloc);
    auto torsion = Field3D{coords->IntShiftTorsion()}.setLocation(resolved_outloc);
    auto torsion_term =
        bout::derivatives::index::DDZ(input, resolved_outloc, method, region);
    torsion_term /= dz;
    torsion_term *= torsion;
    result += torsion_term;
  }

  return result;
}

class DDXDispatchExprTest : public FakeMeshFixture_tmpl<7, 5, 7> {
public:
  DDXDispatchExprTest() {
    for (auto* current_mesh : {bout::globals::mesh, mesh_staggered}) {
      current_mesh->xstart = 2;
      current_mesh->xend = current_mesh->LocalNx - 3;
      current_mesh->addRegion3D(x_safe_region, Region<Ind3D>(2, current_mesh->LocalNx - 3,
                                                             0, current_mesh->LocalNy - 1,
                                                             0, current_mesh->LocalNz - 1,
                                                             current_mesh->LocalNy,
                                                             current_mesh->LocalNz));
    }
  }

  static constexpr auto x_safe_region = "RGN_XSAFE";
};

class DDXDispatchExprParamTest
    : public DDXDispatchExprTest,
      public ::testing::WithParamInterface<std::tuple<CELL_LOC, CELL_LOC, DIFF_METHOD>> {
};

class DDXDispatchExprMethodTest : public DDXDispatchExprTest,
                                  public ::testing::WithParamInterface<DIFF_METHOD> {};

class DDZDispatchExprTest : public FakeMeshFixture {};

class DDZDispatchExprParamTest
    : public DDZDispatchExprTest,
      public ::testing::WithParamInterface<std::tuple<CELL_LOC, CELL_LOC, DIFF_METHOD>> {
};

class DDYDispatchExprTest : public FakeMeshFixture {};

class DDYDispatchExprTwoSliceTest : public FakeMeshFixture_tmpl<5, 7, 7> {
public:
  DDYDispatchExprTwoSliceTest() {
    for (auto* current_mesh : {bout::globals::mesh, mesh_staggered}) {
      current_mesh->ystart = 2;
      current_mesh->yend = current_mesh->LocalNy - 3;
      current_mesh->addRegion3D(y_safe_region, Region<Ind3D>(0, current_mesh->LocalNx - 1,
                                                             2, current_mesh->LocalNy - 3,
                                                             0, current_mesh->LocalNz - 1,
                                                             current_mesh->LocalNy,
                                                             current_mesh->LocalNz));
    }
  }

  static constexpr auto y_safe_region = "RGN_YSAFE";
};

std::string paramToString(
    const ::testing::TestParamInfo<std::tuple<CELL_LOC, CELL_LOC, DIFF_METHOD>>& param) {
  const auto [inloc, outloc, method] = param.param;
  return toString(inloc) + "_to_" + toString(outloc) + "_" + toString(method);
}

std::string methodToString(const ::testing::TestParamInfo<DIFF_METHOD>& param) {
  return toString(param.param);
}

} // namespace

INSTANTIATE_TEST_SUITE_P(
    SupportedLocations, DDXDispatchExprParamTest,
    ::testing::Values(std::make_tuple(CELL_CENTRE, CELL_CENTRE, DIFF_C2),
                      std::make_tuple(CELL_CENTRE, CELL_CENTRE, DIFF_C4),
                      std::make_tuple(CELL_CENTRE, CELL_XLOW, DIFF_C2),
                      std::make_tuple(CELL_CENTRE, CELL_XLOW, DIFF_C4),
                      std::make_tuple(CELL_XLOW, CELL_CENTRE, DIFF_C2),
                      std::make_tuple(CELL_XLOW, CELL_CENTRE, DIFF_C4),
                      std::make_tuple(CELL_XLOW, CELL_XLOW, DIFF_C2),
                      std::make_tuple(CELL_XLOW, CELL_XLOW, DIFF_C4)),
    paramToString);

INSTANTIATE_TEST_SUITE_P(SupportedMethods, DDXDispatchExprMethodTest,
                         ::testing::Values(DIFF_W2, DIFF_W3, DIFF_S2), methodToString);

INSTANTIATE_TEST_SUITE_P(
    SupportedLocations, DDZDispatchExprParamTest,
    ::testing::Values(std::make_tuple(CELL_CENTRE, CELL_CENTRE, DIFF_C2),
                      std::make_tuple(CELL_CENTRE, CELL_CENTRE, DIFF_C4),
                      std::make_tuple(CELL_CENTRE, CELL_ZLOW, DIFF_C2),
                      std::make_tuple(CELL_CENTRE, CELL_ZLOW, DIFF_C4),
                      std::make_tuple(CELL_ZLOW, CELL_CENTRE, DIFF_C2),
                      std::make_tuple(CELL_ZLOW, CELL_CENTRE, DIFF_C4),
                      std::make_tuple(CELL_ZLOW, CELL_ZLOW, DIFF_C2),
                      std::make_tuple(CELL_ZLOW, CELL_ZLOW, DIFF_C4)),
    paramToString);

TEST_P(DDXDispatchExprParamTest, MatchesDDX) {
  const auto [inloc, outloc, method] = GetParam();
  auto input = makeTestField(mesh_staggered, inloc);

  const auto actual =
      Field3D{DDX(input, outloc, method, DDXDispatchExprTest::x_safe_region)};
  const auto expected =
      expectedDDX(input, outloc, toString(method), DDXDispatchExprTest::x_safe_region);

  EXPECT_EQ(actual.getLocation(), outloc);
  EXPECT_TRUE(IsFieldEqual(actual, expected, DDXDispatchExprTest::x_safe_region));
}

TEST_P(DDXDispatchExprMethodTest, MatchesDDX) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  const auto actual =
      Field3D{DDX(input, CELL_CENTRE, GetParam(), DDXDispatchExprTest::x_safe_region)};
  const auto expected = expectedDDX(input, CELL_CENTRE, toString(GetParam()),
                                    DDXDispatchExprTest::x_safe_region);

  EXPECT_EQ(actual.getLocation(), CELL_CENTRE);
  EXPECT_TRUE(IsFieldEqual(actual, expected, DDXDispatchExprTest::x_safe_region));
}

TEST_P(DDZDispatchExprParamTest, MatchesDDZ) {
  const auto [inloc, outloc, method] = GetParam();
  auto input = makeTestField(mesh_staggered, inloc);

  const auto actual = Field3D{DDZ_stencil(input, outloc, method)};
  const auto expected = DDZ(input, outloc, toString(method));

  EXPECT_EQ(actual.getLocation(), outloc);
  EXPECT_TRUE(IsFieldEqual(actual, expected, "RGN_NOBNDRY"));
}

TEST_F(DDXDispatchExprTest, UsesRequestedRegion) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  const auto actual = Field3D{DDX(input, CELL_XLOW, DIFF_C2, x_safe_region)};
  const auto expected = expectedDDX(input, CELL_XLOW, "C2", x_safe_region);

  EXPECT_EQ(actual.getLocation(), CELL_XLOW);
  EXPECT_TRUE(IsFieldEqual(actual, expected, x_safe_region));
}

TEST_F(DDXDispatchExprTest, RejectsUnsupportedMethods) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  EXPECT_THROW((void)DDX(input, CELL_DEFAULT, DIFF_U1), BoutException);
}

TEST_F(DDXDispatchExprTest, RejectsUnsupportedShearMethods) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  mesh_staggered->IncIntShear = true;

  EXPECT_THROW((void)DDX(input, CELL_DEFAULT, DIFF_W2), BoutException);
}

TEST_F(DDXDispatchExprTest, ResolvesCentredDefaultMethodFromMesh) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  Options diff_options{{"ddx", {{"first", "C4"}}}, {"ddxstag", {{"first", "C2"}}}};
  static_cast<FakeMesh*>(mesh_staggered)->initDerivs(&diff_options);

  EXPECT_EQ(mesh_staggered->getDefaultMethod(DIRECTION::X, DERIV::Standard), DIFF_C4);
  const auto centre_actual =
      Field3D{DDX(input, CELL_CENTRE, DIFF_DEFAULT, x_safe_region)};
  const auto centre_expected = expectedDDX(input, CELL_CENTRE, "C4", x_safe_region);
  EXPECT_TRUE(IsFieldEqual(centre_actual, centre_expected, x_safe_region));
}

TEST_F(DDXDispatchExprTest, ResolvesStaggeredDefaultMethodFromMesh) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  Options diff_options{{"ddx", {{"first", "C4"}}}, {"ddxstag", {{"first", "C2"}}}};
  static_cast<FakeMesh*>(mesh_staggered)->initDerivs(&diff_options);

  EXPECT_EQ(mesh_staggered->getDefaultMethod(DIRECTION::X, DERIV::Standard, STAGGER::C2L),
            DIFF_C2);
  const auto staggered_actual =
      Field3D{DDX(input, CELL_XLOW, DIFF_DEFAULT, x_safe_region)};
  const auto staggered_expected = expectedDDX(input, CELL_XLOW, "C2", x_safe_region);
  EXPECT_TRUE(IsFieldEqual(staggered_actual, staggered_expected, x_safe_region));
}

TEST_F(DDXDispatchExprTest, IncludesIntegratedShearCorrection) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  mesh_staggered->IncIntShear = true;
  auto torsion = bout::FieldMetric(0.125, mesh_staggered);
  mesh_staggered->getCoordinates()->setIntShiftTorsion(torsion);
  CoordinatesAccessor::clear(mesh_staggered->getCoordinates());

  const auto actual = Field3D{DDX(input, CELL_CENTRE, DIFF_C2, x_safe_region)};
  const auto expected = expectedDDX(input, CELL_CENTRE, "C2", x_safe_region);

  EXPECT_TRUE(IsFieldEqual(actual, expected, x_safe_region));
}

TEST_F(DDZDispatchExprTest, UsesRequestedRegion) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  const auto actual = Field3D{DDZ_stencil(input, CELL_ZLOW, DIFF_C2, "RGN_ALL")};
  const auto expected = DDZ(input, CELL_ZLOW, "C2", "RGN_ALL");

  EXPECT_EQ(actual.getLocation(), CELL_ZLOW);
  EXPECT_TRUE(IsFieldEqual(actual, expected, "RGN_ALL"));
}

TEST_F(DDZDispatchExprTest, RejectsUnsupportedMethods) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  EXPECT_THROW((void)DDZ_stencil(input, CELL_DEFAULT, DIFF_FFT), BoutException);
}

TEST_F(DDZDispatchExprTest, ResolvesDefaultMethodFromMesh) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  Options diff_options{{"ddz", {{"first", "C4"}}}, {"ddzstag", {{"first", "C2"}}}};
  static_cast<FakeMesh*>(mesh_staggered)->initDerivs(&diff_options);

  const auto centre_actual = Field3D{DDZ_stencil(input, CELL_CENTRE, DIFF_DEFAULT)};
  const auto centre_expected = DDZ(input, CELL_CENTRE, "C4");
  EXPECT_TRUE(IsFieldEqual(centre_actual, centre_expected, "RGN_NOBNDRY"));

  const auto staggered_actual = Field3D{DDZ_stencil(input, CELL_ZLOW, DIFF_DEFAULT)};
  const auto staggered_expected = DDZ(input, CELL_ZLOW, "C2");
  EXPECT_TRUE(IsFieldEqual(staggered_actual, staggered_expected, "RGN_NOBNDRY"));
}

TEST_F(DDYDispatchExprTest, MatchesDDYC2) {
  auto input = makeParallelTestField(mesh_staggered, CELL_CENTRE);

  const auto actual = Field3D{DDY_stencil(input, CELL_CENTRE, DIFF_C2)};
  const auto expected = DDY(input, CELL_CENTRE, "C2");

  EXPECT_EQ(actual.getLocation(), CELL_CENTRE);
  EXPECT_TRUE(IsFieldEqual(actual, expected, "RGN_NOBNDRY"));
}

TEST_F(DDYDispatchExprTest, RejectsWithoutParallelSlices) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  EXPECT_THROW((void)DDY_stencil(input, CELL_CENTRE, DIFF_C2), BoutException);
}

TEST_F(DDYDispatchExprTest, RejectsUnsupportedMethods) {
  auto input = makeParallelTestField(mesh_staggered, CELL_CENTRE);

  EXPECT_THROW((void)DDY_stencil(input, CELL_CENTRE, DIFF_W3), BoutException);
}

TEST_F(DDYDispatchExprTest, RejectsStaggeredOutputs) {
  auto input = makeParallelTestField(mesh_staggered, CELL_CENTRE);

  EXPECT_THROW((void)DDY_stencil(input, CELL_YLOW, DIFF_C2), BoutException);
}

TEST_F(DDYDispatchExprTest, RejectsC4WithOneSlicePair) {
  auto input = makeParallelTestField(mesh_staggered, CELL_CENTRE);

  EXPECT_THROW((void)DDY_stencil(input, CELL_CENTRE, DIFF_C4), BoutException);
}

TEST_F(DDYDispatchExprTest, ResolvesDefaultMethodFromMesh) {
  auto input = makeParallelTestField(mesh_staggered, CELL_CENTRE);

  Options diff_options{{"ddy", {{"first", "C2"}}}};
  static_cast<FakeMesh*>(mesh_staggered)->initDerivs(&diff_options);

  const auto actual = Field3D{DDY_stencil(input, CELL_CENTRE, DIFF_DEFAULT)};
  const auto expected = DDY(input, CELL_CENTRE, "C2");
  EXPECT_TRUE(IsFieldEqual(actual, expected, "RGN_NOBNDRY"));
}

TEST_F(DDYDispatchExprTwoSliceTest, MatchesDDYC4WithTwoSlicePairs) {
  auto input = makeParallelTestField(mesh_staggered, CELL_CENTRE);

  const auto actual = Field3D{DDY_stencil(input, CELL_CENTRE, DIFF_C4, y_safe_region)};
  const auto expected = DDY(input, CELL_CENTRE, "C4", y_safe_region);

  EXPECT_EQ(actual.getLocation(), CELL_CENTRE);
  EXPECT_TRUE(IsFieldEqual(actual, expected, y_safe_region));
}

TEST_F(DDYDispatchExprTwoSliceTest, ResolvesDefaultMethodFromMesh) {
  auto input = makeParallelTestField(mesh_staggered, CELL_CENTRE);

  Options diff_options{{"ddy", {{"first", "C4"}}}};
  static_cast<FakeMesh*>(mesh_staggered)->initDerivs(&diff_options);

  const auto actual =
      Field3D{DDY_stencil(input, CELL_CENTRE, DIFF_DEFAULT, y_safe_region)};
  const auto expected = DDY(input, CELL_CENTRE, "C4", y_safe_region);
  EXPECT_TRUE(IsFieldEqual(actual, expected, y_safe_region));
}
