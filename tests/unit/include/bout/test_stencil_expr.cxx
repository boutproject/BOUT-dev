#include "gtest/gtest.h"

#include "test_extras.hxx"

#include "bout/derivs.hxx"
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

class DDZDispatchExprTest : public FakeMeshFixture {};

class DDZDispatchExprParamTest
    : public DDZDispatchExprTest,
      public ::testing::WithParamInterface<std::tuple<CELL_LOC, CELL_LOC, DIFF_METHOD>> {
};

std::string paramToString(
    const ::testing::TestParamInfo<std::tuple<CELL_LOC, CELL_LOC, DIFF_METHOD>>& param) {
  const auto [inloc, outloc, method] = param.param;
  return toString(inloc) + "_to_" + toString(outloc) + "_" + toString(method);
}

} // namespace

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

TEST_P(DDZDispatchExprParamTest, MatchesDDZ) {
  const auto [inloc, outloc, method] = GetParam();
  auto input = makeTestField(mesh_staggered, inloc);

  const auto actual = Field3D{DDZ_Dispatch(input, outloc, method)};
  const auto expected = DDZ(input, outloc, toString(method));

  EXPECT_EQ(actual.getLocation(), outloc);
  EXPECT_TRUE(IsFieldEqual(actual, expected, "RGN_NOBNDRY"));
}

TEST_F(DDZDispatchExprTest, UsesRequestedRegion) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  const auto actual = Field3D{DDZ_Dispatch(input, CELL_ZLOW, DIFF_C2, "RGN_ALL")};
  const auto expected = DDZ(input, CELL_ZLOW, toString(DIFF_C2), "RGN_ALL");

  EXPECT_EQ(actual.getLocation(), CELL_ZLOW);
  EXPECT_TRUE(IsFieldEqual(actual, expected, "RGN_ALL"));
}

TEST_F(DDZDispatchExprTest, RejectsUnsupportedMethods) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  EXPECT_THROW((void)DDZ_Dispatch(input, CELL_DEFAULT, DIFF_FFT), BoutException);
}

TEST_F(DDZDispatchExprTest, ResolvesDefaultMethodFromMesh) {
  auto input = makeTestField(mesh_staggered, CELL_CENTRE);

  Options diff_options{{"ddz", {{"first", "C4"}}}, {"ddzstag", {{"first", "C2"}}}};
  static_cast<FakeMesh*>(mesh_staggered)->initDerivs(&diff_options);

  const auto centre_actual = Field3D{DDZ_Dispatch(input, CELL_CENTRE, DIFF_DEFAULT)};
  const auto centre_expected = DDZ(input, CELL_CENTRE, toString(DIFF_C4));
  EXPECT_TRUE(IsFieldEqual(centre_actual, centre_expected, "RGN_NOBNDRY"));

  const auto staggered_actual = Field3D{DDZ_Dispatch(input, CELL_ZLOW, DIFF_DEFAULT)};
  const auto staggered_expected = DDZ(input, CELL_ZLOW, toString(DIFF_C2));
  EXPECT_TRUE(IsFieldEqual(staggered_actual, staggered_expected, "RGN_NOBNDRY"));
}
