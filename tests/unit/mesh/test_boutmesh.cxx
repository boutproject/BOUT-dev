#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "options.hxx"
#include "output.hxx"
#include "bout/griddata.hxx"

#include "test_extras.hxx"

#include <ostream>

TEST(BoutMeshTest, NullOptionsCheck) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
}

struct DecompositionTestParameters {
  int total_processors;
  int num_y_processors;
  int ny;
  int num_y_guards;
  int jyseps1_1;
  int jyseps2_1;
  int jyseps1_2;
  int jyseps2_2;
  int ny_inner;
  std::string expected_message; // Expect this fragment to be in the result.reason for bad
                                // decompositions
  std::string name;
};

std::ostream& operator<<(std::ostream& out, const DecompositionTestParameters& value) {
  return out << fmt::format("DecompositionTestParameters{{"
                            "total_processors = {}, "
                            "num_y_processors = {}, "
                            "ny = {}, "
                            "num_y_guards = {}, "
                            "jyseps1_1 = {}, "
                            "jyseps2_1 = {}, "
                            "jyseps1_2 = {}, "
                            "jyseps2_2 = {}, "
                            "ny_inner = {}, "
                            "expected_message = {} }}",
                            value.total_processors, value.num_y_processors, value.ny,
                            value.num_y_guards, value.jyseps1_1, value.jyseps2_1,
                            value.jyseps1_2, value.jyseps2_2, value.ny_inner,
                            value.expected_message);
}

std::string DecompositionTestParametersToString(
    const ::testing::TestParamInfo<DecompositionTestParameters>& param) {
  return param.param.name;
}

struct BoutMeshDecompositionTest
    : public testing::TestWithParam<DecompositionTestParameters> {
  virtual ~BoutMeshDecompositionTest() = default;
};

INSTANTIATE_TEST_SUITE_P(
    GoodDecompositions, BoutMeshDecompositionTest,
    ::testing::Values(
      DecompositionTestParameters{1, 1, 1, 1, -1, 0, 0, 0, 0, "", "OnePoint"},
      DecompositionTestParameters{1, 1, 8, 1, -1, 4, 4, 7, 4, "", "EightPoints"},
      DecompositionTestParameters{2, 1, 8, 1, -1, 4, 4, 7, 4, "", "EightPointsTwoCores"},
      DecompositionTestParameters{2, 2, 8, 1, -1, 4, 4, 7, 4, "", "EightPointsTwoCoresNYPE2"}
      ),
    DecompositionTestParametersToString);

TEST_P(BoutMeshDecompositionTest, SingleCoreYDecomposition) {
  const auto params = GetParam();
  auto result = bout::checkBoutMeshYDecomposition(
      params.total_processors, params.num_y_processors, params.ny, params.num_y_guards,
      params.jyseps1_1, params.jyseps2_1, params.jyseps1_2, params.jyseps2_2,
      params.ny_inner);

  EXPECT_TRUE(result.success);
  EXPECT_TRUE(result.reason.empty());
}

using BadBoutMeshDecompositionTest = BoutMeshDecompositionTest;

INSTANTIATE_TEST_SUITE_P(
    BasicBad, BadBoutMeshDecompositionTest,
    ::testing::Values(
      DecompositionTestParameters{2, 2, 1, 1, -1, 0, 0, 2, 0, "ny/NYPE", "TooManyCores"},
      DecompositionTestParameters{1, 1, 2, 1, 0, 0, 0, 2, 0, "Leg region jyseps1_1+1", "BadLegRegion"}),
    DecompositionTestParametersToString);

INSTANTIATE_TEST_SUITE_P(
    BadDoubleNull, BadBoutMeshDecompositionTest,
    ::testing::Values(
        DecompositionTestParameters{1, 1, 4, 1, 3, 5, 6, 10, 0,
                                    "Core region jyseps2_1-jyseps1_1", "CoreRegion1"},
        DecompositionTestParameters{1, 1, 4, 1, 3, 7, 8, 11, 0,
                                    "Core region jyseps2_2-jyseps1_2", "CoreRegion2"},
        DecompositionTestParameters{1, 1, 4, 1, 3, 7, 8, 12, 11,
                                    "leg region ny_inner-jyseps2_1-1", "UpperLeg1"},
        DecompositionTestParameters{1, 1, 4, 1, 3, 7, 8, 12, 8,
                                    "leg region jyseps1_2-ny_inner+1", "UpperLeg2"},
        DecompositionTestParameters{1, 6, 25, 1, 3, 7, 15, 19, 12,
                                    "leg region ny-jyseps2_2-1", "LegRegion"}),
    DecompositionTestParametersToString);

INSTANTIATE_TEST_SUITE_P(
    BadSingleNull, BadBoutMeshDecompositionTest,
    ::testing::Values(DecompositionTestParameters{1, 1, 4, 1, 3, 4, 4, 6, 0,
                                                  "Core region jyseps2_2-jyseps1_1",
                                                  "CoreRegion"},
                      DecompositionTestParameters{1, 3, 13, 1, 3, 4, 4, 7, 0,
                                                  "leg region ny-jyseps2_2-1",
                                                  "LegRegion"}),
    DecompositionTestParametersToString);

TEST_P(BadBoutMeshDecompositionTest, BadSingleCoreYDecomposition) {
  const auto params = GetParam();
  auto result = bout::checkBoutMeshYDecomposition(
      params.total_processors, params.num_y_processors, params.ny, params.num_y_guards,
      params.jyseps1_1, params.jyseps2_1, params.jyseps1_2, params.jyseps2_2,
      params.ny_inner);

  using ::testing::HasSubstr;

  EXPECT_FALSE(result.success);
  EXPECT_THAT(result.reason, HasSubstr(params.expected_message));
}
