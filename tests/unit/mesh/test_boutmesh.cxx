#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "options.hxx"
#include "output.hxx"
#include "bout/griddata.hxx"

#include "test_extras.hxx"

#include <ostream>

class BoutMeshExposer : public BoutMesh {
public:
  BoutMeshExposer(int input_nx, int input_ny, int input_nz, int mxg, int myg,
                  int input_npes = 1)
      : BoutMesh(input_nx, input_ny, input_nz, mxg, myg, input_npes) {}
  // Make protected methods public for testing
  using BoutMesh::chooseProcessorSplit;
  using BoutMesh::DecompositionIndices;
  using BoutMesh::findProcessorSplit;
  using BoutMesh::setYDecompositionIndices;
};

bool operator==(const BoutMeshExposer::DecompositionIndices& lhs,
                const BoutMeshExposer::DecompositionIndices& rhs) {
  return (lhs.jyseps1_1 == rhs.jyseps1_1) and (lhs.jyseps2_1 == rhs.jyseps2_1)
         and (lhs.jyseps1_2 == rhs.jyseps1_2) and (lhs.jyseps2_2 == rhs.jyseps2_2)
         and (lhs.ny_inner == rhs.ny_inner);
}

std::ostream& operator<<(std::ostream& out,
                         const BoutMeshExposer::DecompositionIndices& value) {
  return out << fmt::format("BoutMesh::DecompositionIndices{{"
                            "jyseps1_1 = {}, "
                            "jyseps2_1 = {}, "
                            "jyseps1_2 = {}, "
                            "jyseps2_2 = {}, "
                            "ny_inner = {}"
                            "}}",
                            value.jyseps1_1, value.jyseps2_1, value.jyseps1_2,
                            value.jyseps2_2, value.ny_inner);
}
TEST(BoutMeshTest, NullOptionsCheck) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
}

TEST(BoutMeshTest, SetYDecompositionIndicesCoreOnly) {
  WithQuietOutput warn{output_warn};
  const BoutMeshExposer::DecompositionIndices expected{-1, 7, 15, 23, 12};

  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  // Should return set indices unchanged
  const auto actual_indices = mesh.setYDecompositionIndices(expected);
  EXPECT_EQ(actual_indices, expected);
  EXPECT_EQ(mesh.numberOfXPoints, 0);
}

TEST(BoutMeshTest, SetYDecompositionIndicesSingleNull) {
  WithQuietOutput warn{output_warn};
  const BoutMeshExposer::DecompositionIndices expected{3, 7, 7, 19, 12};

  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  // Should return set indices unchanged
  const auto actual_indices = mesh.setYDecompositionIndices(expected);
  EXPECT_EQ(actual_indices, expected);
  EXPECT_EQ(mesh.numberOfXPoints, 1);
}

TEST(BoutMeshTest, SetYDecompositionIndicesDoubleNull) {
  WithQuietOutput warn{output_warn};
  const BoutMeshExposer::DecompositionIndices expected{3, 7, 15, 19, 12};

  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  // Should return set indices unchanged
  const auto actual_indices = mesh.setYDecompositionIndices(expected);
  EXPECT_EQ(actual_indices, expected);
  EXPECT_EQ(mesh.numberOfXPoints, 2);
}

TEST(BoutMeshTest, SetYDecompositionIndicesJyseps11Low) {
  WithQuietOutput warn{output_warn};
  const BoutMeshExposer::DecompositionIndices expected{-1, 7, 15, 19, 12};

  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  const auto actual_indices = mesh.setYDecompositionIndices({-12, 7, 15, 19, 12});
  EXPECT_EQ(actual_indices, expected);
}

TEST(BoutMeshTest, SetYDecompositionIndicesJyseps21Low) {
  WithQuietOutput warn{output_warn};
  const BoutMeshExposer::DecompositionIndices expected{3, 4, 15, 19, 12};

  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  const auto actual_indices = mesh.setYDecompositionIndices({3, 1, 15, 19, 12});
  EXPECT_EQ(actual_indices, expected);
}

TEST(BoutMeshTest, SetYDecompositionIndicesJyseps12Low) {
  WithQuietOutput warn{output_warn};
  const BoutMeshExposer::DecompositionIndices expected{3, 7, 7, 19, 12};

  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  const auto actual_indices = mesh.setYDecompositionIndices({3, 7, 5, 19, 12});
  EXPECT_EQ(actual_indices, expected);
}

TEST(BoutMeshTest, SetYDecompositionIndicesJyseps22High) {
  WithQuietOutput warn{output_warn};
  const BoutMeshExposer::DecompositionIndices expected{3, 7, 15, 23, 12};

  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  const auto actual_indices = mesh.setYDecompositionIndices({3, 7, 15, 32, 12});
  EXPECT_EQ(actual_indices, expected);
}

TEST(BoutMeshTest, SetYDecompositionIndicesJyseps22Low) {
  WithQuietOutput warn{output_warn};
  const BoutMeshExposer::DecompositionIndices expected{3, 7, 15, 15, 12};

  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  const auto actual_indices = mesh.setYDecompositionIndices({3, 7, 15, 8, 12});
  EXPECT_EQ(actual_indices, expected);
}

TEST(BoutMeshTest, SetYDecompositionIndicesJyseps22LowInconsistent) {
  WithQuietOutput warn{output_warn};
  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  EXPECT_THROW(mesh.setYDecompositionIndices({3, 7, 32, 8, 12}), BoutException);
}

struct DecompositionTestParameters {
  int total_processors;
  int num_y_processors;
  int ny;
  int num_y_guards;
  BoutMeshExposer::DecompositionIndices indices;
  std::string expected_message; // Expect this fragment to be in the result.reason for bad
                                // decompositions
  std::string name;
};

std::ostream& operator<<(std::ostream& out, const DecompositionTestParameters& value) {
  return out << fmt::format(
             "DecompositionTestParameters{{"
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
             value.total_processors, value.num_y_processors, value.ny, value.num_y_guards,
             value.indices.jyseps1_1, value.indices.jyseps2_1, value.indices.jyseps1_2,
             value.indices.jyseps2_2, value.indices.ny_inner, value.expected_message);
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
        DecompositionTestParameters{1, 1, 1, 1, {-1, 0, 0, 0, 0}, "", "OnePoint"},
        DecompositionTestParameters{1, 1, 8, 1, {-1, 4, 4, 7, 4}, "", "EightPoints"},
        DecompositionTestParameters{
            2, 1, 8, 1, {-1, 4, 4, 7, 4}, "", "EightPointsTwoCores"},
        DecompositionTestParameters{
            2, 2, 8, 1, {-1, 4, 4, 7, 4}, "", "EightPointsTwoCoresNYPE2"}),
    DecompositionTestParametersToString);

TEST_P(BoutMeshDecompositionTest, SingleCoreYDecomposition) {
  const auto params = GetParam();
  auto result = bout::checkBoutMeshYDecomposition(
      params.total_processors, params.num_y_processors, params.ny, params.num_y_guards,
      params.indices.jyseps1_1, params.indices.jyseps2_1, params.indices.jyseps1_2,
      params.indices.jyseps2_2, params.indices.ny_inner);

  EXPECT_TRUE(result.success);
  EXPECT_TRUE(result.reason.empty());
}

using BadBoutMeshDecompositionTest = BoutMeshDecompositionTest;

INSTANTIATE_TEST_SUITE_P(
    BasicBad, BadBoutMeshDecompositionTest,
    ::testing::Values(
        DecompositionTestParameters{
            2, 2, 1, 1, {-1, 0, 0, 2, 0}, "ny/NYPE", "TooManyCores"},
        DecompositionTestParameters{
            1, 1, 2, 1, {0, 0, 0, 2, 0}, "Leg region jyseps1_1+1", "BadLegRegion"}),
    DecompositionTestParametersToString);

INSTANTIATE_TEST_SUITE_P(
    BadDoubleNull, BadBoutMeshDecompositionTest,
    ::testing::Values(
        DecompositionTestParameters{
            1, 1, 4, 1, {3, 5, 6, 10, 0}, "Core region jyseps2_1", "CoreRegion1"},
        DecompositionTestParameters{
            1, 1, 4, 1, {3, 7, 8, 11, 0}, "Core region jyseps2_2", "CoreRegion2"},
        DecompositionTestParameters{
            1, 1, 4, 1, {3, 7, 8, 12, 11}, "leg region ny_inner", "UpperLeg1"},
        DecompositionTestParameters{
            1, 1, 4, 1, {3, 7, 8, 12, 8}, "leg region jyseps1_2-ny_inner+1", "UpperLeg2"},
        DecompositionTestParameters{
            1, 6, 25, 1, {3, 7, 15, 19, 12}, "leg region ny-jyseps2_2-1", "LegRegion"}),
    DecompositionTestParametersToString);

INSTANTIATE_TEST_SUITE_P(
    BadSingleNull, BadBoutMeshDecompositionTest,
    ::testing::Values(
        DecompositionTestParameters{
            1, 1, 4, 1, {3, 4, 4, 6, 0}, "Core region jyseps2_2-jyseps1_1", "CoreRegion"},
        DecompositionTestParameters{
            1, 3, 13, 1, {3, 4, 4, 7, 0}, "leg region ny-jyseps2_2-1", "LegRegion"}),
    DecompositionTestParametersToString);

TEST_P(BadBoutMeshDecompositionTest, BadSingleCoreYDecomposition) {
  const auto params = GetParam();
  auto result = bout::checkBoutMeshYDecomposition(
      params.total_processors, params.num_y_processors, params.ny, params.num_y_guards,
      params.indices.jyseps1_1, params.indices.jyseps2_1, params.indices.jyseps1_2,
      params.indices.jyseps2_2, params.indices.ny_inner);

  using ::testing::HasSubstr;

  EXPECT_FALSE(result.success);
  EXPECT_THAT(result.reason, HasSubstr(params.expected_message));
}

TEST(BoutMeshTest, ChooseProcessorSplitBadNXPE) {
  WithQuietOutput info{output_info};
  Options options{{"NXPE", 3}};

  BoutMeshExposer mesh(1, 24, 1, 1, 1, 8);

  EXPECT_THROW(mesh.chooseProcessorSplit(options), BoutException);
}

TEST(BoutMeshTest, ChooseProcessorSplitBadNYPE) {
  WithQuietOutput info{output_info};
  Options options{{"NYPE", 7}};

  BoutMeshExposer mesh(1, 24, 1, 1, 1, 8);

  EXPECT_THROW(mesh.chooseProcessorSplit(options), BoutException);
}

TEST(BoutMeshTest, ChooseProcessorSplitNXPE) {
  WithQuietOutput info{output_info};
  Options options{{"NXPE", 4}};

  BoutMeshExposer mesh(1, 24, 1, 1, 1, 8);

  EXPECT_NO_THROW(mesh.chooseProcessorSplit(options));

  EXPECT_EQ(mesh.getNXPE(), 4);
  EXPECT_EQ(mesh.getNYPE(), 2);
}

TEST(BoutMeshTest, ChooseProcessorSplitBadNXPENotEnoughGuards) {
  WithQuietOutput info{output_info};
  Options options{{"NXPE", 4}};

  BoutMeshExposer mesh(1, 24, 1, 1, 13, 8);

  EXPECT_THROW(mesh.chooseProcessorSplit(options), BoutException);
}

TEST(BoutMeshTest, ChooseProcessorSplitNYPE) {
  WithQuietOutput info{output_info};
  Options options{{"NYPE", 4}};

  BoutMeshExposer mesh(1, 24, 1, 1, 1, 8);

  EXPECT_NO_THROW(mesh.chooseProcessorSplit(options));

  EXPECT_EQ(mesh.getNXPE(), 2);
  EXPECT_EQ(mesh.getNYPE(), 4);
}

// TODO: Parameterise this test
TEST(BoutMeshTest, FindProcessorSplit1x24x1x8) {
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh(1, 24, 1, 0, 0, 8);

  EXPECT_NO_THROW(mesh.findProcessorSplit());

  EXPECT_EQ(mesh.getNXPE(), 1);
  EXPECT_EQ(mesh.getNYPE(), 8);
}
