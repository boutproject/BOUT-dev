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
  BoutMeshExposer(int nx, int ny, int nz, int nxpe, int nype, int pe_xind, int pe_yind,
                  bool create_topology = true)
      : BoutMesh((nxpe * (nx - 2)) + 2, nype * ny, nz, 1, 1, nxpe, nype, pe_xind, pe_yind,
                 create_topology) {}
  // Make protected methods public for testing
  using BoutMesh::add_target;
  using BoutMesh::addBoundaryRegions;
  using BoutMesh::chooseProcessorSplit;
  using BoutMesh::YDecompositionIndices;
  using BoutMesh::default_connections;
  using BoutMesh::findProcessorSplit;
  using BoutMesh::getConnectionInfo;
  using BoutMesh::PROC_NUM;
  using BoutMesh::set_connection;
  using BoutMesh::setXDecompositionIndices;
  using BoutMesh::setYDecompositionIndices;
  using BoutMesh::topology;
  using BoutMesh::XPROC;
  using BoutMesh::YPROC;
};

bool operator==(const BoutMeshExposer::YDecompositionIndices& lhs,
                const BoutMeshExposer::YDecompositionIndices& rhs) {
  return (lhs.jyseps1_1 == rhs.jyseps1_1) and (lhs.jyseps2_1 == rhs.jyseps2_1)
         and (lhs.jyseps1_2 == rhs.jyseps1_2) and (lhs.jyseps2_2 == rhs.jyseps2_2)
         and (lhs.ny_inner == rhs.ny_inner);
}

std::ostream& operator<<(std::ostream& out,
                         const BoutMeshExposer::YDecompositionIndices& value) {
  return out << fmt::format("BoutMesh::YDecompositionIndices{{"
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

struct SetYDecompositionTestParameters {
  BoutMeshExposer::YDecompositionIndices input;
  BoutMeshExposer::YDecompositionIndices expected;
  int number_of_X_points;
  std::string test_name;
};

std::ostream& operator<<(std::ostream& out,
                         const SetYDecompositionTestParameters& value) {
  return out << "SetYDecompositionTestParameters{input=" << value.input
             << ", expected=" << value.expected
             << ", number_of_X_points=" << value.number_of_X_points << "}";
}

std::string SetYDecompositionTestParametersToString(
    const ::testing::TestParamInfo<SetYDecompositionTestParameters>& param) {
  return param.param.test_name;
}

struct BoutMeshSetYDecompositionTest : public ::testing::TestWithParam<SetYDecompositionTestParameters> {
  virtual ~BoutMeshSetYDecompositionTest() = default;
};

INSTANTIATE_TEST_SUITE_P(GoodDecompositions, BoutMeshSetYDecompositionTest,
                         ::testing::Values(
                           SetYDecompositionTestParameters{{-1, 7, 15, 23, 12}, {-1, 7, 15, 23, 12}, 0, "CoreOnly"},
                           SetYDecompositionTestParameters{{3, 7, 7, 19, 12}, {3, 7, 7, 19, 12}, 1, "SingleNull"},
                           SetYDecompositionTestParameters{{3, 7, 15, 19, 12}, {3, 7, 15, 19, 12}, 2, "DoubleNull"},
                           SetYDecompositionTestParameters{{-12, 7, 15, 19, 12}, {-1, 7, 15, 19, 12}, 2, "Jyseps11Low"},
                           SetYDecompositionTestParameters{{3, 1, 15, 19, 12}, {3, 4, 15, 19, 12}, 2, "Jyseps21Low"},
                           SetYDecompositionTestParameters{{3, 7, 5, 19, 12}, {3, 7, 7, 19, 12}, 1, "Jyseps12Low"},
                           SetYDecompositionTestParameters{{3, 7, 15, 32, 12}, {3, 7, 15, 23, 12}, 2, "Jyseps22High"},
                           SetYDecompositionTestParameters{{3, 7, 15, 8, 12}, {3, 7, 15, 15, 12}, 2, "Jyseps22Low"}
                           ),
                           SetYDecompositionTestParametersToString);

TEST_P(BoutMeshSetYDecompositionTest, BasicTest) {
  WithQuietOutput warn{output_warn};
  const auto params = GetParam();

  BoutMeshExposer mesh(1, 24, 1, 1, 1);
  const auto actual_indices = mesh.setYDecompositionIndices(params.input);
  EXPECT_EQ(actual_indices, params.expected);
  EXPECT_EQ(mesh.numberOfXPoints, params.number_of_X_points);
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
  BoutMeshExposer::YDecompositionIndices indices;
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

struct FindProcessorParameters {
  int total_processors;
  int nx;
  int ny;
  int num_x_guards;
  int num_y_guards;
  BoutMeshExposer::YDecompositionIndices indices;
  int expected_nxpe;
  int expected_nype;
};

std::ostream& operator<<(std::ostream& out, const FindProcessorParameters& value) {
  return out << fmt::format(
             "FindProcessorParameters{{"
             "total_processors = {}, "
             "nx = {}, "
             "ny = {}, "
             "num_x_guards = {}, "
             "num_y_guards = {}, "
             "jyseps1_1 = {}, "
             "jyseps2_1 = {}, "
             "jyseps1_2 = {}, "
             "jyseps2_2 = {}, "
             "ny_inner = {}, "
             "expected_nxpe = {},"
             "expected_nype = {} }}",
             value.total_processors, value.nx, value.ny, value.num_x_guards,
             value.num_y_guards, value.indices.jyseps1_1, value.indices.jyseps2_1,
             value.indices.jyseps1_2, value.indices.jyseps2_2, value.indices.ny_inner,
             value.expected_nxpe, value.expected_nype);
}

struct BoutMeshFindProcessorTest
    : public testing::TestWithParam<FindProcessorParameters> {
  virtual ~BoutMeshFindProcessorTest() = default;
};

INSTANTIATE_TEST_SUITE_P(
    GoodDecompositions, BoutMeshFindProcessorTest,
    ::testing::Values(
        FindProcessorParameters{8, 1, 24, 0, 0, {-1, 12, 12, 23, 12}, 1, 8},
        FindProcessorParameters{12, 1, 24, 0, 0, {-1, 12, 12, 23, 12}, 1, 12},
        FindProcessorParameters{24, 1, 24, 0, 0, {-1, 12, 12, 23, 12}, 1, 24},
        FindProcessorParameters{8, 32, 24, 0, 0, {-1, 12, 12, 23, 12}, 4, 2},
        FindProcessorParameters{12, 32, 24, 0, 0, {-1, 12, 12, 23, 12}, 4, 3},
        FindProcessorParameters{24, 32, 24, 0, 0, {-1, 12, 12, 23, 12}, 4, 6},
        FindProcessorParameters{8, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 1, 8},
        FindProcessorParameters{16, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 2, 8},
        FindProcessorParameters{32, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 4, 8},
        FindProcessorParameters{64, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 8, 8},
        FindProcessorParameters{256, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 16, 16},
        FindProcessorParameters{8192, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 128, 64}));

TEST_P(BoutMeshFindProcessorTest, FindProcessor) {
  WithQuietOutput info{output_info};
  WithQuietOutput progress{output_progress};
  WithQuietOutput warn{output_warn};

  const auto params = GetParam();

  BoutMeshExposer mesh(params.nx, params.ny, 1, params.num_x_guards, params.num_y_guards,
                       params.total_processors);

  mesh.setYDecompositionIndices(params.indices);

  EXPECT_NO_THROW(mesh.findProcessorSplit());

  EXPECT_EQ(mesh.getNXPE(), params.expected_nxpe);
  EXPECT_EQ(mesh.getNYPE(), params.expected_nype);
}

using BadBoutMeshFindProcessorTest = BoutMeshFindProcessorTest;

INSTANTIATE_TEST_SUITE_P(
    BadDecompositions, BadBoutMeshFindProcessorTest,
    ::testing::Values(
        FindProcessorParameters{9, 1, 24, 0, 0, {-1, 12, 12, 23, 12}, 1, 8},
        FindProcessorParameters{25, 1, 24, 0, 0, {-1, 12, 12, 23, 12}, 1, 24},
        FindProcessorParameters{9, 32, 24, 0, 0, {-1, 12, 12, 23, 12}, 4, 2},
        FindProcessorParameters{13, 32, 24, 0, 0, {-1, 12, 12, 23, 12}, 4, 3},
        FindProcessorParameters{23, 32, 24, 0, 0, {-1, 12, 12, 23, 12}, 4, 6},
        FindProcessorParameters{7, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 1, 8},
        FindProcessorParameters{24, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 2, 8},
        FindProcessorParameters{8192, 132, 128, 2, 4, {15, 47, 79, 111, 64}, 16, 16},
        FindProcessorParameters{16384, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 128, 64}));

TEST_P(BadBoutMeshFindProcessorTest, FindProcessor) {
  WithQuietOutput info{output_info};
  WithQuietOutput progress{output_progress};
  WithQuietOutput warn{output_warn};

  const auto params = GetParam();

  BoutMeshExposer mesh(params.nx, params.ny, 1, params.num_x_guards, params.num_y_guards,
                       params.total_processors);

  mesh.setYDecompositionIndices(params.indices);

  EXPECT_THROW(mesh.findProcessorSplit(), BoutException);
}

struct ProcNumParameters {
  int nxpe;
  int xind;
  int yind;
  int expected_result;
};

std::ostream& operator<<(std::ostream& out, const ProcNumParameters& value) {
  return out << fmt::format("NXPE = {}, processor index = ({}, {}), expected_result = {}",
                            value.nxpe, value.xind, value.yind, value.expected_result);
}

struct BoutMeshProcNumTest : public testing::TestWithParam<ProcNumParameters> {
  virtual ~BoutMeshProcNumTest() = default;
};

// Square domain with 4 processors:
//     +-+-+
//     |0|1|
//     +-+-+
//     |2|3|
//     +-+-+
INSTANTIATE_TEST_SUITE_P(
    Square, BoutMeshProcNumTest,
    ::testing::Values(ProcNumParameters{2, -8, 1, -1}, ProcNumParameters{2, 1, -8, -1},
                      ProcNumParameters{2, 0, 0, 0}, ProcNumParameters{2, 1, 0, 1},
                      ProcNumParameters{2, 0, 1, 2}, ProcNumParameters{2, 1, 1, 3},
                      ProcNumParameters{2, 2, 1, -1}, ProcNumParameters{2, 1, 2, -1}));

// Rectangular domain with 4 processors:
//     +-+-+-+-+
//     |0|1|2|3|
//     +-+-+-+-+
INSTANTIATE_TEST_SUITE_P(
    Rectangle, BoutMeshProcNumTest,
    ::testing::Values(ProcNumParameters{4, -8, 1, -1}, ProcNumParameters{4, 1, -8, -1},
                      ProcNumParameters{4, 0, 0, 0}, ProcNumParameters{4, 1, 0, 1},
                      ProcNumParameters{4, 2, 0, 2}, ProcNumParameters{4, 3, 0, 3},
                      ProcNumParameters{4, 2, 1, -1}, ProcNumParameters{4, 1, 2, -1}));

TEST_P(BoutMeshProcNumTest, ProcNum) {
  WithQuietOutput info{output_info};
  BoutMeshExposer mesh(4, 4, 1, 1, 1, 4);

  const auto params = GetParam();
  Options options{{"NXPE", params.nxpe}};
  mesh.chooseProcessorSplit(options);

  const int result = mesh.PROC_NUM(params.xind, params.yind);
  EXPECT_EQ(result, params.expected_result);
}

TEST(BoutMeshTest, YProc) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor
  BoutMeshExposer mesh(5, 3, 1, 2, 2, 0, 0);

  // YPROC is defined over the range (0, ny=6)
  EXPECT_EQ(mesh.YPROC(-4), -1);
  EXPECT_EQ(mesh.YPROC(0), 0);
  EXPECT_EQ(mesh.YPROC(1), 0);
  EXPECT_EQ(mesh.YPROC(2), 0);
  EXPECT_EQ(mesh.YPROC(3), 1);
  EXPECT_EQ(mesh.YPROC(4), 1);
  EXPECT_EQ(mesh.YPROC(5), 1);
  EXPECT_EQ(mesh.YPROC(6), -1);
  EXPECT_EQ(mesh.YPROC(7), -1);
}

TEST(BoutMeshTest, XProc) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor
  BoutMeshExposer mesh(5, 3, 1, 2, 2, 0, 0);

  EXPECT_EQ(mesh.XPROC(-4), 0);
  EXPECT_EQ(mesh.XPROC(0), 0);
  EXPECT_EQ(mesh.XPROC(1), 0);
  EXPECT_EQ(mesh.XPROC(2), 0);
  EXPECT_EQ(mesh.XPROC(3), 0);
  EXPECT_EQ(mesh.XPROC(4), 1);
  EXPECT_EQ(mesh.XPROC(5), 1);
  EXPECT_EQ(mesh.XPROC(6), 1);
  // BoutMesh::XPROC doesn't have an upper-bound, but also is only
  // used in one function which itself is only (optionally) used in
  // one example, so probably fine
}

TEST(BoutMeshTest, GetGlobalXIndex) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Boundaries are included in the global index

  // |<--  1st X-proc -->|
  //             |<--  2nd X-proc -->|
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st X-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd X-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getGlobalXIndex(0), 0);
  EXPECT_EQ(mesh00.getGlobalXIndex(1), 1);
  EXPECT_EQ(mesh00.getGlobalXIndex(2), 2);
  EXPECT_EQ(mesh00.getGlobalXIndex(3), 3);
  EXPECT_EQ(mesh00.getGlobalXIndex(4), 4);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getGlobalXIndex(0), 0);
  EXPECT_EQ(mesh01.getGlobalXIndex(1), 1);
  EXPECT_EQ(mesh01.getGlobalXIndex(2), 2);
  EXPECT_EQ(mesh01.getGlobalXIndex(3), 3);
  EXPECT_EQ(mesh01.getGlobalXIndex(4), 4);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getGlobalXIndex(0), 3);
  EXPECT_EQ(mesh10.getGlobalXIndex(1), 4);
  EXPECT_EQ(mesh10.getGlobalXIndex(2), 5);
  EXPECT_EQ(mesh10.getGlobalXIndex(3), 6);
  EXPECT_EQ(mesh10.getGlobalXIndex(4), 7);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getGlobalXIndex(0), 3);
  EXPECT_EQ(mesh11.getGlobalXIndex(1), 4);
  EXPECT_EQ(mesh11.getGlobalXIndex(2), 5);
  EXPECT_EQ(mesh11.getGlobalXIndex(3), 6);
  EXPECT_EQ(mesh11.getGlobalXIndex(4), 7);
}

TEST(BoutMeshTest, GetGlobalXIndexNoBoundaries) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Global indices start counting from the first non-boundary point

  // |<--  1st X-proc -->|
  //             |<--  2nd X-proc -->|
  // +---+---+---+---+---+---+---+---+
  // |-1*| 0 | 1 | 2 | 3 | 4 | 5 | 6*| <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st X-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd X-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getGlobalXIndexNoBoundaries(0), -1);
  EXPECT_EQ(mesh00.getGlobalXIndexNoBoundaries(1), 0);
  EXPECT_EQ(mesh00.getGlobalXIndexNoBoundaries(2), 1);
  EXPECT_EQ(mesh00.getGlobalXIndexNoBoundaries(3), 2);
  EXPECT_EQ(mesh00.getGlobalXIndexNoBoundaries(4), 3);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getGlobalXIndexNoBoundaries(0), -1);
  EXPECT_EQ(mesh01.getGlobalXIndexNoBoundaries(1), 0);
  EXPECT_EQ(mesh01.getGlobalXIndexNoBoundaries(2), 1);
  EXPECT_EQ(mesh01.getGlobalXIndexNoBoundaries(3), 2);
  EXPECT_EQ(mesh01.getGlobalXIndexNoBoundaries(4), 3);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getGlobalXIndexNoBoundaries(0), 2);
  EXPECT_EQ(mesh10.getGlobalXIndexNoBoundaries(1), 3);
  EXPECT_EQ(mesh10.getGlobalXIndexNoBoundaries(2), 4);
  EXPECT_EQ(mesh10.getGlobalXIndexNoBoundaries(3), 5);
  EXPECT_EQ(mesh10.getGlobalXIndexNoBoundaries(4), 6);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getGlobalXIndexNoBoundaries(0), 2);
  EXPECT_EQ(mesh11.getGlobalXIndexNoBoundaries(1), 3);
  EXPECT_EQ(mesh11.getGlobalXIndexNoBoundaries(2), 4);
  EXPECT_EQ(mesh11.getGlobalXIndexNoBoundaries(3), 5);
  EXPECT_EQ(mesh11.getGlobalXIndexNoBoundaries(4), 6);
}

TEST(BoutMeshTest, GetLocalXIndex) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Boundaries are included in the local index

  // |<--  1st X-proc -->|
  //             |<--  2nd X-proc -->|
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st X-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd X-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getLocalXIndex(0), 0);
  EXPECT_EQ(mesh00.getLocalXIndex(1), 1);
  EXPECT_EQ(mesh00.getLocalXIndex(2), 2);
  EXPECT_EQ(mesh00.getLocalXIndex(3), 3);
  EXPECT_EQ(mesh00.getLocalXIndex(4), 4);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getLocalXIndex(0), 0);
  EXPECT_EQ(mesh01.getLocalXIndex(1), 1);
  EXPECT_EQ(mesh01.getLocalXIndex(2), 2);
  EXPECT_EQ(mesh01.getLocalXIndex(3), 3);
  EXPECT_EQ(mesh01.getLocalXIndex(4), 4);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getLocalXIndex(3), 0);
  EXPECT_EQ(mesh10.getLocalXIndex(4), 1);
  EXPECT_EQ(mesh10.getLocalXIndex(5), 2);
  EXPECT_EQ(mesh10.getLocalXIndex(6), 3);
  EXPECT_EQ(mesh10.getLocalXIndex(7), 4);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getLocalXIndex(3), 0);
  EXPECT_EQ(mesh11.getLocalXIndex(4), 1);
  EXPECT_EQ(mesh11.getLocalXIndex(5), 2);
  EXPECT_EQ(mesh11.getLocalXIndex(6), 3);
  EXPECT_EQ(mesh11.getLocalXIndex(7), 4);
}

TEST(BoutMeshTest, GetLocalXIndexNoBoundaries) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Local indices start counting from the first non-boundary point

  // |<--  1st X-proc -->|
  //             |<--  2nd X-proc -->|
  // +---+---+---+---+---+---+---+---+
  // |-1*| 0 | 1 | 2 | 3 | 4 | 5 | 6*| <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st X-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd X-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getLocalXIndexNoBoundaries(-1), 0);
  EXPECT_EQ(mesh00.getLocalXIndexNoBoundaries(0), 1);
  EXPECT_EQ(mesh00.getLocalXIndexNoBoundaries(1), 2);
  EXPECT_EQ(mesh00.getLocalXIndexNoBoundaries(2), 3);
  EXPECT_EQ(mesh00.getLocalXIndexNoBoundaries(3), 4);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getLocalXIndexNoBoundaries(-1), 0);
  EXPECT_EQ(mesh01.getLocalXIndexNoBoundaries(0), 1);
  EXPECT_EQ(mesh01.getLocalXIndexNoBoundaries(1), 2);
  EXPECT_EQ(mesh01.getLocalXIndexNoBoundaries(2), 3);
  EXPECT_EQ(mesh01.getLocalXIndexNoBoundaries(3), 4);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getLocalXIndexNoBoundaries(2), 0);
  EXPECT_EQ(mesh10.getLocalXIndexNoBoundaries(3), 1);
  EXPECT_EQ(mesh10.getLocalXIndexNoBoundaries(4), 2);
  EXPECT_EQ(mesh10.getLocalXIndexNoBoundaries(5), 3);
  EXPECT_EQ(mesh10.getLocalXIndexNoBoundaries(6), 4);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getLocalXIndexNoBoundaries(2), 0);
  EXPECT_EQ(mesh11.getLocalXIndexNoBoundaries(3), 1);
  EXPECT_EQ(mesh11.getLocalXIndexNoBoundaries(4), 2);
  EXPECT_EQ(mesh11.getLocalXIndexNoBoundaries(5), 3);
  EXPECT_EQ(mesh11.getLocalXIndexNoBoundaries(6), 4);
}

TEST(BoutMeshTest, GetGlobalYIndexSingleNull) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Boundaries are included in the global index

  // |<--  1st Y-proc -->|
  //             |<--  2nd Y-proc -->|
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st Y-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd Y-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getGlobalYIndex(0), 0);
  EXPECT_EQ(mesh00.getGlobalYIndex(1), 1);
  EXPECT_EQ(mesh00.getGlobalYIndex(2), 2);
  EXPECT_EQ(mesh00.getGlobalYIndex(3), 3);
  EXPECT_EQ(mesh00.getGlobalYIndex(4), 4);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getGlobalYIndex(0), 3);
  EXPECT_EQ(mesh01.getGlobalYIndex(1), 4);
  EXPECT_EQ(mesh01.getGlobalYIndex(2), 5);
  EXPECT_EQ(mesh01.getGlobalYIndex(3), 6);
  EXPECT_EQ(mesh01.getGlobalYIndex(4), 7);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getGlobalYIndex(0), 0);
  EXPECT_EQ(mesh10.getGlobalYIndex(1), 1);
  EXPECT_EQ(mesh10.getGlobalYIndex(2), 2);
  EXPECT_EQ(mesh10.getGlobalYIndex(3), 3);
  EXPECT_EQ(mesh10.getGlobalYIndex(4), 4);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getGlobalYIndex(0), 3);
  EXPECT_EQ(mesh11.getGlobalYIndex(1), 4);
  EXPECT_EQ(mesh11.getGlobalYIndex(2), 5);
  EXPECT_EQ(mesh11.getGlobalYIndex(3), 6);
  EXPECT_EQ(mesh11.getGlobalYIndex(4), 7);
}

TEST(BoutMeshTest, GetGlobalYIndexDoubleNull) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Boundaries are included in the global index
  // Double-null, so extra boundary in middle of domain

  // |<--  1st Y-proc -->|
  //                     |<--  2nd Y-proc -->|
  // +---+---+---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | <- Global indices
  // +---+---+---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | <- 1st Y-processor
  // +---+---+---+---+---+---+---+---+---+---+
  // |-5 |-4 |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd Y-processor
  // +---+---+---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  mesh00.setYDecompositionIndices({-1, 2, 5, 5, 4});
  EXPECT_EQ(mesh00.getGlobalYIndex(0), 0);
  EXPECT_EQ(mesh00.getGlobalYIndex(1), 1);
  EXPECT_EQ(mesh00.getGlobalYIndex(2), 2);
  EXPECT_EQ(mesh00.getGlobalYIndex(3), 3);
  EXPECT_EQ(mesh00.getGlobalYIndex(4), 4);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  mesh01.setYDecompositionIndices({-1, 2, 5, 5, 4});
  EXPECT_EQ(mesh01.getGlobalYIndex(0), 5);
  EXPECT_EQ(mesh01.getGlobalYIndex(1), 6);
  EXPECT_EQ(mesh01.getGlobalYIndex(2), 7);
  EXPECT_EQ(mesh01.getGlobalYIndex(3), 8);
  EXPECT_EQ(mesh01.getGlobalYIndex(4), 9);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  mesh10.setYDecompositionIndices({-1, 2, 5, 5, 4});
  EXPECT_EQ(mesh10.getGlobalYIndex(0), 0);
  EXPECT_EQ(mesh10.getGlobalYIndex(1), 1);
  EXPECT_EQ(mesh10.getGlobalYIndex(2), 2);
  EXPECT_EQ(mesh10.getGlobalYIndex(3), 3);
  EXPECT_EQ(mesh10.getGlobalYIndex(4), 4);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  mesh11.setYDecompositionIndices({-1, 2, 5, 5, 4});
  EXPECT_EQ(mesh11.getGlobalYIndex(0), 5);
  EXPECT_EQ(mesh11.getGlobalYIndex(1), 6);
  EXPECT_EQ(mesh11.getGlobalYIndex(2), 7);
  EXPECT_EQ(mesh11.getGlobalYIndex(3), 8);
  EXPECT_EQ(mesh11.getGlobalYIndex(4), 9);
}

TEST(BoutMeshTest, GetGlobalYIndexNoBoundaries) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Global indices start counting from the first non-boundary point

  // |<--  1st Y-proc -->|
  //             |<--  2nd Y-proc -->|
  // +---+---+---+---+---+---+---+---+
  // |-1*| 0 | 1 | 2 | 3 | 4 | 5 | 6*| <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st Y-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd Y-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getGlobalYIndexNoBoundaries(0), -1);
  EXPECT_EQ(mesh00.getGlobalYIndexNoBoundaries(1), 0);
  EXPECT_EQ(mesh00.getGlobalYIndexNoBoundaries(2), 1);
  EXPECT_EQ(mesh00.getGlobalYIndexNoBoundaries(3), 2);
  EXPECT_EQ(mesh00.getGlobalYIndexNoBoundaries(4), 3);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getGlobalYIndexNoBoundaries(0), 2);
  EXPECT_EQ(mesh01.getGlobalYIndexNoBoundaries(1), 3);
  EXPECT_EQ(mesh01.getGlobalYIndexNoBoundaries(2), 4);
  EXPECT_EQ(mesh01.getGlobalYIndexNoBoundaries(3), 5);
  EXPECT_EQ(mesh01.getGlobalYIndexNoBoundaries(4), 6);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getGlobalYIndexNoBoundaries(0), -1);
  EXPECT_EQ(mesh10.getGlobalYIndexNoBoundaries(1), 0);
  EXPECT_EQ(mesh10.getGlobalYIndexNoBoundaries(2), 1);
  EXPECT_EQ(mesh10.getGlobalYIndexNoBoundaries(3), 2);
  EXPECT_EQ(mesh10.getGlobalYIndexNoBoundaries(4), 3);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getGlobalYIndexNoBoundaries(0), 2);
  EXPECT_EQ(mesh11.getGlobalYIndexNoBoundaries(1), 3);
  EXPECT_EQ(mesh11.getGlobalYIndexNoBoundaries(2), 4);
  EXPECT_EQ(mesh11.getGlobalYIndexNoBoundaries(3), 5);
  EXPECT_EQ(mesh11.getGlobalYIndexNoBoundaries(4), 6);
}

TEST(BoutMeshTest, GetLocalYIndexSingleNull) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Boundaries are included in the local index

  // |<--  1st Y-proc -->|
  //             |<--  2nd Y-proc -->|
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st Y-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd Y-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getLocalYIndex(0), 0);
  EXPECT_EQ(mesh00.getLocalYIndex(1), 1);
  EXPECT_EQ(mesh00.getLocalYIndex(2), 2);
  EXPECT_EQ(mesh00.getLocalYIndex(3), 3);
  EXPECT_EQ(mesh00.getLocalYIndex(4), 4);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getLocalYIndex(3), 0);
  EXPECT_EQ(mesh01.getLocalYIndex(4), 1);
  EXPECT_EQ(mesh01.getLocalYIndex(5), 2);
  EXPECT_EQ(mesh01.getLocalYIndex(6), 3);
  EXPECT_EQ(mesh01.getLocalYIndex(7), 4);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getLocalYIndex(0), 0);
  EXPECT_EQ(mesh10.getLocalYIndex(1), 1);
  EXPECT_EQ(mesh10.getLocalYIndex(2), 2);
  EXPECT_EQ(mesh10.getLocalYIndex(3), 3);
  EXPECT_EQ(mesh10.getLocalYIndex(4), 4);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getLocalYIndex(3), 0);
  EXPECT_EQ(mesh11.getLocalYIndex(4), 1);
  EXPECT_EQ(mesh11.getLocalYIndex(5), 2);
  EXPECT_EQ(mesh11.getLocalYIndex(6), 3);
  EXPECT_EQ(mesh11.getLocalYIndex(7), 4);
}

TEST(BoutMeshTest, GetLocalYIndexDoubleNull) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Boundaries are included in the global index
  // Double-null, so extra boundary in middle of domain

  // |<--  1st Y-proc -->|
  //                     |<--  2nd Y-proc -->|
  // +---+---+---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | <- Global indices
  // +---+---+---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | <- 1st Y-processor
  // +---+---+---+---+---+---+---+---+---+---+
  // |-5 |-4 |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd Y-processor
  // +---+---+---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  mesh00.setYDecompositionIndices({-1, 2, 5, 5, 4});
  EXPECT_EQ(mesh00.getLocalYIndex(0), 0);
  EXPECT_EQ(mesh00.getLocalYIndex(1), 1);
  EXPECT_EQ(mesh00.getLocalYIndex(2), 2);
  EXPECT_EQ(mesh00.getLocalYIndex(3), 3);
  EXPECT_EQ(mesh00.getLocalYIndex(4), 4);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  mesh01.setYDecompositionIndices({-1, 2, 5, 5, 4});
  EXPECT_EQ(mesh01.getLocalYIndex(5), 0);
  EXPECT_EQ(mesh01.getLocalYIndex(6), 1);
  EXPECT_EQ(mesh01.getLocalYIndex(7), 2);
  EXPECT_EQ(mesh01.getLocalYIndex(8), 3);
  EXPECT_EQ(mesh01.getLocalYIndex(9), 4);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  mesh10.setYDecompositionIndices({-1, 2, 5, 5, 4});
  EXPECT_EQ(mesh10.getLocalYIndex(0), 0);
  EXPECT_EQ(mesh10.getLocalYIndex(1), 1);
  EXPECT_EQ(mesh10.getLocalYIndex(2), 2);
  EXPECT_EQ(mesh10.getLocalYIndex(3), 3);
  EXPECT_EQ(mesh10.getLocalYIndex(4), 4);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  mesh11.setYDecompositionIndices({-1, 2, 5, 5, 4});
  EXPECT_EQ(mesh11.getLocalYIndex(5), 0);
  EXPECT_EQ(mesh11.getLocalYIndex(6), 1);
  EXPECT_EQ(mesh11.getLocalYIndex(7), 2);
  EXPECT_EQ(mesh11.getLocalYIndex(8), 3);
  EXPECT_EQ(mesh11.getLocalYIndex(9), 4);
}

TEST(BoutMeshTest, GetLocalYIndexNoBoundaries) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Local indices start counting from the first non-boundary point

  // |<--  1st Y-proc -->|
  //             |<--  2nd Y-proc -->|
  // +---+---+---+---+---+---+---+---+
  // |-1*| 0 | 1 | 2 | 3 | 4 | 5 | 6*| <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st Y-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd Y-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getLocalYIndexNoBoundaries(-1), 0);
  EXPECT_EQ(mesh00.getLocalYIndexNoBoundaries(0), 1);
  EXPECT_EQ(mesh00.getLocalYIndexNoBoundaries(1), 2);
  EXPECT_EQ(mesh00.getLocalYIndexNoBoundaries(2), 3);
  EXPECT_EQ(mesh00.getLocalYIndexNoBoundaries(3), 4);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getLocalYIndexNoBoundaries(2), 0);
  EXPECT_EQ(mesh01.getLocalYIndexNoBoundaries(3), 1);
  EXPECT_EQ(mesh01.getLocalYIndexNoBoundaries(4), 2);
  EXPECT_EQ(mesh01.getLocalYIndexNoBoundaries(5), 3);
  EXPECT_EQ(mesh01.getLocalYIndexNoBoundaries(6), 4);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getLocalYIndexNoBoundaries(-1), 0);
  EXPECT_EQ(mesh10.getLocalYIndexNoBoundaries(0), 1);
  EXPECT_EQ(mesh10.getLocalYIndexNoBoundaries(1), 2);
  EXPECT_EQ(mesh10.getLocalYIndexNoBoundaries(2), 3);
  EXPECT_EQ(mesh10.getLocalYIndexNoBoundaries(3), 4);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getLocalYIndexNoBoundaries(2), 0);
  EXPECT_EQ(mesh11.getLocalYIndexNoBoundaries(3), 1);
  EXPECT_EQ(mesh11.getLocalYIndexNoBoundaries(4), 2);
  EXPECT_EQ(mesh11.getLocalYIndexNoBoundaries(5), 3);
  EXPECT_EQ(mesh11.getLocalYIndexNoBoundaries(6), 4);
}

TEST(BoutMeshTest, GetGlobalZIndex) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Boundaries are included in the global index

  // No parallelisation in Z, so function is just the identity

  BoutMeshExposer mesh00(5, 3, 4, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getGlobalZIndex(0), 0);
  EXPECT_EQ(mesh00.getGlobalZIndex(1), 1);
  EXPECT_EQ(mesh00.getGlobalZIndex(2), 2);
  EXPECT_EQ(mesh00.getGlobalZIndex(3), 3);
  EXPECT_EQ(mesh00.getGlobalZIndex(4), 4);

  BoutMeshExposer mesh01(5, 3, 4, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getGlobalZIndex(0), 0);
  EXPECT_EQ(mesh01.getGlobalZIndex(1), 1);
  EXPECT_EQ(mesh01.getGlobalZIndex(2), 2);
  EXPECT_EQ(mesh01.getGlobalZIndex(3), 3);
  EXPECT_EQ(mesh01.getGlobalZIndex(4), 4);

  BoutMeshExposer mesh10(5, 3, 4, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getGlobalZIndex(0), 0);
  EXPECT_EQ(mesh10.getGlobalZIndex(1), 1);
  EXPECT_EQ(mesh10.getGlobalZIndex(2), 2);
  EXPECT_EQ(mesh10.getGlobalZIndex(3), 3);
  EXPECT_EQ(mesh10.getGlobalZIndex(4), 4);

  BoutMeshExposer mesh11(5, 3, 4, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getGlobalZIndex(0), 0);
  EXPECT_EQ(mesh11.getGlobalZIndex(1), 1);
  EXPECT_EQ(mesh11.getGlobalZIndex(2), 2);
  EXPECT_EQ(mesh11.getGlobalZIndex(3), 3);
  EXPECT_EQ(mesh11.getGlobalZIndex(4), 4);
}

TEST(BoutMeshTest, GetGlobalZIndexNoBoundaries) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  BoutMeshExposer mesh00(5, 3, 4, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getGlobalZIndexNoBoundaries(0), 0);
  EXPECT_EQ(mesh00.getGlobalZIndexNoBoundaries(1), 1);
  EXPECT_EQ(mesh00.getGlobalZIndexNoBoundaries(2), 2);
  EXPECT_EQ(mesh00.getGlobalZIndexNoBoundaries(3), 3);
  EXPECT_EQ(mesh00.getGlobalZIndexNoBoundaries(4), 4);

  BoutMeshExposer mesh01(5, 3, 4, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getGlobalZIndexNoBoundaries(0), 0);
  EXPECT_EQ(mesh01.getGlobalZIndexNoBoundaries(1), 1);
  EXPECT_EQ(mesh01.getGlobalZIndexNoBoundaries(2), 2);
  EXPECT_EQ(mesh01.getGlobalZIndexNoBoundaries(3), 3);
  EXPECT_EQ(mesh01.getGlobalZIndexNoBoundaries(4), 4);

  BoutMeshExposer mesh10(5, 3, 4, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getGlobalZIndexNoBoundaries(0), 0);
  EXPECT_EQ(mesh10.getGlobalZIndexNoBoundaries(1), 1);
  EXPECT_EQ(mesh10.getGlobalZIndexNoBoundaries(2), 2);
  EXPECT_EQ(mesh10.getGlobalZIndexNoBoundaries(3), 3);
  EXPECT_EQ(mesh10.getGlobalZIndexNoBoundaries(4), 4);

  BoutMeshExposer mesh11(5, 3, 4, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getGlobalZIndexNoBoundaries(0), 0);
  EXPECT_EQ(mesh11.getGlobalZIndexNoBoundaries(1), 1);
  EXPECT_EQ(mesh11.getGlobalZIndexNoBoundaries(2), 2);
  EXPECT_EQ(mesh11.getGlobalZIndexNoBoundaries(3), 3);
  EXPECT_EQ(mesh11.getGlobalZIndexNoBoundaries(4), 4);
}

TEST(BoutMeshTest, GetLocalZIndex) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Boundaries are included in the local index

  // |<--  1st Z-proc -->|
  //             |<--  2nd Z-proc -->|
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st Z-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd Z-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 4, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getLocalZIndex(0), 0);
  EXPECT_EQ(mesh00.getLocalZIndex(1), 1);
  EXPECT_EQ(mesh00.getLocalZIndex(2), 2);
  EXPECT_EQ(mesh00.getLocalZIndex(3), 3);
  EXPECT_EQ(mesh00.getLocalZIndex(4), 4);

  BoutMeshExposer mesh01(5, 3, 4, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getLocalZIndex(0), 0);
  EXPECT_EQ(mesh01.getLocalZIndex(1), 1);
  EXPECT_EQ(mesh01.getLocalZIndex(2), 2);
  EXPECT_EQ(mesh01.getLocalZIndex(3), 3);
  EXPECT_EQ(mesh01.getLocalZIndex(4), 4);

  BoutMeshExposer mesh10(5, 3, 4, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getLocalZIndex(0), 0);
  EXPECT_EQ(mesh10.getLocalZIndex(1), 1);
  EXPECT_EQ(mesh10.getLocalZIndex(2), 2);
  EXPECT_EQ(mesh10.getLocalZIndex(3), 3);
  EXPECT_EQ(mesh10.getLocalZIndex(4), 4);

  BoutMeshExposer mesh11(5, 3, 4, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getLocalZIndex(0), 0);
  EXPECT_EQ(mesh11.getLocalZIndex(1), 1);
  EXPECT_EQ(mesh11.getLocalZIndex(2), 2);
  EXPECT_EQ(mesh11.getLocalZIndex(3), 3);
  EXPECT_EQ(mesh11.getLocalZIndex(4), 4);
}

TEST(BoutMeshTest, GetLocalZIndexNoBoundaries) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  // 2x2 processors, 3x3x1 (not including guards) on each processor

  // Local indices start counting from the first non-boundary point

  // |<--  1st Z-proc -->|
  //             |<--  2nd Z-proc -->|
  // +---+---+---+---+---+---+---+---+
  // |-1*| 0 | 1 | 2 | 3 | 4 | 5 | 6*| <- Global indices
  // +---+---+---+---+---+---+---+---+
  // | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | <- 1st Z-processor
  // +---+---+---+---+---+---+---+---+
  // |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | <- 2nd Z-processor
  // +---+---+---+---+---+---+---+---+

  BoutMeshExposer mesh00(5, 3, 4, 2, 2, 0, 0);
  EXPECT_EQ(mesh00.getLocalZIndexNoBoundaries(0), 0);
  EXPECT_EQ(mesh00.getLocalZIndexNoBoundaries(1), 1);
  EXPECT_EQ(mesh00.getLocalZIndexNoBoundaries(2), 2);
  EXPECT_EQ(mesh00.getLocalZIndexNoBoundaries(3), 3);
  EXPECT_EQ(mesh00.getLocalZIndexNoBoundaries(4), 4);

  BoutMeshExposer mesh01(5, 3, 4, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.getLocalZIndexNoBoundaries(0), 0);
  EXPECT_EQ(mesh01.getLocalZIndexNoBoundaries(1), 1);
  EXPECT_EQ(mesh01.getLocalZIndexNoBoundaries(2), 2);
  EXPECT_EQ(mesh01.getLocalZIndexNoBoundaries(3), 3);
  EXPECT_EQ(mesh01.getLocalZIndexNoBoundaries(4), 4);

  BoutMeshExposer mesh10(5, 3, 4, 2, 2, 1, 0);
  EXPECT_EQ(mesh10.getLocalZIndexNoBoundaries(0), 0);
  EXPECT_EQ(mesh10.getLocalZIndexNoBoundaries(1), 1);
  EXPECT_EQ(mesh10.getLocalZIndexNoBoundaries(2), 2);
  EXPECT_EQ(mesh10.getLocalZIndexNoBoundaries(3), 3);
  EXPECT_EQ(mesh10.getLocalZIndexNoBoundaries(4), 4);

  BoutMeshExposer mesh11(5, 3, 4, 2, 2, 1, 1);
  EXPECT_EQ(mesh11.getLocalZIndexNoBoundaries(0), 0);
  EXPECT_EQ(mesh11.getLocalZIndexNoBoundaries(1), 1);
  EXPECT_EQ(mesh11.getLocalZIndexNoBoundaries(2), 2);
  EXPECT_EQ(mesh11.getLocalZIndexNoBoundaries(3), 3);
  EXPECT_EQ(mesh11.getLocalZIndexNoBoundaries(4), 4);
}

TEST(BoutMeshTest, FirstX) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh00(5, 3, 4, 3, 3, 0, 0);
  EXPECT_TRUE(mesh00.firstX());
  BoutMeshExposer mesh10(5, 3, 4, 3, 3, 1, 0);
  EXPECT_FALSE(mesh10.firstX());
  BoutMeshExposer mesh20(5, 3, 4, 3, 3, 2, 0);
  EXPECT_FALSE(mesh20.firstX());
  BoutMeshExposer mesh01(5, 3, 4, 3, 3, 0, 1);
  EXPECT_TRUE(mesh01.firstX());
  BoutMeshExposer mesh11(5, 3, 4, 3, 3, 1, 1);
  EXPECT_FALSE(mesh11.firstX());
  BoutMeshExposer mesh21(5, 3, 4, 3, 3, 2, 1);
  EXPECT_FALSE(mesh21.firstX());
  BoutMeshExposer mesh02(5, 3, 4, 3, 3, 0, 2);
  EXPECT_TRUE(mesh02.firstX());
  BoutMeshExposer mesh12(5, 3, 4, 3, 3, 1, 2);
  EXPECT_FALSE(mesh12.firstX());
  BoutMeshExposer mesh22(5, 3, 4, 3, 3, 2, 2);
  EXPECT_FALSE(mesh22.firstX());
}

TEST(BoutMeshTest, LastX) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh00(5, 3, 4, 3, 3, 0, 0);
  EXPECT_FALSE(mesh00.lastX());
  BoutMeshExposer mesh10(5, 3, 4, 3, 3, 1, 0);
  EXPECT_FALSE(mesh10.lastX());
  BoutMeshExposer mesh20(5, 3, 4, 3, 3, 2, 0);
  EXPECT_TRUE(mesh20.lastX());
  BoutMeshExposer mesh01(5, 3, 4, 3, 3, 0, 1);
  EXPECT_FALSE(mesh01.lastX());
  BoutMeshExposer mesh11(5, 3, 4, 3, 3, 1, 1);
  EXPECT_FALSE(mesh11.lastX());
  BoutMeshExposer mesh21(5, 3, 4, 3, 3, 2, 1);
  EXPECT_TRUE(mesh21.lastX());
  BoutMeshExposer mesh02(5, 3, 4, 3, 3, 0, 2);
  EXPECT_FALSE(mesh02.lastX());
  BoutMeshExposer mesh12(5, 3, 4, 3, 3, 1, 2);
  EXPECT_FALSE(mesh12.lastX());
  BoutMeshExposer mesh22(5, 3, 4, 3, 3, 2, 2);
  EXPECT_TRUE(mesh22.lastX());
}

TEST(BoutMeshTest, FirstY) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh00(5, 3, 4, 3, 3, 0, 0);
  EXPECT_TRUE(mesh00.firstY());
  BoutMeshExposer mesh10(5, 3, 4, 3, 3, 1, 0);
  EXPECT_TRUE(mesh10.firstY());
  BoutMeshExposer mesh20(5, 3, 4, 3, 3, 2, 0);
  EXPECT_TRUE(mesh20.firstY());
  BoutMeshExposer mesh01(5, 3, 4, 3, 3, 0, 1);
  EXPECT_FALSE(mesh01.firstY());
  BoutMeshExposer mesh11(5, 3, 4, 3, 3, 1, 1);
  EXPECT_FALSE(mesh11.firstY());
  BoutMeshExposer mesh21(5, 3, 4, 3, 3, 2, 1);
  EXPECT_FALSE(mesh21.firstY());
  BoutMeshExposer mesh02(5, 3, 4, 3, 3, 0, 2);
  EXPECT_FALSE(mesh02.firstY());
  BoutMeshExposer mesh12(5, 3, 4, 3, 3, 1, 2);
  EXPECT_FALSE(mesh12.firstY());
  BoutMeshExposer mesh22(5, 3, 4, 3, 3, 2, 2);
  EXPECT_FALSE(mesh22.firstY());
}

TEST(BoutMeshTest, LastY) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh00(5, 3, 4, 3, 3, 0, 0);
  EXPECT_FALSE(mesh00.lastY());
  BoutMeshExposer mesh10(5, 3, 4, 3, 3, 1, 0);
  EXPECT_FALSE(mesh10.lastY());
  BoutMeshExposer mesh20(5, 3, 4, 3, 3, 2, 0);
  EXPECT_FALSE(mesh20.lastY());
  BoutMeshExposer mesh01(5, 3, 4, 3, 3, 0, 1);
  EXPECT_FALSE(mesh01.lastY());
  BoutMeshExposer mesh11(5, 3, 4, 3, 3, 1, 1);
  EXPECT_FALSE(mesh11.lastY());
  BoutMeshExposer mesh21(5, 3, 4, 3, 3, 2, 1);
  EXPECT_FALSE(mesh21.lastY());
  BoutMeshExposer mesh02(5, 3, 4, 3, 3, 0, 2);
  EXPECT_TRUE(mesh02.lastY());
  BoutMeshExposer mesh12(5, 3, 4, 3, 3, 1, 2);
  EXPECT_TRUE(mesh12.lastY());
  BoutMeshExposer mesh22(5, 3, 4, 3, 3, 2, 2);
  EXPECT_TRUE(mesh22.lastY());
}

TEST(BoutMeshTest, DefaultConnectionsSingleNull1x1) {
  WithQuietOutput info{output_info};
  // 5x3x1 grid on 1 processor, 1 boundary point. Boundaries should be
  // simple 1D rectangles, with 4 boundaries on this processor
  BoutMeshExposer mesh00(5, 3, 1, 1, 1, 0, 0, false);

  mesh00.default_connections();

  auto connection00 = mesh00.getConnectionInfo();
  EXPECT_EQ(connection00.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection00.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection00.DDATA_INDEST, -1);
  EXPECT_EQ(connection00.UDATA_INDEST, -1);
  EXPECT_EQ(connection00.DDATA_OUTDEST, -1);
  EXPECT_EQ(connection00.UDATA_OUTDEST, -1);
  EXPECT_EQ(connection00.IDATA_DEST, -1);
  EXPECT_EQ(connection00.ODATA_DEST, -1);

  mesh00.createDefaultRegions();
  mesh00.addBoundaryRegions();

  EXPECT_EQ( mesh00.getRegion("RGN_LOWER_INNER_Y").size(), 5);
  EXPECT_EQ( mesh00.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ( mesh00.getRegion("RGN_LOWER_Y").size(), 5);

  EXPECT_EQ( mesh00.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ( mesh00.getRegion("RGN_UPPER_OUTER_Y").size(), 5);
  EXPECT_EQ( mesh00.getRegion("RGN_UPPER_Y").size(), 5);

  EXPECT_EQ( mesh00.getRegion("RGN_INNER_X").size(), 3);
  EXPECT_EQ( mesh00.getRegion("RGN_OUTER_X").size(), 3);
}

TEST(BoutMeshTest, DefaultConnectionsSingleNull2x2) {
  // This test checks both default_connections and the Region
  // creation, as these are quite tightly linked.

  WithQuietOutput info{output_info};
  // 5x3x1 local grid on 4 processors, 1 boundary point. Boundaries should
  // be simple 1D rectangles, with 2 boundaries on each processor.

  // +---+---+  ^
  // | 0 | 1 |  | DDATA_OUTDEST
  // +---+---+
  // | 2 | 3 |  | UDATA_OUTDEST
  // +---+---+  v
  // <--
  // IDATA_DEST
  //       -->
  // ODATA_DEST
  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0, false);
  mesh00.default_connections();

  auto connection00 = mesh00.getConnectionInfo();
  EXPECT_EQ(connection00.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection00.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection00.DDATA_INDEST, -1);
  EXPECT_EQ(connection00.UDATA_INDEST, -1);
  EXPECT_EQ(connection00.DDATA_OUTDEST, -1);
  EXPECT_EQ(connection00.UDATA_OUTDEST, 2);
  EXPECT_EQ(connection00.IDATA_DEST, -1);
  EXPECT_EQ(connection00.ODATA_DEST, 1);

  mesh00.createDefaultRegions();
  mesh00.addBoundaryRegions();

  EXPECT_EQ(mesh00.getRegion("RGN_LOWER_INNER_Y").size(), 4);
  EXPECT_EQ(mesh00.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_LOWER_Y").size(), 4);

  EXPECT_EQ(mesh00.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_UPPER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_UPPER_Y").size(), 0);

  EXPECT_EQ(mesh00.getRegion("RGN_INNER_X").size(), 3);
  EXPECT_EQ(mesh00.getRegion("RGN_OUTER_X").size(), 0);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1, false);
  mesh01.default_connections();

  auto connection01 = mesh01.getConnectionInfo();
  EXPECT_EQ(connection01.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection01.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection01.DDATA_INDEST, -1);
  EXPECT_EQ(connection01.UDATA_INDEST, -1);
  EXPECT_EQ(connection01.DDATA_OUTDEST, 0);
  EXPECT_EQ(connection01.UDATA_OUTDEST, -1);
  EXPECT_EQ(connection01.IDATA_DEST, -1);
  EXPECT_EQ(connection01.ODATA_DEST, 3);

  mesh01.createDefaultRegions();
  mesh01.addBoundaryRegions();

  EXPECT_EQ(mesh01.getRegion("RGN_LOWER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_LOWER_Y").size(), 0);

  EXPECT_EQ(mesh01.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_UPPER_OUTER_Y").size(), 4);
  EXPECT_EQ(mesh01.getRegion("RGN_UPPER_Y").size(), 4);

  EXPECT_EQ(mesh01.getRegion("RGN_INNER_X").size(), 3);
  EXPECT_EQ(mesh01.getRegion("RGN_OUTER_X").size(), 0);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0, false);
  mesh10.default_connections();

  auto connection10 = mesh10.getConnectionInfo();
  EXPECT_EQ(connection10.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection10.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection10.DDATA_INDEST, -1);
  EXPECT_EQ(connection10.UDATA_INDEST, -1);
  EXPECT_EQ(connection10.DDATA_OUTDEST, -1);
  EXPECT_EQ(connection10.UDATA_OUTDEST, 3);
  EXPECT_EQ(connection10.IDATA_DEST, 0);
  EXPECT_EQ(connection10.ODATA_DEST, -1);

  mesh10.createDefaultRegions();
  mesh10.addBoundaryRegions();

  EXPECT_EQ(mesh10.getRegion("RGN_LOWER_INNER_Y").size(), 4);
  EXPECT_EQ(mesh10.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_LOWER_Y").size(), 4);

  EXPECT_EQ(mesh10.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_UPPER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_UPPER_Y").size(), 0);

  EXPECT_EQ(mesh10.getRegion("RGN_INNER_X").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_OUTER_X").size(), 3);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1, false);
  mesh11.default_connections();

  auto connection11 = mesh11.getConnectionInfo();
  EXPECT_EQ(connection11.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection11.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection11.DDATA_INDEST, -1);
  EXPECT_EQ(connection11.UDATA_INDEST, -1);
  EXPECT_EQ(connection11.DDATA_OUTDEST, 1);
  EXPECT_EQ(connection11.UDATA_OUTDEST, -1);
  EXPECT_EQ(connection11.IDATA_DEST, 2);
  EXPECT_EQ(connection11.ODATA_DEST, -1);

  mesh11.createDefaultRegions();
  mesh11.addBoundaryRegions();

  EXPECT_EQ(mesh11.getRegion("RGN_LOWER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_LOWER_Y").size(), 0);

  EXPECT_EQ(mesh11.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_UPPER_OUTER_Y").size(), 4);
  EXPECT_EQ(mesh11.getRegion("RGN_UPPER_Y").size(), 4);

  EXPECT_EQ(mesh11.getRegion("RGN_INNER_X").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_OUTER_X").size(), 3);
}

TEST(BoutMeshTest, DefaultConnectionsSingleNullPeriodicX2x2) {
  // This test checks both default_connections and the Region
  // creation, as these are quite tightly linked.

  WithQuietOutput info{output_info};
  // 5x3x1 local grid on 4 processors, 1 boundary point. Boundaries should
  // be simple 1D rectangles, with 2 boundaries on each processor.
  // Periodic in X, so {I,O}DATA_DEST wrap around, and no boundaries in X

  // +---+---+  ^
  // | 0 | 1 |  | DDATA_OUTDEST
  // +---+---+
  // | 2 | 3 |  | UDATA_OUTDEST
  // +---+---+  v
  // <--
  // IDATA_DEST
  //       -->
  // ODATA_DEST
  BoutMeshExposer mesh00(5, 3, 1, 2, 2, 0, 0, false);
  mesh00.periodicX = true;
  mesh00.default_connections();

  auto connection00 = mesh00.getConnectionInfo();
  EXPECT_EQ(connection00.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection00.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection00.DDATA_INDEST, -1);
  EXPECT_EQ(connection00.UDATA_INDEST, -1);
  EXPECT_EQ(connection00.DDATA_OUTDEST, -1);
  EXPECT_EQ(connection00.UDATA_OUTDEST, 2);
  EXPECT_EQ(connection00.IDATA_DEST, 1);
  EXPECT_EQ(connection00.ODATA_DEST, 1);

  mesh00.createDefaultRegions();
  mesh00.addBoundaryRegions();

  EXPECT_EQ(mesh00.getRegion("RGN_LOWER_INNER_Y").size(), 4);
  EXPECT_EQ(mesh00.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_LOWER_Y").size(), 4);

  EXPECT_EQ(mesh00.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_UPPER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_UPPER_Y").size(), 0);

  EXPECT_EQ(mesh00.getRegion("RGN_INNER_X").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_OUTER_X").size(), 0);

  BoutMeshExposer mesh01(5, 3, 1, 2, 2, 0, 1, false);
  mesh01.periodicX = true;
  mesh01.default_connections();

  auto connection01 = mesh01.getConnectionInfo();
  EXPECT_EQ(connection01.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection01.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection01.DDATA_INDEST, -1);
  EXPECT_EQ(connection01.UDATA_INDEST, -1);
  EXPECT_EQ(connection01.DDATA_OUTDEST, 0);
  EXPECT_EQ(connection01.UDATA_OUTDEST, -1);
  EXPECT_EQ(connection01.IDATA_DEST, 3);
  EXPECT_EQ(connection01.ODATA_DEST, 3);

  mesh01.createDefaultRegions();
  mesh01.addBoundaryRegions();

  EXPECT_EQ(mesh01.getRegion("RGN_LOWER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_LOWER_Y").size(), 0);

  EXPECT_EQ(mesh01.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_UPPER_OUTER_Y").size(), 4);
  EXPECT_EQ(mesh01.getRegion("RGN_UPPER_Y").size(), 4);

  EXPECT_EQ(mesh01.getRegion("RGN_INNER_X").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_OUTER_X").size(), 0);

  BoutMeshExposer mesh10(5, 3, 1, 2, 2, 1, 0, false);
  mesh10.periodicX = true;
  mesh10.default_connections();

  auto connection10 = mesh10.getConnectionInfo();
  EXPECT_EQ(connection10.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection10.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection10.DDATA_INDEST, -1);
  EXPECT_EQ(connection10.UDATA_INDEST, -1);
  EXPECT_EQ(connection10.DDATA_OUTDEST, -1);
  EXPECT_EQ(connection10.UDATA_OUTDEST, 3);
  EXPECT_EQ(connection10.IDATA_DEST, 0);
  EXPECT_EQ(connection10.ODATA_DEST, 0);

  mesh10.createDefaultRegions();
  mesh10.addBoundaryRegions();

  EXPECT_EQ(mesh10.getRegion("RGN_LOWER_INNER_Y").size(), 4);
  EXPECT_EQ(mesh10.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_LOWER_Y").size(), 4);

  EXPECT_EQ(mesh10.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_UPPER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_UPPER_Y").size(), 0);

  EXPECT_EQ(mesh10.getRegion("RGN_INNER_X").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_OUTER_X").size(), 0);

  BoutMeshExposer mesh11(5, 3, 1, 2, 2, 1, 1, false);
  mesh11.periodicX = true;
  mesh11.default_connections();

  auto connection11 = mesh11.getConnectionInfo();
  EXPECT_EQ(connection11.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection11.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection11.DDATA_INDEST, -1);
  EXPECT_EQ(connection11.UDATA_INDEST, -1);
  EXPECT_EQ(connection11.DDATA_OUTDEST, 1);
  EXPECT_EQ(connection11.UDATA_OUTDEST, -1);
  EXPECT_EQ(connection11.IDATA_DEST, 2);
  EXPECT_EQ(connection11.ODATA_DEST, 2);

  mesh11.createDefaultRegions();
  mesh11.addBoundaryRegions();

  EXPECT_EQ(mesh11.getRegion("RGN_LOWER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_LOWER_Y").size(), 0);

  EXPECT_EQ(mesh11.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_UPPER_OUTER_Y").size(), 4);
  EXPECT_EQ(mesh11.getRegion("RGN_UPPER_Y").size(), 4);

  EXPECT_EQ(mesh11.getRegion("RGN_INNER_X").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_OUTER_X").size(), 0);
}

TEST(BoutMeshTest, TopologySingleNull2x2) {
  WithQuietOutput info{output_info};
  // 32x24x1 local grid on 4 processors, 1 boundary point. Boundaries
  // are now a bit more complicated

  // +---+---+  ^
  // | 0 | 1 |  | DDATA_OUTDEST
  // +---+---+
  // | 2 | 3 |  | UDATA_OUTDEST
  // +---+---+  v
  // <--
  // IDATA_DEST
  //       -->
  // ODATA_DEST
  BoutMeshExposer mesh00(32, 24, 1, 2, 2, 0, 0, false);
  mesh00.setYDecompositionIndices(-1, 12, 12, 23, 12);
  mesh00.topology();

  auto connection00 = mesh00.getConnectionInfo();
  EXPECT_EQ(connection00.DDATA_XSPLIT, 32);
  EXPECT_EQ(connection00.UDATA_XSPLIT, 32);
  EXPECT_EQ(connection00.DDATA_INDEST, 0);
  EXPECT_EQ(connection00.UDATA_INDEST, 0);
  EXPECT_EQ(connection00.DDATA_OUTDEST, -1);
  EXPECT_EQ(connection00.UDATA_OUTDEST, -1);
  EXPECT_EQ(connection00.IDATA_DEST, -1);
  EXPECT_EQ(connection00.ODATA_DEST, 1);

  mesh00.createDefaultRegions();
  mesh00.addBoundaryRegions();

  EXPECT_EQ(mesh00.getRegion("RGN_LOWER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_LOWER_Y").size(), 0);

  EXPECT_EQ(mesh00.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_UPPER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh00.getRegion("RGN_UPPER_Y").size(), 0);

  EXPECT_EQ(mesh00.getRegion("RGN_INNER_X").size(), 24);
  EXPECT_EQ(mesh00.getRegion("RGN_OUTER_X").size(), 0);

  BoutMeshExposer mesh01(32, 24, 1, 2, 2, 0, 1, false);
  mesh01.setYDecompositionIndices(-1, 12, 12, 23, 12);
  mesh01.topology();

  auto connection01 = mesh01.getConnectionInfo();
  EXPECT_EQ(connection01.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection01.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection01.DDATA_INDEST, -1);
  EXPECT_EQ(connection01.UDATA_INDEST, -1);
  EXPECT_EQ(connection01.DDATA_OUTDEST, 0);
  EXPECT_EQ(connection01.UDATA_OUTDEST, -1);
  EXPECT_EQ(connection01.IDATA_DEST, -1);
  EXPECT_EQ(connection01.ODATA_DEST, 3);

  mesh01.createDefaultRegions();
  mesh01.addBoundaryRegions();

  EXPECT_EQ(mesh01.getRegion("RGN_LOWER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_LOWER_Y").size(), 0);

  EXPECT_EQ(mesh01.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh01.getRegion("RGN_UPPER_OUTER_Y").size(), 31);
  EXPECT_EQ(mesh01.getRegion("RGN_UPPER_Y").size(), 31);

  EXPECT_EQ(mesh01.getRegion("RGN_INNER_X").size(), 24);
  EXPECT_EQ(mesh01.getRegion("RGN_OUTER_X").size(), 0);

  BoutMeshExposer mesh10(32, 24, 1, 2, 2, 1, 0, false);
  mesh10.setYDecompositionIndices(-1, 12, 12, 23, 12);
  mesh10.topology();

  auto connection10 = mesh10.getConnectionInfo();
  EXPECT_EQ(connection10.DDATA_XSPLIT, 32);
  EXPECT_EQ(connection10.UDATA_XSPLIT, 32);
  EXPECT_EQ(connection10.DDATA_INDEST, 1);
  EXPECT_EQ(connection10.UDATA_INDEST, 1);
  EXPECT_EQ(connection10.DDATA_OUTDEST, -1);
  EXPECT_EQ(connection10.UDATA_OUTDEST, -1);
  EXPECT_EQ(connection10.IDATA_DEST, 0);
  EXPECT_EQ(connection10.ODATA_DEST, -1);

  mesh10.createDefaultRegions();
  mesh10.addBoundaryRegions();

  EXPECT_EQ(mesh10.getRegion("RGN_LOWER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_LOWER_Y").size(), 0);

  EXPECT_EQ(mesh10.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_UPPER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_UPPER_Y").size(), 0);

  EXPECT_EQ(mesh10.getRegion("RGN_INNER_X").size(), 0);
  EXPECT_EQ(mesh10.getRegion("RGN_OUTER_X").size(), 24);

  BoutMeshExposer mesh11(32, 24, 1, 2, 2, 1, 1, false);
  mesh11.setYDecompositionIndices(-1, 12, 12, 23, 12);
  mesh11.topology();

  auto connection11 = mesh11.getConnectionInfo();
  EXPECT_EQ(connection11.DDATA_XSPLIT, 0);
  EXPECT_EQ(connection11.UDATA_XSPLIT, 0);
  EXPECT_EQ(connection11.DDATA_INDEST, -1);
  EXPECT_EQ(connection11.UDATA_INDEST, -1);
  EXPECT_EQ(connection11.DDATA_OUTDEST, 1);
  EXPECT_EQ(connection11.UDATA_OUTDEST, -1);
  EXPECT_EQ(connection11.IDATA_DEST, 2);
  EXPECT_EQ(connection11.ODATA_DEST, -1);

  mesh11.createDefaultRegions();
  mesh11.addBoundaryRegions();

  EXPECT_EQ(mesh11.getRegion("RGN_LOWER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_LOWER_OUTER_Y").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_LOWER_Y").size(), 0);

  EXPECT_EQ(mesh11.getRegion("RGN_UPPER_INNER_Y").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_UPPER_OUTER_Y").size(), 31);
  EXPECT_EQ(mesh11.getRegion("RGN_UPPER_Y").size(), 31);

  EXPECT_EQ(mesh11.getRegion("RGN_INNER_X").size(), 0);
  EXPECT_EQ(mesh11.getRegion("RGN_OUTER_X").size(), 24);
}
