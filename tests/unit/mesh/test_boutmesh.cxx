#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "options.hxx"
#include "output.hxx"
#include "bout/griddata.hxx"

#include "test_extras.hxx"

#include <array>
#include <ostream>

/// Forward declaration so we can construct a `BoutMeshExposer` from this
struct BoutMeshParameters;

/// Inherits from `BoutMesh` so that we can make some protected things
/// public to aid testing
class BoutMeshExposer : public BoutMesh {
public:
  BoutMeshExposer(int input_nx, int input_ny, int input_nz, int mxg, int myg,
                  int input_npes = 1)
      : BoutMesh(input_nx, input_ny, input_nz, mxg, myg, input_npes) {}
  BoutMeshExposer(int nx, int ny, int nz, int nxpe, int nype, int pe_xind, int pe_yind,
                  bool create_topology = true, bool symmetric_X = true,
                  bool symmetric_Y = true)
      : BoutMesh((nxpe * (nx - 2)) + 2, nype * ny, nz, 1, 1, nxpe, nype, pe_xind, pe_yind,
                 create_topology, symmetric_X, symmetric_Y) {}
  BoutMeshExposer(const BoutMeshParameters& inputs, bool periodicX_ = false);
  // Make protected methods public for testing
  using BoutMesh::add_target;
  using BoutMesh::addBoundaryRegions;
  using BoutMesh::chooseProcessorSplit;
  using BoutMesh::ConnectionInfo;
  using BoutMesh::createXBoundaries;
  using BoutMesh::createYBoundaries;
  using BoutMesh::default_connections;
  using BoutMesh::findProcessorSplit;
  using BoutMesh::getConnectionInfo;
  using BoutMesh::PROC_NUM;
  using BoutMesh::set_connection;
  using BoutMesh::setShiftAngle;
  using BoutMesh::setXDecompositionIndices;
  using BoutMesh::setYDecompositionIndices;
  using BoutMesh::topology;
  using BoutMesh::XDecompositionIndices;
  using BoutMesh::XPROC;
  using BoutMesh::YDecompositionIndices;
  using BoutMesh::YPROC;
};

/// Minimal parameters need to construct a grid useful for testing
struct BoutMeshGridInfo {
  int local_nx; // Does _not_ include guard cells
  int local_ny; // Does _not_ include guard cells
  int num_x_guards;
  int num_y_guards;
  int nxpe;
  int nype;
  int pe_xind;
  int pe_yind;
  bool symmetric_X;
  bool symmetric_Y;
  // The below are constructed consistently with the above
  int total_nx; // _Does_ include guard cells
  int total_ny; // Does _not_ include guard cells
  int total_processors;
  BoutMeshGridInfo(int local_nx_, int local_ny_, int num_x_guards_, int num_y_guards_,
                   int nxpe_, int nype_, int pe_xind_ = 0, int pe_yind_ = 0,
                   bool symmetric_X_ = true, bool symmetric_Y_ = true)
      : local_nx(local_nx_), local_ny(local_ny_), num_x_guards(num_x_guards_),
        num_y_guards(num_y_guards_), nxpe(nxpe_), nype(nype_), pe_xind(pe_xind_),
        pe_yind(pe_yind_), symmetric_X(symmetric_X_), symmetric_Y(symmetric_Y_),
        total_nx((nxpe * local_nx) + (2 * num_x_guards)), total_ny(nype * local_ny),
        total_processors(nxpe * nype) {}
};

/// Grid and topology information to make a `BoutMesh`
struct BoutMeshParameters {
  BoutMeshGridInfo grid;
  BoutMeshExposer::XDecompositionIndices x_indices;
  BoutMeshExposer::YDecompositionIndices y_indices;
};

/// Now we've got the definition of `BoutMeshParameters`, we can
/// actually make a `BoutMeshExposer`
BoutMeshExposer::BoutMeshExposer(const BoutMeshParameters& inputs, bool periodicX_)
    : BoutMesh(inputs.grid.total_nx, inputs.grid.total_ny, 1, inputs.grid.num_x_guards,
               inputs.grid.num_y_guards, inputs.grid.nxpe, inputs.grid.nype,
               inputs.grid.pe_xind, inputs.grid.pe_yind, inputs.grid.symmetric_X,
               inputs.grid.symmetric_Y, periodicX_, inputs.x_indices.ixseps1,
               inputs.x_indices.ixseps2, inputs.y_indices.jyseps1_1,
               inputs.y_indices.jyseps2_1, inputs.y_indices.jyseps1_2,
               inputs.y_indices.jyseps2_2, inputs.y_indices.ny_inner) {}

/// Equality operator to help testing
bool operator==(const BoutMeshExposer::YDecompositionIndices& lhs,
                const BoutMeshExposer::YDecompositionIndices& rhs) {
  return (lhs.jyseps1_1 == rhs.jyseps1_1) and (lhs.jyseps2_1 == rhs.jyseps2_1)
         and (lhs.jyseps1_2 == rhs.jyseps1_2) and (lhs.jyseps2_2 == rhs.jyseps2_2)
         and (lhs.ny_inner == rhs.ny_inner);
}

bool operator==(const BoutMeshExposer::ConnectionInfo& lhs,
                const BoutMeshExposer::ConnectionInfo& rhs) {
  return (lhs.TS_up_in == rhs.TS_up_in) and (lhs.TS_up_out == rhs.TS_up_out)
         and (lhs.TS_down_in == rhs.TS_down_in) and (lhs.TS_down_out == rhs.TS_down_out)
         and (lhs.UDATA_INDEST == rhs.UDATA_INDEST)
         and (lhs.UDATA_OUTDEST == rhs.UDATA_OUTDEST)
         and (lhs.UDATA_XSPLIT == rhs.UDATA_XSPLIT)
         and (lhs.DDATA_INDEST == rhs.DDATA_INDEST)
         and (lhs.DDATA_OUTDEST == rhs.DDATA_OUTDEST)
         and (lhs.DDATA_XSPLIT == rhs.DDATA_XSPLIT) and (lhs.IDATA_DEST == rhs.IDATA_DEST)
         and (lhs.ODATA_DEST == rhs.ODATA_DEST);
}

/// Stream operator to print a nice message instead of bytes if a test fails
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

std::ostream& operator<<(std::ostream& out,
                         const BoutMeshExposer::ConnectionInfo& value) {
  return out << fmt::format("BoutMesh::ConnectionInfo{{"
                            "TS_up_in={}, "
                            "TS_up_out={}, "
                            "TS_down_in={}, "
                            "TS_down_out={}, "
                            "UDATA_INDEST={}, "
                            "UDATA_OUTDEST={}, "
                            "UDATA_XSPLIT={}, "
                            "DDATA_INDEST={}, "
                            "DDATA_OUTDEST={}, "
                            "DDATA_XSPLIT={}, "
                            "IDATA_DEST={}, "
                            "ODATA_DEST={}"
                            "}}",
                            value.TS_up_in, value.TS_up_out, value.TS_down_in,
                            value.TS_down_out, value.UDATA_INDEST, value.UDATA_OUTDEST,
                            value.UDATA_XSPLIT, value.DDATA_INDEST, value.DDATA_OUTDEST,
                            value.DDATA_XSPLIT, value.IDATA_DEST, value.ODATA_DEST);
}

////////////////////////////////////////////////////////////
// A bunch of functions for creating consistent configurations in
// different topologies. We don't just return a `BoutMeshExposer`,
// because we can reuse the `BoutMeshParameters` for other tests where
// we don't want a full `Mesh` object
BoutMeshParameters createCore(const BoutMeshGridInfo& grid) {
  return {grid,
          {grid.total_nx, grid.total_nx},
          {-1, (grid.total_ny / 2) - 1, (grid.total_ny / 2) - 1, grid.total_ny - 1,
           grid.total_ny / 2}};
}

BoutMeshParameters createSOL(const BoutMeshGridInfo& grid) {
  return {grid,
          {0, 0},
          {-1, (grid.total_ny / 2) - 1, (grid.total_ny / 2) - 1, grid.total_ny - 1,
           grid.total_ny / 2}};
}

BoutMeshParameters createLimiter(const BoutMeshGridInfo& grid) {
  return {grid,
          {grid.total_nx / 2, grid.total_nx},
          {-1, (grid.total_ny / 2) - 1, (grid.total_ny / 2) - 1, grid.total_ny - 1,
           grid.total_ny / 2}};
}

BoutMeshParameters createXPoint(const BoutMeshGridInfo& grid) {
  if (grid.nype < 4) {
    throw BoutException(
        "createXPoint: Not enough processors for x-point topology (nype={}, needs 4)",
        grid.nype);
  }

  return {grid,
          {grid.total_nx / 2, grid.total_nx / 2},
          {grid.local_ny - 1, grid.local_ny - 1, grid.total_ny - grid.local_ny - 1,
           grid.total_ny - grid.local_ny - 1, 2 * grid.local_ny}};
}

BoutMeshParameters createSingleNull(const BoutMeshGridInfo& grid) {
  if (grid.nype < 3) {
    throw BoutException(
        "createXPoint: Not enough processors for single-null topology (nype={}, needs 3)",
        grid.nype);
  }

  return {grid,
          {grid.total_nx / 2, grid.total_nx},
          {grid.local_ny - 1, (grid.total_ny / 2) - 1, (grid.total_ny / 2) - 1,
           grid.total_ny - grid.local_ny - 1, grid.total_ny / 2}};
}

BoutMeshParameters createDoubleNull(const BoutMeshGridInfo& grid) {
  if (grid.nype < 6) {
    throw BoutException("createDoubleNull: Not enough processors for double-null "
                        "topology (nype={}, needs 6)",
                        grid.nype);
  }

  const int ny_inner = 3 * grid.local_ny;
  return {grid,
          {grid.total_nx / 2, grid.total_nx / 2},
          {grid.local_ny - 1, ny_inner - grid.local_ny - 1, ny_inner + grid.local_ny - 1,
           grid.total_ny - grid.local_ny - 1, ny_inner}};
}

BoutMeshParameters createDisconnectedDoubleNull(const BoutMeshGridInfo& grid) {
  if (grid.nype < 6) {
    throw BoutException(
        "createDisconnectedDoubleNull: Not enough processors for disconnected "
        "double-null topology (nype={}, needs 6)",
        grid.nype);
  }

  if ((grid.total_nx / 2) + 4 > grid.total_nx) {
    throw BoutException(
        "createDisconnectedDoubleNull: Not enough points in x-direction "
        "(need ixseps2 = ((nxpe * (local_nx - 2)) + 2) / 2 + 4 = {} to "
        "be less than total_nx = (nxpe * (local_nx - 2)) + 2 = {}; nxpe={}, local_nx={}",
        (grid.total_nx / 2) + 4, grid.total_nx, grid.nxpe, grid.local_nx);
  }

  const int ny_inner = 3 * grid.local_ny;
  return {grid,
          {grid.total_nx / 2, grid.total_nx / 2 + 4},
          {grid.local_ny - 1, ny_inner - grid.local_ny - 1, ny_inner + grid.local_ny - 1,
           grid.total_ny - grid.local_ny - 1, ny_inner}};
}

////////////////////////////////////////////////////////////
// Start of tests

TEST(BoutMeshTest, NullOptionsCheck) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
}

// Not a great test as it's not specific to the thing we want to test,
// and can also take a whopping ~300ms!
TEST(BoutMeshTest, SingleCoreDecomposition) {
  WithQuietOutput debug{output_debug};
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  WithQuietOutput progress{output_progress};

  Options options{};
  options["ny"] = 1;
  options["nx"] = 4;
  options["nz"] = 1;
  options["MXG"] = 1;
  options["MYG"] = 0;

  bout::globals::mpi = new MpiWrapper();
  BoutMesh mesh{new GridFromOptions{&options}, &options};
  EXPECT_NO_THROW(mesh.load());
  delete bout::globals::mpi;
  bout::globals::mpi = nullptr;
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

struct BoutMeshSetYDecompositionTest
    : public ::testing::TestWithParam<SetYDecompositionTestParameters> {
  virtual ~BoutMeshSetYDecompositionTest() = default;
};

INSTANTIATE_TEST_SUITE_P(
    GoodDecompositions, BoutMeshSetYDecompositionTest,
    ::testing::Values(SetYDecompositionTestParameters{{-1, 7, 15, 23, 12},
                                                      {-1, 7, 15, 23, 12},
                                                      0,
                                                      "CoreOnly"},
                      SetYDecompositionTestParameters{
                          {3, 7, 7, 19, 12}, {3, 7, 7, 19, 12}, 1, "SingleNull"},
                      SetYDecompositionTestParameters{
                          {3, 7, 15, 19, 12}, {3, 7, 15, 19, 12}, 2, "DoubleNull"},
                      SetYDecompositionTestParameters{
                          {-12, 7, 15, 19, 12}, {-1, 7, 15, 19, 12}, 2, "Jyseps11Low"},
                      SetYDecompositionTestParameters{
                          {3, 1, 15, 19, 12}, {3, 4, 15, 19, 12}, 2, "Jyseps21Low"},
                      SetYDecompositionTestParameters{
                          {3, 7, 5, 19, 12}, {3, 7, 7, 19, 12}, 1, "Jyseps12Low"},
                      SetYDecompositionTestParameters{
                          {3, 7, 15, 32, 12}, {3, 7, 15, 23, 12}, 2, "Jyseps22High"},
                      SetYDecompositionTestParameters{
                          {3, 7, 15, 8, 12}, {3, 7, 15, 15, 12}, 2, "Jyseps22Low"}),
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

DecompositionTestParameters
makeDecompositionTestParameters(const BoutMeshParameters& inputs,
                                const std::string& name) {
  return {inputs.grid.total_processors,
          inputs.grid.nype,
          inputs.grid.total_ny,
          inputs.grid.num_y_guards,
          inputs.y_indices,
          "",
          name};
}

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
            2, 2, 8, 1, {-1, 4, 4, 7, 4}, "", "EightPointsTwoCoresNYPE2"},
        // The following should basically all work by construction
        makeDecompositionTestParameters(createCore({4, 4, 2, 2, 1, 1}), "Core"),
        makeDecompositionTestParameters(createSOL({4, 4, 2, 2, 1, 1}), "SOL"),
        makeDecompositionTestParameters(createLimiter({4, 4, 2, 2, 1, 1}), "Limiter"),
        makeDecompositionTestParameters(createXPoint({4, 4, 2, 2, 1, 4}), "XPoint"),
        makeDecompositionTestParameters(createSingleNull({4, 4, 2, 2, 1, 3}),
                                        "SingleNull"),
        makeDecompositionTestParameters(createDoubleNull({4, 4, 2, 2, 1, 6}),
                                        "DoubleNull"),
        makeDecompositionTestParameters(createDisconnectedDoubleNull({12, 4, 2, 2, 1, 6}),
                                        "DisconnectedDoubleNull")),
    DecompositionTestParametersToString);

TEST_P(BoutMeshDecompositionTest, CheckYDecomposition) {
  const auto params = GetParam();
  auto result = bout::checkBoutMeshYDecomposition(
      params.num_y_processors, params.ny, 1, params.indices.jyseps1_1,
      params.indices.jyseps2_1, params.indices.jyseps1_2, params.indices.jyseps2_2,
      params.indices.ny_inner);

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
      params.num_y_processors, params.ny, params.num_y_guards, params.indices.jyseps1_1,
      params.indices.jyseps2_1, params.indices.jyseps1_2, params.indices.jyseps2_2,
      params.indices.ny_inner);

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

FindProcessorParameters makeFindProcessorParameters(const BoutMeshParameters& inputs) {
  return {inputs.grid.total_processors,
          inputs.grid.total_nx,
          inputs.grid.total_ny,
          inputs.grid.num_x_guards,
          inputs.grid.num_x_guards,
          inputs.y_indices,
          inputs.grid.nxpe,
          inputs.grid.nype};
}

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
        FindProcessorParameters{8192, 132, 128, 2, 2, {15, 47, 79, 111, 64}, 128, 64},
        // The following should work basically by construction
        makeFindProcessorParameters(createCore({4, 4, 2, 2, 1, 1})),
        makeFindProcessorParameters(createCore({4, 4, 2, 2, 2, 2})),
        makeFindProcessorParameters(createCore({4, 4, 2, 2, 4, 4})),
        makeFindProcessorParameters(createSOL({4, 4, 2, 2, 1, 1})),
        makeFindProcessorParameters(createSOL({4, 4, 2, 2, 17, 13})),
        makeFindProcessorParameters(createSOL({4, 4, 2, 2, 37, 67})),
        makeFindProcessorParameters(createLimiter({4, 4, 2, 2, 1, 1})),
        makeFindProcessorParameters(createLimiter({4, 4, 2, 2, 5, 6})),
        makeFindProcessorParameters(createXPoint({4, 4, 2, 2, 1, 4})),
        makeFindProcessorParameters(createXPoint({4, 4, 2, 2, 89, 32})),
        makeFindProcessorParameters(createSingleNull({4, 4, 2, 2, 1, 3})),
        makeFindProcessorParameters(createSingleNull({4, 4, 2, 2, 23, 31})),
        makeFindProcessorParameters(createDoubleNull({4, 4, 2, 2, 1, 6})),
        makeFindProcessorParameters(createDoubleNull({4, 4, 2, 2, 7, 7})),
        makeFindProcessorParameters(createDisconnectedDoubleNull({12, 4, 2, 2, 1, 6})),
        makeFindProcessorParameters(createDisconnectedDoubleNull({12, 4, 2, 2, 6, 66}))));

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

TEST(BoutMeshTest, GlobalXIntSymmetricX) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh01(4, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.GlobalX(0), -0.125);
  EXPECT_EQ(mesh01.GlobalX(1), 0.125);
  EXPECT_EQ(mesh01.GlobalX(2), 0.375);
  EXPECT_EQ(mesh01.GlobalX(3), 0.625);
  EXPECT_EQ(mesh01.GlobalX(4), 0.875);
}

TEST(BoutMeshTest, GlobalXIntAsymmetricX) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh01(4, 3, 1, 2, 2, 0, 1, false, false);
  EXPECT_EQ(mesh01.GlobalX(0), 0.);
  EXPECT_EQ(mesh01.GlobalX(1), 0.25);
  EXPECT_EQ(mesh01.GlobalX(2), 0.5);
  EXPECT_EQ(mesh01.GlobalX(3), 0.75);
  EXPECT_EQ(mesh01.GlobalX(4), 1.0);
}

TEST(BoutMeshTest, GlobalXRealSymmetricX) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh01(4, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.GlobalX(0.5), 0.);
  EXPECT_EQ(mesh01.GlobalX(1.5), 0.25);
  EXPECT_EQ(mesh01.GlobalX(2.5), 0.5);
  EXPECT_EQ(mesh01.GlobalX(3.5), 0.75);
  EXPECT_EQ(mesh01.GlobalX(4.5), 1.0);
}

TEST(BoutMeshTest, GlobalXRealAsymmetricX) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh01(4, 3, 1, 2, 2, 0, 1, false, false);
  EXPECT_EQ(mesh01.GlobalX(0.5), 0.125);
  EXPECT_EQ(mesh01.GlobalX(1.5), 0.375);
  EXPECT_EQ(mesh01.GlobalX(2.5), 0.625);
  EXPECT_EQ(mesh01.GlobalX(3.5), 0.875);
  EXPECT_EQ(mesh01.GlobalX(4.5), 1.125);
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

TEST(BoutMeshTest, GlobalYIntSymmetricY) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh_inner_pf(createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 0}));
  EXPECT_EQ(mesh_inner_pf.GlobalY(0), -0.5625);
  EXPECT_EQ(mesh_inner_pf.GlobalY(1), -0.4375);
  EXPECT_EQ(mesh_inner_pf.GlobalY(2), -0.3125);
  EXPECT_EQ(mesh_inner_pf.GlobalY(3), -0.1875);

  BoutMeshExposer mesh_inner_core(
      createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 1}));
  EXPECT_EQ(mesh_inner_core.GlobalY(0), -0.0625);
  EXPECT_EQ(mesh_inner_core.GlobalY(1), 0.0625);
  EXPECT_EQ(mesh_inner_core.GlobalY(2), 0.1875);
  EXPECT_EQ(mesh_inner_core.GlobalY(3), 0.3125);

  BoutMeshExposer mesh_outer_core(
      createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 4}));
  EXPECT_EQ(mesh_outer_core.GlobalY(0), 0.4375);
  EXPECT_EQ(mesh_outer_core.GlobalY(1), 0.5625);
  EXPECT_EQ(mesh_outer_core.GlobalY(2), 0.6875);
  EXPECT_EQ(mesh_outer_core.GlobalY(3), 0.8125);

  BoutMeshExposer mesh_outer_pf(createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 5}));
  EXPECT_EQ(mesh_outer_pf.GlobalY(0), 0.9375);
  EXPECT_EQ(mesh_outer_pf.GlobalY(1), 1.0625);
  EXPECT_EQ(mesh_outer_pf.GlobalY(2), 1.1875);
  EXPECT_EQ(mesh_outer_pf.GlobalY(3), 1.3125);
}

TEST(BoutMeshTest, GlobalYIntAsymmetricY) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  auto grid_inner_pf = createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 0});
  grid_inner_pf.grid.symmetric_Y = false;
  BoutMeshExposer mesh_inner_pf(grid_inner_pf);
  EXPECT_EQ(mesh_inner_pf.GlobalY(0), 0);
  EXPECT_EQ(mesh_inner_pf.GlobalY(1), 0);
  EXPECT_EQ(mesh_inner_pf.GlobalY(2), 0);
  EXPECT_EQ(mesh_inner_pf.GlobalY(3), 0);

  auto grid_inner_core = createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 1});
  grid_inner_core.grid.symmetric_Y = false;
  BoutMeshExposer mesh_inner_core(grid_inner_core);
  EXPECT_EQ(mesh_inner_core.GlobalY(0), -0.125);
  EXPECT_EQ(mesh_inner_core.GlobalY(1), 0.0);
  EXPECT_EQ(mesh_inner_core.GlobalY(2), 0.125);
  EXPECT_EQ(mesh_inner_core.GlobalY(3), 0.25);

  auto grid_outer_core = createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 4});
  grid_outer_core.grid.symmetric_Y = false;
  BoutMeshExposer mesh_outer_core(grid_outer_core);
  EXPECT_EQ(mesh_outer_core.GlobalY(0), 0.375);
  EXPECT_EQ(mesh_outer_core.GlobalY(1), 0.5);
  EXPECT_EQ(mesh_outer_core.GlobalY(2), 0.625);
  EXPECT_EQ(mesh_outer_core.GlobalY(3), 0.75);

  auto grid_outer_pf = createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 5});
  grid_outer_pf.grid.symmetric_Y = false;
  BoutMeshExposer mesh_outer_pf(grid_outer_pf);
  // EXPECT_EQ(mesh_outer_pf.GlobalY(0), 2.375); // Does this make sense?
  EXPECT_EQ(mesh_outer_pf.GlobalY(1), 1);
  EXPECT_EQ(mesh_outer_pf.GlobalY(2), 1);
  EXPECT_EQ(mesh_outer_pf.GlobalY(3), 1);
}

TEST(BoutMeshTest, GlobalYRealSymmetricY) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh_inner_pf(createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 0}));
  EXPECT_EQ(mesh_inner_pf.GlobalY(0.5), -0.5);
  EXPECT_EQ(mesh_inner_pf.GlobalY(1.5), -0.375);
  EXPECT_EQ(mesh_inner_pf.GlobalY(2.5), -0.25);
  EXPECT_EQ(mesh_inner_pf.GlobalY(3.5), -0.125);

  BoutMeshExposer mesh_inner_core(
      createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 1}));
  EXPECT_EQ(mesh_inner_core.GlobalY(0.5), 0.0);
  EXPECT_EQ(mesh_inner_core.GlobalY(1.5), 0.125);
  EXPECT_EQ(mesh_inner_core.GlobalY(2.5), 0.25);
  EXPECT_EQ(mesh_inner_core.GlobalY(3.5), 0.375);

  BoutMeshExposer mesh_outer_core(
      createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 4}));
  EXPECT_EQ(mesh_outer_core.GlobalY(0.5), 0.5);
  EXPECT_EQ(mesh_outer_core.GlobalY(1.5), 0.625);
  EXPECT_EQ(mesh_outer_core.GlobalY(2.5), 0.75);
  EXPECT_EQ(mesh_outer_core.GlobalY(3.5), 0.875);

  BoutMeshExposer mesh_outer_pf(createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 5}));
  EXPECT_EQ(mesh_outer_pf.GlobalY(0.5), 1.0);
  EXPECT_EQ(mesh_outer_pf.GlobalY(1.5), 1.125);
  EXPECT_EQ(mesh_outer_pf.GlobalY(2.5), 1.25);
  EXPECT_EQ(mesh_outer_pf.GlobalY(3.5), 1.375);
}

TEST(BoutMeshTest, GlobalYRealAsymmetricY) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  auto grid_inner_pf = createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 0});
  grid_inner_pf.grid.symmetric_Y = false;
  BoutMeshExposer mesh_inner_pf(grid_inner_pf);
  EXPECT_EQ(mesh_inner_pf.GlobalY(0.5), 0);
  EXPECT_EQ(mesh_inner_pf.GlobalY(1.5), 0);
  EXPECT_EQ(mesh_inner_pf.GlobalY(2.5), 0);
  EXPECT_EQ(mesh_inner_pf.GlobalY(3.5), 0);

  auto grid_inner_core = createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 1});
  grid_inner_core.grid.symmetric_Y = false;
  BoutMeshExposer mesh_inner_core(grid_inner_core);
  EXPECT_EQ(mesh_inner_core.GlobalY(0.5), -0.0625);
  EXPECT_EQ(mesh_inner_core.GlobalY(1.5), 0.0625);
  EXPECT_EQ(mesh_inner_core.GlobalY(2.5), 0.1875);
  EXPECT_EQ(mesh_inner_core.GlobalY(3.5), 0.3125);

  auto grid_outer_core = createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 4});
  grid_outer_core.grid.symmetric_Y = false;
  BoutMeshExposer mesh_outer_core(grid_outer_core);
  EXPECT_EQ(mesh_outer_core.GlobalY(0.5), 0.4375);
  EXPECT_EQ(mesh_outer_core.GlobalY(1.5), 0.5625);
  EXPECT_EQ(mesh_outer_core.GlobalY(2.5), 0.6875);
  EXPECT_EQ(mesh_outer_core.GlobalY(3.5), 0.8125);

  auto grid_outer_pf = createDisconnectedDoubleNull({12, 4, 1, 1, 1, 6, 0, 5});
  grid_outer_pf.grid.symmetric_Y = false;
  BoutMeshExposer mesh_outer_pf(grid_outer_pf);
  EXPECT_EQ(mesh_outer_pf.GlobalY(0.5), 1);
  EXPECT_EQ(mesh_outer_pf.GlobalY(1.5), 1);
  EXPECT_EQ(mesh_outer_pf.GlobalY(2.5), 1);
  EXPECT_EQ(mesh_outer_pf.GlobalY(3.5), 1);
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

void checkRegionSizes(const BoutMeshExposer& mesh, std::array<int, 3> rgn_lower_y,
                      std::array<int, 3> rgn_upper_y, std::array<int, 2> rgn_x) {
  EXPECT_EQ(mesh.getRegion("RGN_LOWER_INNER_Y").size(), rgn_lower_y[0]);
  EXPECT_EQ(mesh.getRegion("RGN_LOWER_OUTER_Y").size(), rgn_lower_y[1]);
  EXPECT_EQ(mesh.getRegion("RGN_LOWER_Y").size(), rgn_lower_y[2]);

  EXPECT_EQ(mesh.getRegion("RGN_UPPER_INNER_Y").size(), rgn_upper_y[0]);
  EXPECT_EQ(mesh.getRegion("RGN_UPPER_OUTER_Y").size(), rgn_upper_y[1]);
  EXPECT_EQ(mesh.getRegion("RGN_UPPER_Y").size(), rgn_upper_y[2]);

  EXPECT_EQ(mesh.getRegion("RGN_INNER_X").size(), rgn_x[0]);
  EXPECT_EQ(mesh.getRegion("RGN_OUTER_X").size(), rgn_x[1]);
}

// These next few tests check both default_connections and the Region
// creation, as these are quite tightly linked.

TEST(BoutMeshTest, DefaultConnectionsCore1x1) {
  WithQuietOutput info{output_info};
  // 5x3x1 grid on 1 processor, 1 boundary point. Boundaries should be
  // simple 1D rectangles, with 4 boundaries on this processor
  BoutMeshExposer mesh00(5, 3, 1, 1, 1, 0, 0, false);

  mesh00.default_connections();

  BoutMeshExposer::ConnectionInfo expected{false, false, false, false, -1, -1,
                                           0,     -1,    -1,    0,     -1, -1};
  EXPECT_EQ(mesh00.getConnectionInfo(), expected);

  mesh00.createDefaultRegions();
  mesh00.addBoundaryRegions();

  SCOPED_TRACE("DefaultConnectionsCore1x1");
  checkRegionSizes(mesh00, {5, 0, 5}, {0, 5, 5}, {3, 3});
}

TEST(BoutMeshTest, TopologySOL2x2) {
  WithQuietOutput info{output_info};

  {
    SCOPED_TRACE("TopologySOL2x2, mesh00");
    BoutMeshExposer mesh00(createSOL({3, 3, 1, 1, 2, 2, 0, 0}));
    BoutMeshExposer::ConnectionInfo expected00{false, false, false, false, -1, 2,
                                               0,     -1,    -1,    0,     -1, 1};
    EXPECT_EQ(mesh00.getConnectionInfo(), expected00);
    checkRegionSizes(mesh00, {4, 0, 4}, {0, 0, 0}, {3, 0});
  }

  {
    SCOPED_TRACE("TopologySOL2x2, mesh01");
    BoutMeshExposer mesh01(createSOL({3, 3, 1, 1, 2, 2, 0, 1}));
    BoutMeshExposer::ConnectionInfo expected01{false, false, false, false, -1, -1,
                                               0,     -1,    0,     0,     -1, 3};
    EXPECT_EQ(mesh01.getConnectionInfo(), expected01);
    checkRegionSizes(mesh01, {0, 0, 0}, {0, 4, 4}, {3, 0});
  }

  {
    SCOPED_TRACE("TopologySOL2x2, mesh10");
    BoutMeshExposer mesh10(createSOL({3, 3, 1, 1, 2, 2, 1, 0}));
    BoutMeshExposer::ConnectionInfo expected10{false, false, false, false, -1, 3,
                                               0,     -1,    -1,    0,     0,  -1};
    EXPECT_EQ(mesh10.getConnectionInfo(), expected10);
    checkRegionSizes(mesh10, {4, 0, 4}, {0, 0, 0}, {0, 3});
  }

  {
    SCOPED_TRACE("TopologySOL2x2, mesh11");
    BoutMeshExposer mesh11(createSOL({3, 3, 1, 1, 2, 2, 1, 1}));
    BoutMeshExposer::ConnectionInfo expected11{false, false, false, false, -1, -1,
                                               0,     -1,    1,     0,     2,  -1};
    EXPECT_EQ(mesh11.getConnectionInfo(), expected11);
    checkRegionSizes(mesh11, {0, 0, 0}, {0, 4, 4}, {0, 3});
  }
}

TEST(BoutMeshTest, TopologySOLPeriodicX2x2) {
  WithQuietOutput info{output_info};

  {
    SCOPED_TRACE("TopologySOLPeriodicX2x2, mesh00");

    BoutMeshExposer mesh00(createSOL({3, 3, 1, 1, 2, 2, 0, 0}), true);
    BoutMeshExposer::ConnectionInfo expected00{false, false, false, false, -1, 2,
                                               0,     -1,    -1,    0,     1,  1};
    EXPECT_EQ(mesh00.getConnectionInfo(), expected00);
    checkRegionSizes(mesh00, {4, 0, 4}, {0, 0, 0}, {0, 0});
  }

  {
    SCOPED_TRACE("TopologySOLPeriodicX2x2, mesh01");
    BoutMeshExposer mesh01(createSOL({3, 3, 1, 1, 2, 2, 0, 1}), true);
    BoutMeshExposer::ConnectionInfo expected01{false, false, false, false, -1, -1,
                                               0,     -1,    0,     0,     3,  3};
    EXPECT_EQ(mesh01.getConnectionInfo(), expected01);
    checkRegionSizes(mesh01, {0, 0, 0}, {0, 4, 4}, {0, 0});
  }

  {
    SCOPED_TRACE("TopologySOLPeriodicX2x2, mesh10");
    BoutMeshExposer mesh10(createSOL({3, 3, 1, 1, 2, 2, 1, 0}), true);
    BoutMeshExposer::ConnectionInfo expected10{false, false, false, false, -1, 3,
                                               0,     -1,    -1,    0,     0,  0};
    EXPECT_EQ(mesh10.getConnectionInfo(), expected10);
    checkRegionSizes(mesh10, {4, 0, 4}, {0, 0, 0}, {0, 0});
  }

  {
    SCOPED_TRACE("TopologySOLPeriodicX2x2, mesh11");
    BoutMeshExposer mesh11(createSOL({3, 3, 1, 1, 2, 2, 1, 1}), true);
    BoutMeshExposer::ConnectionInfo expected11{false, false, false, false, -1, -1,
                                               0,     -1,    1,     0,     2,  2};
    EXPECT_EQ(mesh11.getConnectionInfo(), expected11);
    checkRegionSizes(mesh11, {0, 0, 0}, {0, 4, 4}, {0, 0});
  }
}

TEST(BoutMeshTest, TopologySingleNull2x3) {
  WithQuietOutput info{output_info};

  {
    SCOPED_TRACE("TopologySingleNull2x3, mesh00");
    BoutMeshExposer mesh00(createSingleNull({3, 3, 1, 1, 2, 3, 0, 0}));
    BoutMeshExposer::ConnectionInfo expected00{false, false, false, false, 4,  2,
                                               4,     -1,    -1,    0,     -1, 1};
    EXPECT_EQ(mesh00.getConnectionInfo(), expected00);
    checkRegionSizes(mesh00, {4, 0, 4}, {0, 0, 0}, {3, 0});
  }

  {
    SCOPED_TRACE("TopologySingleNull2x3, mesh01");
    BoutMeshExposer mesh01(createSingleNull({3, 3, 1, 1, 2, 3, 0, 1}));
    BoutMeshExposer::ConnectionInfo expected01{true, false, true, false, 2,  4,
                                               4,    2,     0,    4,     -1, 3};
    EXPECT_EQ(mesh01.getConnectionInfo(), expected01);
    checkRegionSizes(mesh01, {0, 0, 0}, {0, 0, 0}, {3, 0});
  }

  {
    SCOPED_TRACE("TopologySingleNull2x3, mesh02");
    BoutMeshExposer mesh02(createSingleNull({3, 3, 1, 1, 2, 3, 0, 2}));
    BoutMeshExposer::ConnectionInfo expected02{false, false, false, false, -1, -1,
                                               0,     0,     2,     4,     -1, 5};
    EXPECT_EQ(mesh02.getConnectionInfo(), expected02);
    checkRegionSizes(mesh02, {0, 0, 0}, {0, 4, 4}, {3, 0});
  }

  {
    SCOPED_TRACE("TopologySingleNull2x3, mesh10");
    BoutMeshExposer mesh10(createSingleNull({3, 3, 1, 1, 2, 3, 1, 0}));
    BoutMeshExposer::ConnectionInfo expected10{false, false, false, false, 5, 3,
                                               1,     -1,    -1,    0,     0, -1};
    EXPECT_EQ(mesh10.getConnectionInfo(), expected10);
    checkRegionSizes(mesh10, {4, 0, 4}, {0, 0, 0}, {0, 3});
  }

  {
    SCOPED_TRACE("TopologySingleNull2x3, mesh11");
    BoutMeshExposer mesh11(createSingleNull({3, 3, 1, 1, 2, 3, 1, 1}));
    BoutMeshExposer::ConnectionInfo expected11{true, false, true, false, 3, 5,
                                               1,    3,     1,    1,     2, -1};
    EXPECT_EQ(mesh11.getConnectionInfo(), expected11);
    checkRegionSizes(mesh11, {0, 0, 0}, {0, 0, 0}, {0, 3});
  }

  {
    SCOPED_TRACE("TopologySingleNull2x3, mesh12");
    BoutMeshExposer mesh12(createSingleNull({3, 3, 1, 1, 2, 3, 1, 2}));
    BoutMeshExposer::ConnectionInfo expected11{false, false, false, false, -1, -1,
                                               0,     1,     3,     1,     4,  -1};
    EXPECT_EQ(mesh12.getConnectionInfo(), expected11);
    checkRegionSizes(mesh12, {0, 0, 0}, {0, 4, 4}, {0, 3});
  }
}

TEST(BoutMeshTest, TopologyDisconnectedDoubleNull1x6) {
  WithQuietOutput info{output_info};

  {
    SCOPED_TRACE("TopologyDisconnectedDoubleNull1x6, mesh00"); // Inner lower leg
    BoutMeshExposer mesh00(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 0}));
    BoutMeshExposer::ConnectionInfo expected00{false, false, false, false, 5,  1,
                                               7,     -1,    -1,    0,     -1, -1};
    EXPECT_EQ(mesh00.getConnectionInfo(), expected00);
    checkRegionSizes(mesh00, {14, 0, 14}, {0, 0, 0}, {3, 3});
  }

  {
    SCOPED_TRACE("TopologyDisconnectedDoubleNull1x6, mesh01"); // Inner core
    BoutMeshExposer mesh01(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
    BoutMeshExposer::ConnectionInfo expected01{false, false, true, false, 4,  2,
                                               11,    4,     0,    7,     -1, -1};
    EXPECT_EQ(mesh01.getConnectionInfo(), expected01);
    checkRegionSizes(mesh01, {0, 0, 0}, {0, 0, 0}, {3, 3});
  }

  {
    SCOPED_TRACE("TopologyDisconnectedDoubleNull1x6, mesh02"); // Inner upper leg
    BoutMeshExposer mesh02(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 2}));
    BoutMeshExposer::ConnectionInfo expected01{false, false, false, false, -1, -1,
                                               14,    3,     1,     11,    -1, -1};
    EXPECT_EQ(mesh02.getConnectionInfo(), expected01);
    checkRegionSizes(mesh02, {0, 0, 0}, {14, 0, 14}, {3, 3});
  }

  {
    SCOPED_TRACE("TopologyDisconnectedDoubleNull1x6, mesh03"); // Outer upper leg
    BoutMeshExposer mesh03(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 3}));
    BoutMeshExposer::ConnectionInfo expected10{false, false, false, false, 2,  4,
                                               11,    -1,    -1,    14,    -1, -1};
    EXPECT_EQ(mesh03.getConnectionInfo(), expected10);
    checkRegionSizes(mesh03, {0, 14, 14}, {0, 0, 0}, {3, 3});
  }

  {
    SCOPED_TRACE("TopologyDisconnectedDoubleNull1x6, mesh04"); // Outer core
    BoutMeshExposer mesh04(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 4}));
    BoutMeshExposer::ConnectionInfo expected11{true, false, false, false, 1,  5,
                                               7,    1,     3,     11,    -1, -1};
    EXPECT_EQ(mesh04.getConnectionInfo(), expected11);
    checkRegionSizes(mesh04, {0, 0, 0}, {0, 0, 0}, {3, 3});
  }

  {
    SCOPED_TRACE("TopologyDisconnectedDoubleNull1x6, mesh05"); // Outer lower leg
    BoutMeshExposer mesh05(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 5}));
    BoutMeshExposer::ConnectionInfo expected11{false, false, false, false, -1, -1,
                                               0,     0,     4,     7,     -1, -1};
    EXPECT_EQ(mesh05.getConnectionInfo(), expected11);
    checkRegionSizes(mesh05, {0, 0, 0}, {0, 14, 14}, {3, 3});
  }
}

TEST(BoutMeshTest, SetDerivedGridSizes) {
  WithQuietOutput info{output_info};
  BoutMeshGridInfo grid{12, 3, 1, 2, 3, 6, 2, 2};
  BoutMeshExposer mesh(createDisconnectedDoubleNull(grid));

  EXPECT_EQ(mesh.GlobalNx, grid.total_nx);
  EXPECT_EQ(mesh.GlobalNy, grid.total_ny + 8);
  EXPECT_EQ(mesh.GlobalNz, 1);

  EXPECT_EQ(mesh.GlobalNxNoBoundaries, grid.total_nx - 2);
  EXPECT_EQ(mesh.GlobalNyNoBoundaries, grid.total_ny);
  EXPECT_EQ(mesh.GlobalNzNoBoundaries, 1);

  EXPECT_EQ(mesh.OffsetX, 2 * grid.local_nx);
  EXPECT_EQ(mesh.OffsetY, 2 * grid.local_ny);
  EXPECT_EQ(mesh.OffsetZ, 0);

  EXPECT_EQ(mesh.LocalNx, grid.local_nx + 2);
  EXPECT_EQ(mesh.LocalNy, grid.local_ny + 4);
  EXPECT_EQ(mesh.LocalNz, 1);

  EXPECT_EQ(mesh.xstart, 1);
  EXPECT_EQ(mesh.xend, 12);
  EXPECT_EQ(mesh.ystart, 2);
  EXPECT_EQ(mesh.yend, 4);
  EXPECT_EQ(mesh.zstart, 0);
  EXPECT_EQ(mesh.zend, 0);
}

TEST(BoutMeshTest, CreateXBoundariesPeriodicX) {
  WithQuietOutput info{output_info};
  // Periodic in X, so no boundaries
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 3, 6, 1, 0}));
  mesh.periodicX = true;
  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_TRUE(boundaries.empty());
}

TEST(BoutMeshTest, CreateXBoundariesNoGuards) {
  WithQuietOutput info{output_info};
  // No guards in X, so no boundaries
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 0, 1, 3, 6, 1, 0}));
  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_TRUE(boundaries.empty());
}

TEST(BoutMeshTest, CreateXBoundariesDoubleNullInsidePF) {
  WithQuietOutput info{output_info};
  // Three cores in X, inside core, one boundary
  BoutMeshExposer mesh_inside(createDisconnectedDoubleNull({12, 3, 1, 1, 3, 6, 0, 0}));
  mesh_inside.createXBoundaries();

  auto boundaries_inside = mesh_inside.getBoundaries();
  EXPECT_EQ(boundaries_inside.size(), 1);
  EXPECT_EQ(boundaries_inside[0]->label, "pf");
}

TEST(BoutMeshTest, CreateXBoundariesDoubleNullMiddlePF) {
  WithQuietOutput info{output_info};
  // Three cores in X, middle core, so no boundaries
  BoutMeshExposer mesh_middle(createDisconnectedDoubleNull({12, 3, 1, 1, 3, 6, 1, 0}));
  mesh_middle.createXBoundaries();

  auto boundaries_middle = mesh_middle.getBoundaries();
  EXPECT_TRUE(boundaries_middle.empty());
}

TEST(BoutMeshTest, CreateXBoundariesDoubleNullOutsidePF) {
  WithQuietOutput info{output_info};
  // Three cores in X, outside core, one boundary
  BoutMeshExposer mesh_inside(createDisconnectedDoubleNull({12, 3, 1, 1, 3, 6, 0, 0}));
  mesh_inside.createXBoundaries();

  auto boundaries_inside = mesh_inside.getBoundaries();
  EXPECT_EQ(boundaries_inside.size(), 1);
  EXPECT_EQ(boundaries_inside[0]->label, "pf");
}

TEST(BoutMeshTest, CreateXBoundariesDoubleNullInsideOutsideCore) {
  WithQuietOutput info{output_info};
  // One core in X, so we expect two boundaries
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 2);
  EXPECT_EQ(boundaries[0]->label, "core");
  EXPECT_EQ(boundaries[1]->label, "sol");
}

TEST(BoutMeshTest, CreateYBoundariesNoGuards) {
  WithQuietOutput info{output_info};

  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 0, 1, 6, 0, 0}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_TRUE(boundaries.empty());
}

TEST(BoutMeshTest, CreateYBoundariesClosedFieldLines) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh(createCore({4, 4, 2, 2, 4, 4}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_TRUE(boundaries.empty());
}

TEST(BoutMeshTest, CreateYBoundariesInnerLower) {
  WithQuietOutput info{output_info};

  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 0}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "lower_target");
}

TEST(BoutMeshTest, CreateYBoundariesInnerUpper) {
  WithQuietOutput info{output_info};

  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 2}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "upper_target");
}

TEST(BoutMeshTest, CreateYBoundariesOuterUpper) {
  WithQuietOutput info{output_info};

  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 5}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "upper_target");
}

TEST(BoutMeshTest, CreateYBoundariesOuterLower) {
  WithQuietOutput info{output_info};

  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 3}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "lower_target");
}

TEST(BoutMestTest, PeriodicY) {
  WithQuietOutput info{output_info};

  BoutMeshExposer mesh00(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 0}));
  EXPECT_FALSE(mesh00.periodicY(2));
  EXPECT_FALSE(mesh00.periodicY(10));

  BoutMeshExposer mesh01(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
  EXPECT_TRUE(mesh01.periodicY(2));
  EXPECT_FALSE(mesh01.periodicY(10));
}

TEST(BoutMestTest, PeriodicYWithShiftAngle) {
  WithQuietOutput info{output_info};

  const std::vector<BoutReal> shift_angle = {-1., 11., 10., 9., 8., 7., 6.,
                                             5.,  4.,  3.,  2., 1., 0., -1.};

  BoutMeshExposer mesh00(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 0}));
  mesh00.setShiftAngle(shift_angle);
  BoutReal twist_shift00;
  EXPECT_FALSE(mesh00.periodicY(2, twist_shift00));
  EXPECT_EQ(twist_shift00, 0.);
  EXPECT_FALSE(mesh00.periodicY(10, twist_shift00));
  EXPECT_EQ(twist_shift00, 0.);

  BoutMeshExposer mesh01(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
  mesh01.setShiftAngle(shift_angle);
  BoutReal twist_shift01;
  EXPECT_TRUE(mesh01.periodicY(2, twist_shift01));
  EXPECT_EQ(twist_shift01, 10.);
  EXPECT_FALSE(mesh01.periodicY(10, twist_shift01));
  EXPECT_EQ(twist_shift01, 0.);
}

TEST(BoutMeshTest, NumberOfYBoundaries) {
  WithQuietOutput info{output_info};

  BoutMeshExposer mesh_SOL(createSOL({3, 3, 1, 1, 2, 2, 1, 1}));
  EXPECT_EQ(mesh_SOL.numberOfYBoundaries(), 1);

  BoutMeshExposer mesh_DND(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 0}));
  EXPECT_EQ(mesh_DND.numberOfYBoundaries(), 2);
}

TEST(BoutMeshTest, HasBranchCutLower) {
  WithQuietOutput info{output_info};

  BoutMeshExposer mesh_SOL(createSOL({3, 3, 1, 1, 2, 2, 1, 1}));
  EXPECT_EQ(mesh_SOL.hasBranchCutLower(2), std::make_pair(false, 0.));

  const std::vector<BoutReal> shift_angle = {-1., 11., 10., 9., 8., 7., 6.,
                                             5.,  4.,  3.,  2., 1., 0., -1.};
  BoutMeshExposer mesh_DND01(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
  mesh_DND01.setShiftAngle(shift_angle);
  EXPECT_EQ(mesh_DND01.hasBranchCutLower(3), std::make_pair(true, 9.));

  BoutMeshExposer mesh_DND04(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 4}));
  mesh_DND04.setShiftAngle(shift_angle);
  EXPECT_EQ(mesh_DND04.hasBranchCutLower(2), std::make_pair(false, 0.));
}

TEST(BoutMeshTest, HasBranchCutUpper) {
  WithQuietOutput info{output_info};

  BoutMeshExposer mesh_SOL(createSOL({3, 3, 1, 1, 2, 2, 1, 1}));
  EXPECT_EQ(mesh_SOL.hasBranchCutUpper(2), std::make_pair(false, 0.));

  const std::vector<BoutReal> shift_angle = {-1., 11., 10., 9., 8., 7., 6.,
                                             5.,  4.,  3.,  2., 1., 0., -1.};
  BoutMeshExposer mesh_DND01(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
  mesh_DND01.setShiftAngle(shift_angle);
  EXPECT_EQ(mesh_DND01.hasBranchCutUpper(3), std::make_pair(false, 0.));

  BoutMeshExposer mesh_DND04(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 4}));
  mesh_DND04.setShiftAngle(shift_angle);
  EXPECT_EQ(mesh_DND04.hasBranchCutUpper(2), std::make_pair(true, 10.));
}

TEST(BoutMeshTest, GetPossibleBoundariesCore) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh_core_1x1(createCore({12, 3, 1, 1, 1, 1, 0, 0}));
  BoutMeshExposer mesh_core_32x64(createCore({12, 3, 1, 1, 32, 64, 7, 4}));

  std::set<std::string> boundaries{"core", "sol"};

  EXPECT_EQ(mesh_core_1x1.getPossibleBoundaries(), boundaries);
  EXPECT_EQ(mesh_core_32x64.getPossibleBoundaries(), boundaries);
}

TEST(BoutMeshTest, GetPossibleBoundariesCorePeriodicX) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh_core_1x1(createCore({12, 3, 1, 1, 1, 1, 0, 0}), true);
  BoutMeshExposer mesh_core_32x64(createCore({12, 3, 1, 1, 32, 64, 7, 4}), true);

  EXPECT_TRUE(mesh_core_1x1.getPossibleBoundaries().empty());
  EXPECT_TRUE(mesh_core_32x64.getPossibleBoundaries().empty());
}

TEST(BoutMeshTest, GetPossibleBoundariesDND) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  BoutMeshExposer mesh_DND_1x6(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
  BoutMeshExposer mesh_DND_32x64(createDisconnectedDoubleNull({12, 3, 1, 1, 32, 64, 0, 4}));

  std::set<std::string> boundaries{"core", "pf", "sol", "upper_target", "lower_target"};

  EXPECT_EQ(mesh_DND_1x6.getPossibleBoundaries(), boundaries);
  EXPECT_EQ(mesh_DND_32x64.getPossibleBoundaries(), boundaries);
}
