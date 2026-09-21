#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "bout/griddata.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"

#include "fake_mesh.hxx"

#include <array>
#include <cstdio>
#include <iostream>
#include <ostream>
#include <sstream>
#include <streambuf>
#include <string>

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
  /// Construct from a grid data source, for testing the `Mesh::get`-based
  /// readers
  explicit BoutMeshExposer(GridDataSource* source, Options* options = nullptr)
      : BoutMesh(source, options) {}
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
  using BoutMesh::getMeshTopology;
  using BoutMesh::IngridTopology;
  using BoutMesh::mesh_topology;
  using BoutMesh::snowflake_type;
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

BoutMeshParameters createSnowflake(const BoutMeshGridInfo& grid) {
  // Need at least 6 y-subdomains for a minimal snowflake
  if (grid.nype < 6) {
    throw BoutException(
        "createSnowflake: Not enough processors for snowflake topology "
        "(nype={}, needs at least 6)",
        grid.nype);
  }

  if ((grid.total_nx / 2) + 4 > grid.total_nx) {
  throw BoutException(
      "createSnowflake: Not enough points in x-direction "
      "(need ixseps2 = ((nxpe * (local_nx - 2)) + 2) / 2 + 4 = {} to "
      "be less than total_nx = (nxpe * (local_nx - 2)) + 2 = {}; nxpe={}, local_nx={}",
      (grid.total_nx / 2) + 4, grid.total_nx, grid.nxpe, grid.local_nx);
}

  const int ny_inner = 4 * grid.local_ny;
  // Separatrix indices
  const int jyseps1_1 = grid.local_ny - 1;
  const int jyseps2_1 = ny_inner - 2 * grid.local_ny - 1;
  const int jyseps1_2 = ny_inner - grid.local_ny - 1;
  const int jyseps2_2 = grid.total_ny - grid.local_ny - 1;

  return {
      grid,
      // X separatrices (same as standard snowflake assumption)
      {grid.total_nx / 2, grid.total_nx / 2 + 4},
      // Y separatrices + ny_inner
      {jyseps1_1,
       jyseps2_1,
       jyseps1_2,
       jyseps2_2,
       ny_inner}
  };
}


////////////////////////////////////////////////////////////
// Start of tests

struct BoutMeshTest : public ::testing::Test {
  WithQuietOutput debug{output_debug};
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  WithQuietOutput progress{output_progress};
};

TEST_F(BoutMeshTest, NullOptionsCheck) {
  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
}

// Not a great test as it's not specific to the thing we want to test,
// and can also take a whopping ~300ms!
TEST_F(BoutMeshTest, SingleCoreDecomposition) {
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

TEST_F(BoutMeshTest, SetYDecompositionIndicesJyseps22LowInconsistent) {
  BoutMeshExposer mesh(1, 24, 1, 1, 1);

  EXPECT_THROW(mesh.setYDecompositionIndices({3, 7, 32, 8, 12}), BoutException);
}

//New bit: 

struct DecompositionTestParameters {
  int total_processors;
  int num_y_processors;
  int ny;
  int num_y_guards;
  BoutMeshExposer::YDecompositionIndices indices;
  std::string expected_message; // Expect this fragment to be in the result.reason for bad
                                // decompositions
  std::string name;
  MeshTopology mesh_topology; // New: topology enum
};

DecompositionTestParameters
makeDecompositionTestParameters(const BoutMeshParameters& inputs,
                                const std::string& name,
                                MeshTopology mesh_topology = MeshTopology::snowflake) { // default to snowflake
  return {inputs.grid.total_processors,
          inputs.grid.nype,
          inputs.grid.total_ny,
          inputs.grid.num_y_guards,
          inputs.y_indices,
          "",
          name,
          mesh_topology};
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
             "topology = {}, "
             "expected_message = {} }}",
             value.total_processors, value.num_y_processors, value.ny, value.num_y_guards,
             value.indices.jyseps1_1, value.indices.jyseps2_1, value.indices.jyseps1_2,
             value.indices.jyseps2_2, value.indices.ny_inner, toString(value.mesh_topology),
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
        DecompositionTestParameters{1, 1, 1, 1, {-1, 0, 0, 0, 0}, "", "OnePoint", MeshTopology::single_null},
        DecompositionTestParameters{1, 1, 8, 1, {-1, 4, 4, 7, 4}, "", "EightPoints", MeshTopology::single_null},
        DecompositionTestParameters{
            2, 1, 8, 1, {-1, 4, 4, 7, 4}, "", "EightPointsTwoCores", MeshTopology::single_null},
        DecompositionTestParameters{
            2, 2, 8, 1, {-1, 4, 4, 7, 4}, "", "EightPointsTwoCoresNYPE2", MeshTopology::single_null},
        // The following should basically all work by construction
        makeDecompositionTestParameters(createCore({4, 4, 2, 2, 1, 1}), "Core", MeshTopology::single_null),
        makeDecompositionTestParameters(createSOL({4, 4, 2, 2, 1, 1}), "SOL", MeshTopology::single_null),
        makeDecompositionTestParameters(createLimiter({4, 4, 2, 2, 1, 1}), "Limiter", MeshTopology::single_null),
        makeDecompositionTestParameters(createXPoint({4, 4, 2, 2, 1, 4}), "XPoint", MeshTopology::single_null),
        makeDecompositionTestParameters(createSingleNull({4, 4, 2, 2, 1, 3}),
                                        "SingleNull", MeshTopology::single_null),
        makeDecompositionTestParameters(createDoubleNull({4, 4, 2, 2, 1, 6}),
                                        "DoubleNull", MeshTopology::connected_double_null),
        makeDecompositionTestParameters(createDisconnectedDoubleNull({12, 4, 2, 2, 1, 6}),
                                        "DisconnectedDoubleNull", MeshTopology::unconnected_double_null)),
    DecompositionTestParametersToString);


//New tests for snowflake topology
INSTANTIATE_TEST_SUITE_P(
    GoodSnowflake, BoutMeshDecompositionTest,
    ::testing::Values(
        // Snowflake with 6 procesors in y, which is the minimum for a snowflake with the standard assumptions about where the separatrices are.
        DecompositionTestParameters{6, 6, 48, 1, {-1, 15, 23, 31, 24}, "", "SF48PointsNYPE6", MeshTopology::snowflake},

        // A slightly more realistic snowflake grid
        makeDecompositionTestParameters(createSnowflake({4, 6, 2, 2, 1, 9}), 
                                            "Snowflake", MeshTopology::snowflake)
    ),
    DecompositionTestParametersToString);



TEST_P(BoutMeshDecompositionTest, CheckYDecomposition) {
  const auto params = GetParam();
  auto result = bout::checkBoutMeshYDecomposition(
      params.num_y_processors, params.ny, 1,
      params.indices.jyseps1_1, params.indices.jyseps2_1,
      params.indices.jyseps1_2, params.indices.jyseps2_2,
      params.indices.ny_inner,
      params.mesh_topology); // <- pass topology

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
            1, 1, 4, 1, {3, 5, 6, 10, 0}, "Core region jyseps2_1", "CoreRegion1", MeshTopology::unconnected_double_null},
        DecompositionTestParameters{
            1, 1, 4, 1, {3, 7, 8, 11, 0}, "Core region jyseps2_2", "CoreRegion2", MeshTopology::unconnected_double_null},
        DecompositionTestParameters{
            1, 1, 4, 1, {3, 7, 8, 12, 11}, "leg region ny_inner", "UpperLeg1", MeshTopology::unconnected_double_null},
        DecompositionTestParameters{
            1, 1, 4, 1, {3, 7, 8, 12, 8}, "leg region jyseps1_2-ny_inner+1", "UpperLeg2", MeshTopology::unconnected_double_null},
        DecompositionTestParameters{
            1, 6, 25, 1, {3, 7, 15, 19, 12}, "leg region ny-jyseps2_2-1", "LegRegion", MeshTopology::unconnected_double_null}),
    DecompositionTestParametersToString);

INSTANTIATE_TEST_SUITE_P(
    BadSingleNull, BadBoutMeshDecompositionTest,
    ::testing::Values(
        DecompositionTestParameters{
            1, 1, 4, 1, {3, 4, 4, 6, 0}, "Core region jyseps2_2-jyseps1_1", "CoreRegion"},
        DecompositionTestParameters{
            1, 3, 13, 1, {3, 4, 4, 7, 0}, "leg region ny-jyseps2_2-1", "LegRegion"}),
    DecompositionTestParametersToString);


// New bad snowflake tests
INSTANTIATE_TEST_SUITE_P(
    BadSnowflake, BadBoutMeshDecompositionTest,
    ::testing::Values(
        // Core region
        DecompositionTestParameters{
          2, 2, 16, 1, {7, 14, 6, 9, 7}, "Core region jyseps2_1", "SF_CoreRegion1", MeshTopology::snowflake},

        // E-leg region
        DecompositionTestParameters{
          6, 6, 48, 1, {7, 15, 19, 25, 23}, "leg region jyseps2_2", "SF_ELeg1", MeshTopology::snowflake},
        DecompositionTestParameters{
          8, 8, 32, 1, {3, 7, 9, 15, 12}, "leg region ny_inner - 1", "SF_ELeg2", MeshTopology::snowflake},
        // W-leg region
        DecompositionTestParameters{
          6, 6, 48, 1, {7, 15, 24, 40, 33}, "leg region jyseps1_2", "SF_WLeg1", MeshTopology::snowflake}


        // Central region violation not possible if the others are correct
        // South-leg region violation not possible if the others are correct
    ),
    DecompositionTestParametersToString);



TEST_P(BadBoutMeshDecompositionTest, BadSingleCoreYDecomposition) {
  const auto params = GetParam();
  auto result = bout::checkBoutMeshYDecomposition(
      params.num_y_processors, params.ny, params.num_y_guards,
      params.indices.jyseps1_1, params.indices.jyseps2_1,
      params.indices.jyseps1_2, params.indices.jyseps2_2,
      params.indices.ny_inner,
      params.mesh_topology); // <- pass topology

  using ::testing::HasSubstr;

  EXPECT_FALSE(result.success);
  //Ask Peter about baddecomtest
  EXPECT_THAT(result.reason, HasSubstr(params.expected_message));
}

////////////////////////////////////////////////////////////
// Y decomposition for each member of the snowflake family

/// Y indices laid out so that every region of an SF+ topology is exactly
/// `mysub` points long, with `nype` = 6 and `ny` = 6 * `mysub`.
///
///   [0, jyseps1_1]            W target
///   (jyseps1_1, jyseps2_1]    core / W PFR middle
///   (jyseps2_1, jyseps1_2]    core / W PFR middle
///   (jyseps1_2, ny_inner-1]   E PFR above the E target
///   [ny_inner, jyseps2_2]     E PFR below the E target
///   (jyseps2_2, ny-1]         S target
BoutMeshExposer::YDecompositionIndices snowflakePlusIndices(int mysub) {
  return {mysub - 1, (2 * mysub) - 1, (3 * mysub) - 1, (5 * mysub) - 1, 4 * mysub};
}

/// As `snowflakePlusIndices`, but for SF-: the second X-point is in the SOL, so
/// `ny_inner` sits between `jyseps2_1` and `jyseps1_2` instead.
BoutMeshExposer::YDecompositionIndices snowflakeMinusIndices(int mysub) {
  return {mysub - 1, (2 * mysub) - 1, (4 * mysub) - 1, (5 * mysub) - 1, 3 * mysub};
}

struct SnowflakeDecompositionParameters {
  SnowflakeType snowflake_type;
  BoutMeshExposer::YDecompositionIndices indices;
  std::string test_name;
};

std::ostream& operator<<(std::ostream& out,
                         const SnowflakeDecompositionParameters& value) {
  return out << "SnowflakeDecompositionParameters{snowflake_type="
             << toString(value.snowflake_type) << ", indices=" << value.indices << "}";
}

std::string SnowflakeDecompositionParametersToString(
    const ::testing::TestParamInfo<SnowflakeDecompositionParameters>& param) {
  return param.param.test_name;
}

struct SnowflakeFamilyDecompositionTest
    : public ::testing::TestWithParam<SnowflakeDecompositionParameters> {};

/// `mysub` and `nype` used by every case below
constexpr int snowflake_mysub = 4;
constexpr int snowflake_nype = 6;
constexpr int snowflake_ny = snowflake_mysub * snowflake_nype;

INSTANTIATE_TEST_SUITE_P(
    EveryRegionExactlyOneProcessor, SnowflakeFamilyDecompositionTest,
    ::testing::Values(
        SnowflakeDecompositionParameters{SnowflakeType::SF_plus_low_field_side,
                                         snowflakePlusIndices(snowflake_mysub),
                                         "SFplusLFS"},
        SnowflakeDecompositionParameters{SnowflakeType::SF_plus_high_field_side,
                                         snowflakePlusIndices(snowflake_mysub),
                                         "SFplusHFS"},
        SnowflakeDecompositionParameters{SnowflakeType::SF_minus_low_field_side,
                                         snowflakeMinusIndices(snowflake_mysub),
                                         "SFminusLFS"},
        SnowflakeDecompositionParameters{SnowflakeType::SF_minus_high_field_side,
                                         snowflakeMinusIndices(snowflake_mysub),
                                         "SFminusHFS"},
        SnowflakeDecompositionParameters{SnowflakeType::SF,
                                         snowflakePlusIndices(snowflake_mysub),
                                         "GenericSF"}),
    SnowflakeDecompositionParametersToString);

/// Every branch cut in these layouts falls exactly on a processor boundary, so
/// the decomposition has to be accepted whichever family member it is.
TEST_P(SnowflakeFamilyDecompositionTest, EveryRegionIsAWholeNumberOfProcessors) {
  const auto params = GetParam();

  const auto result = bout::checkBoutMeshYDecomposition(
      snowflake_nype, snowflake_ny, 1, params.indices.jyseps1_1,
      params.indices.jyseps2_1, params.indices.jyseps1_2, params.indices.jyseps2_2,
      params.indices.ny_inner, MeshTopology::snowflake, params.snowflake_type);

  EXPECT_TRUE(result.success) << result.reason;
  EXPECT_TRUE(result.reason.empty()) << result.reason;
}

/// Shifting a single branch cut off a processor boundary has to be rejected.
/// `jyseps2_2` bounds a leg region in every member of the family.
TEST_P(SnowflakeFamilyDecompositionTest, RejectsBranchCutOffProcessorBoundary) {
  const auto params = GetParam();

  const auto result = bout::checkBoutMeshYDecomposition(
      snowflake_nype, snowflake_ny, 1, params.indices.jyseps1_1,
      params.indices.jyseps2_1, params.indices.jyseps1_2, params.indices.jyseps2_2 - 1,
      params.indices.ny_inner, MeshTopology::snowflake, params.snowflake_type);

  EXPECT_FALSE(result.success);
  EXPECT_FALSE(result.reason.empty());
}

/// The core is bounded by `jyseps2_1` in SF+ LFS and SF-, and by `jyseps1_2` in
/// SF+ HFS, so moving `jyseps2_1` has to be rejected by all of them: it also
/// breaks the W PFR middle segment for SF+ HFS.
TEST_P(SnowflakeFamilyDecompositionTest, RejectsCoreOffProcessorBoundary) {
  const auto params = GetParam();

  const auto result = bout::checkBoutMeshYDecomposition(
      snowflake_nype, snowflake_ny, 1, params.indices.jyseps1_1,
      params.indices.jyseps2_1 + 1, params.indices.jyseps1_2, params.indices.jyseps2_2,
      params.indices.ny_inner, MeshTopology::snowflake, params.snowflake_type);

  EXPECT_FALSE(result.success);
  EXPECT_FALSE(result.reason.empty());
}

/// findValidYDecomposition must only ever suggest an index set that
/// checkBoutMeshYDecomposition then accepts, for the same family member.
TEST_P(SnowflakeFamilyDecompositionTest, SuggestedDecompositionIsSelfConsistent) {
  const auto params = GetParam();

  const auto suggestion = bout::findValidYDecomposition(
      snowflake_ny, snowflake_nype, 1, 0, 1, 2, 3, 1, MeshTopology::snowflake,
      params.snowflake_type);

  ASSERT_TRUE(suggestion.success) << suggestion.reason;

  int jyseps1_1 = -1;
  int jyseps2_1 = -1;
  int jyseps1_2 = -1;
  int jyseps2_2 = -1;
  int ny_inner = -1;
  const auto found = std::sscanf(
      suggestion.reason.c_str(),
      "\t -> A valid decomposition in Y close to the one given in the grid would be: "
      "jyseps1_1=%d, jyseps2_1=%d, jyseps1_2=%d, jyseps2_2=%d, ny_inner=%d",
      &jyseps1_1, &jyseps2_1, &jyseps1_2, &jyseps2_2, &ny_inner);
  ASSERT_EQ(found, 5) << suggestion.reason;

  const auto recheck = bout::checkBoutMeshYDecomposition(
      snowflake_nype, snowflake_ny, 1, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2,
      ny_inner, MeshTopology::snowflake, params.snowflake_type);

  EXPECT_TRUE(recheck.success) << recheck.reason;
}

/// The branch-cut ordering findValidYDecomposition searches over has to match
/// the one the topology actually uses: SF+ puts both core branch cuts below
/// `ny_inner`, SF- puts `ny_inner` between them.
TEST_P(SnowflakeFamilyDecompositionTest, SuggestedDecompositionHasTheRightOrdering) {
  const auto params = GetParam();

  const auto suggestion = bout::findValidYDecomposition(
      snowflake_ny, snowflake_nype, 1, 0, 1, 2, 3, 1, MeshTopology::snowflake,
      params.snowflake_type);

  ASSERT_TRUE(suggestion.success) << suggestion.reason;

  int jyseps1_1 = -1;
  int jyseps2_1 = -1;
  int jyseps1_2 = -1;
  int jyseps2_2 = -1;
  int ny_inner = -1;
  const auto found = std::sscanf(
      suggestion.reason.c_str(),
      "\t -> A valid decomposition in Y close to the one given in the grid would be: "
      "jyseps1_1=%d, jyseps2_1=%d, jyseps1_2=%d, jyseps2_2=%d, ny_inner=%d",
      &jyseps1_1, &jyseps2_1, &jyseps1_2, &jyseps2_2, &ny_inner);
  ASSERT_EQ(found, 5) << suggestion.reason;

  EXPECT_LT(jyseps1_1, jyseps2_1) << suggestion.reason;

  const bool is_snowflake_minus =
      params.snowflake_type == SnowflakeType::SF_minus_low_field_side
      or params.snowflake_type == SnowflakeType::SF_minus_high_field_side;

  if (is_snowflake_minus) {
    EXPECT_LT(jyseps2_1, ny_inner) << suggestion.reason;
    EXPECT_LE(ny_inner, jyseps1_2) << suggestion.reason;
    EXPECT_LT(jyseps1_2, jyseps2_2) << suggestion.reason;
  } else {
    EXPECT_LT(jyseps2_1, jyseps1_2) << suggestion.reason;
    EXPECT_LT(jyseps1_2, ny_inner) << suggestion.reason;
    EXPECT_LT(ny_inner, jyseps2_2) << suggestion.reason;
  }
}

/// Regression: the double-null search must keep `ny_inner` between the two core
/// branch cuts. If it searches the snowflake ordering (`jyseps1_2 < ny_inner`)
/// it will suggest an index set that is not a double null at all.
TEST(BoutMeshDecompositionTest, ValidDoubleNullDecompositionKeepsDoubleNullOrdering) {
  const auto suggestion = bout::findValidYDecomposition(
      24, 6, 1, 0, 1, 2, 3, 1, MeshTopology::unconnected_double_null);

  ASSERT_TRUE(suggestion.success) << suggestion.reason;

  int jyseps1_1 = -1;
  int jyseps2_1 = -1;
  int jyseps1_2 = -1;
  int jyseps2_2 = -1;
  int ny_inner = -1;
  const auto found = std::sscanf(
      suggestion.reason.c_str(),
      "\t -> A valid decomposition in Y close to the one given in the grid would be: "
      "jyseps1_1=%d, jyseps2_1=%d, jyseps1_2=%d, jyseps2_2=%d, ny_inner=%d",
      &jyseps1_1, &jyseps2_1, &jyseps1_2, &jyseps2_2, &ny_inner);
  ASSERT_EQ(found, 5) << suggestion.reason;

  EXPECT_LT(jyseps1_1, jyseps2_1) << suggestion.reason;
  EXPECT_LT(jyseps2_1, ny_inner) << suggestion.reason;
  EXPECT_LE(ny_inner, jyseps1_2) << suggestion.reason;
  EXPECT_LT(jyseps1_2, jyseps2_2) << suggestion.reason;
}


  TEST(BoutMeshDecompositionTest, ValidYDecomposition) {
  int ny = 16;
  int num_y_processors = 4;
  int num_y_guards = 1;

  int jyseps1_1_start = 1;
  int jyseps2_1_start = 3;
  int jyseps1_2_start = 6;
  int jyseps2_2_start = 12;
  int ny_inner_start = 8;

  MeshTopology mesh_topology = MeshTopology::snowflake;

  auto result = bout::findValidYDecomposition(ny, num_y_processors, num_y_guards,
                                        jyseps1_1_start, jyseps2_1_start,
                                        jyseps1_2_start, jyseps2_2_start,
                                        ny_inner_start, mesh_topology);

  EXPECT_TRUE(result.success);
}

TEST(BoutMeshDecompositionTest, InvalidYDecompositionBecuaseofNYPE) {
  int ny = 8;
  int num_y_processors = 3; // deliberately incompatible
  int num_y_guards = 1;

  int jyseps1_1_start = 0;
  int jyseps2_1_start = 0;
  int jyseps1_2_start = 0;
  int jyseps2_2_start = 0;
  int ny_inner_start = 0;

  MeshTopology mesh_topology = MeshTopology::snowflake;

  auto result = bout::findValidYDecomposition(ny, num_y_processors, num_y_guards,
                                        jyseps1_1_start, jyseps2_1_start,
                                        jyseps1_2_start, jyseps2_2_start,
                                        ny_inner_start, mesh_topology);

  EXPECT_FALSE(result.success);
}

TEST(BoutMeshDecompositionTest, InvalidYDecompositionBecuaseofTopologyUDN) {
  int ny = 18;
  int num_y_processors = 9;
  int num_y_guards = 1;

  int jyseps1_1_start = 1;
  int jyseps2_1_start = 1;
  int jyseps1_2_start = 17;
  int jyseps2_2_start = 1;
  int ny_inner_start = 1;

  MeshTopology mesh_topology = MeshTopology::unconnected_double_null;

  auto result = bout::findValidYDecomposition(ny, num_y_processors, num_y_guards,
                                        jyseps1_1_start, jyseps2_1_start,
                                        jyseps1_2_start, jyseps2_2_start,
                                        ny_inner_start, mesh_topology);
  EXPECT_FALSE(result.success);
}

TEST(BoutMeshDecompositionTest, InvalidYDecompositionBecuaseofTopologySF) {
  int ny = 18;
  int num_y_processors = 9;
  int num_y_guards = 1;

  int jyseps1_1_start = 1;
  int jyseps2_1_start = 1;
  int jyseps1_2_start = 1;
  int jyseps2_2_start = 1;
  int ny_inner_start = 17;

  MeshTopology mesh_topology = MeshTopology::snowflake;

  auto result = bout::findValidYDecomposition(ny, num_y_processors, num_y_guards,
                                        jyseps1_1_start, jyseps2_1_start,
                                        jyseps1_2_start, jyseps2_2_start,
                                        ny_inner_start, mesh_topology);
  EXPECT_FALSE(result.success);
}


TEST(BoutMeshDecompositionTest, BasicValidProcessDecompositionDefaults) {
  // 8x6 grid, up to 16 processors
  auto result = bout::findValidProcessorNum(/*ny=*/8, /*nx=*/6, /*NPES=*/16);
  using ::testing::HasSubstr;
  EXPECT_TRUE(result.success);
  EXPECT_THAT(result.reason, HasSubstr("NPES=16"));
  EXPECT_THAT(result.reason, HasSubstr("NXPE=2"));
  EXPECT_THAT(result.reason, HasSubstr("NYPE=8"));
}

TEST(BoutMeshDecompositionTest, RespectsNXPE) {
  int NXPE=2;
  auto result = bout::findValidProcessorNum(/*ny=*/8, /*nx=*/8, /*NPES=*/16,
                    NXPE);
  using ::testing::HasSubstr;
  EXPECT_TRUE(result.success);
  EXPECT_THAT(result.reason, HasSubstr("NPES=16"));
  EXPECT_THAT(result.reason, HasSubstr("NXPE=2"));
  EXPECT_THAT(result.reason, HasSubstr("NYPE=8"));
}

TEST(BoutMeshDecompositionTest, RespectsNYPE) {
  int NYPE=16;
  auto result = bout::findValidProcessorNum(/*ny=*/16, /*nx=*/8, /*NPES=*/16,
                    NYPE);
  using ::testing::HasSubstr;
  EXPECT_TRUE(result.success);
  EXPECT_THAT(result.reason, HasSubstr("NPES=16"));
  EXPECT_THAT(result.reason, HasSubstr("NXPE=1"));
  EXPECT_THAT(result.reason, HasSubstr("NYPE=16"));
}


TEST(BoutMeshDecompositionTest, NoValidDecomposition) {
  // Prime sizes, limited processors
  auto result = bout::findValidProcessorNum(/*ny=*/7, /*nx=*/8, /*NPES=*/5);
  using ::testing::HasSubstr;
  EXPECT_FALSE(result.success);
  EXPECT_THAT(result.reason, HasSubstr("No valid processor decomposition found"));
}

TEST(BoutMeshDecompositionTest, SingleProcessorOnly) {
  auto result = bout::findValidProcessorNum(/*ny=*/10, /*nx=*/10, /*NPES=*/1);
  using ::testing::HasSubstr;
  EXPECT_TRUE(result.success);
  EXPECT_THAT(result.reason, HasSubstr("NPES=1"));
}

  //End of new bit

//End of the test
TEST_F(BoutMeshTest, ChooseProcessorSplitBadNXPE) {
  WithQuietOutput info{output_info};
  Options options{{"NXPE", 3}};

  BoutMeshExposer mesh(1, 24, 1, 1, 1, 8);

  EXPECT_THROW(mesh.chooseProcessorSplit(options), BoutException);
}

TEST_F(BoutMeshTest, ChooseProcessorSplitBadNYPETooManyYProcs) {
  WithQuietOutput info{output_info};
  Options options{{"NYPE", 7}};

  BoutMeshExposer mesh(1, 24, 1, 1, 1, 8);

  EXPECT_THROW(mesh.chooseProcessorSplit(options), BoutException);
}

TEST_F(BoutMeshTest, ChooseProcessorSplitBadNXPENotDivisibleByNYPE) {
  WithQuietOutput info{output_info};
  Options options{{"NXPE", 5}};

  BoutMeshExposer mesh(4, 24, 1, 1, 1, 8);

  EXPECT_THROW(mesh.chooseProcessorSplit(options), BoutException);
}

TEST_F(BoutMeshTest, ChooseProcessorSplitBadNYPENotDivisibleByNYPE) {
  WithQuietOutput info{output_info};
  Options options{{"NYPE", 5}};

  BoutMeshExposer mesh(5, 5, 1, 1, 1, 8);

  EXPECT_THROW(mesh.chooseProcessorSplit(options), BoutException);
}

TEST_F(BoutMeshTest, ChooseProcessorSplitNXPE) {
  Options options{{"NXPE", 4}};

  BoutMeshExposer mesh(4, 24, 1, 1, 1, 8);

  EXPECT_NO_THROW(mesh.chooseProcessorSplit(options));

  EXPECT_EQ(mesh.getNXPE(), 4);
  EXPECT_EQ(mesh.getNYPE(), 2);
}

TEST_F(BoutMeshTest, ChooseProcessorSplitBadNXPENotEnoughGuards) {
  Options options{{"NXPE", 4}};

  BoutMeshExposer mesh(1, 24, 1, 1, 13, 8);

  EXPECT_THROW(mesh.chooseProcessorSplit(options), BoutException);
}

TEST_F(BoutMeshTest, ChooseProcessorSplitNYPE) {
  Options options{{"NYPE", 4}};

  BoutMeshExposer mesh(1, 24, 1, 1, 1, 8);

  EXPECT_NO_THROW(mesh.chooseProcessorSplit(options));

  EXPECT_EQ(mesh.getNXPE(), 2);
  EXPECT_EQ(mesh.getNYPE(), 4);
}

TEST(getMeshTopologyTest, ReturnsCFLWhenNoXPoints) {
  BoutMeshExposer mesh(8, 8, 1, 1, 1);
  mesh.numberOfXPoints = 0;
  EXPECT_EQ(mesh.getMeshTopology(-1, 2, 3, 10, 5, 6, 7, ""), MeshTopology::closed_field_line);
}

TEST(getMeshTopologyTest, ReturnsSNWhenOneXPoint) {
  BoutMeshExposer mesh(8, 8, 1, 1, 1);
  mesh.numberOfXPoints = 1;
  EXPECT_EQ(mesh.getMeshTopology(1, 2, 2, 4, 5, 6, 7, ""), MeshTopology::single_null);
}

TEST(getMeshTopologyTest, ReturnsSFWhenSnowflakeConditionMet) {
  BoutMeshExposer mesh(8, 8, 1, 1, 1);
  mesh.numberOfXPoints = 2;
  // ny_inner between jyseps1_2 and jyseps2_2
  EXPECT_EQ(mesh.getMeshTopology(7, 39, 45, 63, 56, 8, 5, ""), MeshTopology::snowflake);
}

TEST(getMeshTopologyTest, ReturnsUDNWhenTwoXPointsDifferentIndices) {
  BoutMeshExposer mesh(8, 8, 1, 1, 1);
  mesh.numberOfXPoints = 2;
  // ny_inner not between jyseps1_2 and jyseps2_2
  EXPECT_EQ(mesh.getMeshTopology(0, 0, 10, 20, 25, 1, 2, ""), MeshTopology::unconnected_double_null);
}

TEST(getMeshTopologyTest, ReturnsCDNWhenTwoXPointsSameIndices) {
  BoutMeshExposer mesh(8, 8, 1, 1, 1);
  mesh.numberOfXPoints = 2;
  // ny_inner not between jyseps1_2 and jyseps2_2 but ixseps1 == ixseps2
  EXPECT_EQ(mesh.getMeshTopology(0, 0, 10, 20, 25, 1, 1, ""), MeshTopology::connected_double_null);
}

struct SnowflakeTypeParameters {
  MeshTopology mesh_topology;
  std::string ingrid_topology;
  SnowflakeType expected;
  std::string test_name;
};

std::ostream& operator<<(std::ostream& out, const SnowflakeTypeParameters& value) {
  return out << "SnowflakeTypeParameters{mesh_topology=" << toString(value.mesh_topology)
             << ", ingrid_topology='" << value.ingrid_topology
             << "', expected=" << toString(value.expected) << "}";
}

std::string SnowflakeTypeParametersToString(
    const ::testing::TestParamInfo<SnowflakeTypeParameters>& param) {
  return param.param.test_name;
}

struct GetSnowflakeTypeTest : public ::testing::TestWithParam<SnowflakeTypeParameters> {};

INSTANTIATE_TEST_SUITE_P(
    SnowflakeFamily, GetSnowflakeTypeTest,
    ::testing::Values(
        SnowflakeTypeParameters{MeshTopology::snowflake, "SF15",
                                SnowflakeType::SF_minus_low_field_side, "SF15"},
        SnowflakeTypeParameters{MeshTopology::snowflake, "SF45",
                                SnowflakeType::SF_plus_low_field_side, "SF45"},
        SnowflakeTypeParameters{MeshTopology::snowflake, "SF75",
                                SnowflakeType::SF_plus_low_field_side, "SF75"},
        SnowflakeTypeParameters{MeshTopology::snowflake, "SF105",
                                SnowflakeType::SF_plus_high_field_side, "SF105"},
        SnowflakeTypeParameters{MeshTopology::snowflake, "SF135",
                                SnowflakeType::SF_plus_high_field_side, "SF135"},
        SnowflakeTypeParameters{MeshTopology::snowflake, "SF165",
                                SnowflakeType::SF_minus_high_field_side, "SF165"},
        // The INGRID label is upper-cased by readIngridTopology, but
        // getSnowflakeType is also called directly with raw strings
        SnowflakeTypeParameters{MeshTopology::snowflake, "sf165",
                                SnowflakeType::SF_minus_high_field_side, "SF165LowerCase"},
        SnowflakeTypeParameters{MeshTopology::snowflake, " SF45 ",
                                SnowflakeType::SF_plus_low_field_side, "SF45Padded"}),
    SnowflakeTypeParametersToString);

INSTANTIATE_TEST_SUITE_P(
    XPointTarget, GetSnowflakeTypeTest,
    ::testing::Values(
        SnowflakeTypeParameters{MeshTopology::XPoint_target, "XPOINT_TARGET",
                                SnowflakeType::XPT, "NormalisedFromGrid"},
        SnowflakeTypeParameters{MeshTopology::XPoint_target, "XPoint_target",
                                SnowflakeType::XPT, "AsWrittenInGrid"},
        SnowflakeTypeParameters{MeshTopology::XPoint_target, "xpoint_target",
                                SnowflakeType::XPT, "LowerCase"}),
    SnowflakeTypeParametersToString);

INSTANTIATE_TEST_SUITE_P(
    GenericSnowflake, GetSnowflakeTypeTest,
    ::testing::Values(
        SnowflakeTypeParameters{MeshTopology::snowflake, "SF", SnowflakeType::SF, "PlainSF"},
        SnowflakeTypeParameters{MeshTopology::snowflake, "", SnowflakeType::SF, "EmptyString"},
        SnowflakeTypeParameters{MeshTopology::snowflake, "SF-IDEAL", SnowflakeType::SF,
                                "SFideal"},
        SnowflakeTypeParameters{MeshTopology::snowflake, "SF999", SnowflakeType::SF,
                                "UnknownFamilyAngle"}),
    SnowflakeTypeParametersToString);

INSTANTIATE_TEST_SUITE_P(
    NotASnowflake, GetSnowflakeTypeTest,
    ::testing::Values(
        SnowflakeTypeParameters{MeshTopology::closed_field_line, "closed_field_line", SnowflakeType::SF,
                                "ClosedFieldLine"},
        SnowflakeTypeParameters{MeshTopology::single_null, "single_null", SnowflakeType::SF, "SingleNull"},
        SnowflakeTypeParameters{MeshTopology::unconnected_double_null, "unconnected_double_null", SnowflakeType::SF,
                                "UnconnectedDoubleNull"},
        SnowflakeTypeParameters{MeshTopology::connected_double_null, "connected_double_null", SnowflakeType::SF,
                                "ConnectedDoubleNull"}),
    SnowflakeTypeParametersToString);

TEST_P(GetSnowflakeTypeTest, ClassifiesSnowflakeFamilyMember) {
  const auto params = GetParam();

  BoutMeshExposer mesh(8, 8, 1, 1, 1);

  EXPECT_EQ(mesh.getSnowflakeType(params.mesh_topology, params.ingrid_topology),
            params.expected);
}

// getMeshTopology driven by the INGRID `topology` label rather than the
// separatrix indices.

struct MeshTopologyFromLabelParameters {
  std::string ingrid_topology;
  MeshTopology expected;
  std::string test_name;
};

std::ostream& operator<<(std::ostream& out, const MeshTopologyFromLabelParameters& value) {
  return out << "MeshTopologyFromLabelParameters{ingrid_topology='"
             << value.ingrid_topology << "', expected=" << toString(value.expected) << "}";
}

std::string MeshTopologyFromLabelParametersToString(
    const ::testing::TestParamInfo<MeshTopologyFromLabelParameters>& param) {
  return param.param.test_name;
}

struct GetMeshTopologyFromLabelTest
    : public ::testing::TestWithParam<MeshTopologyFromLabelParameters> {};

INSTANTIATE_TEST_SUITE_P(
    AcronymLabels, GetMeshTopologyFromLabelTest,
    ::testing::Values(
        MeshTopologyFromLabelParameters{"CFL", MeshTopology::closed_field_line, "CFL"},
        MeshTopologyFromLabelParameters{"SN", MeshTopology::single_null, "SN"},
        MeshTopologyFromLabelParameters{"UDN", MeshTopology::unconnected_double_null,
                                        "UDN"},
        MeshTopologyFromLabelParameters{"CDN", MeshTopology::connected_double_null, "CDN"},
        MeshTopologyFromLabelParameters{"SF45", MeshTopology::snowflake, "SF45"},
        MeshTopologyFromLabelParameters{"SF165", MeshTopology::snowflake, "SF165"},
        MeshTopologyFromLabelParameters{"XPT", MeshTopology::XPoint_target, "XPT"}),
    MeshTopologyFromLabelParametersToString);

INSTANTIATE_TEST_SUITE_P(
    LongLabels, GetMeshTopologyFromLabelTest,
    ::testing::Values(
        MeshTopologyFromLabelParameters{"CLOSED_FIELD_LINE",
                                        MeshTopology::closed_field_line,
                                        "ClosedFieldLine"},
        MeshTopologyFromLabelParameters{"SINGLE_NULL", MeshTopology::single_null,
                                        "SingleNull"},
        MeshTopologyFromLabelParameters{"UNCONNECTED_DOUBLE_NULL",
                                        MeshTopology::unconnected_double_null,
                                        "UnconnectedDoubleNull"},
        MeshTopologyFromLabelParameters{"CONNECTED_DOUBLE_NULL",
                                        MeshTopology::connected_double_null,
                                        "ConnectedDoubleNull"},
        MeshTopologyFromLabelParameters{"SNOWFLAKE", MeshTopology::snowflake, "Snowflake"},
        MeshTopologyFromLabelParameters{"XPOINT_TARGET", MeshTopology::XPoint_target,
                                        "XPointTarget"}),
    MeshTopologyFromLabelParametersToString);

/// The label wins over the separatrix indices. The indices used here would be
/// classified as a snowflake by the index-based fallback (ny_inner lies between
/// jyseps1_2 and jyseps2_2), so any answer other than the label's is the
/// fallback leaking through.
TEST_P(GetMeshTopologyFromLabelTest, LabelTakesPrecedenceOverIndices) {
  WithQuietOutput warn{output_warn};
  const auto params = GetParam();

  BoutMeshExposer mesh(8, 8, 1, 1, 1);

  EXPECT_EQ(mesh.getMeshTopology(7, 39, 45, 63, 56, 8, 5, params.ingrid_topology),
            params.expected);
}

TEST_F(BoutMeshTest, GetMeshTopologyUnrecognisedLabelFallsBackToIndices) {
  BoutMeshExposer mesh(8, 8, 1, 1, 1);

  // Not a label the reader knows about, so the separatrix indices decide. These
  // indices are a snowflake.
  EXPECT_EQ(mesh.getMeshTopology(7, 39, 45, 63, 56, 8, 5, "BN"),
            MeshTopology::snowflake);
}

TEST(GetMeshTopologyTest, UnrecognisedLabelWarns) {
  WithQuietOutput info{output_info};
  WithQuietOutput debug{output_debug};
  WithQuietOutput progress{output_progress};

  std::stringstream buffer;
  auto* old_buffer = std::cout.rdbuf(buffer.rdbuf());
  {
    BoutMeshExposer mesh(8, 8, 1, 1, 1);
    mesh.getMeshTopology(7, 39, 45, 63, 56, 8, 5, "BN");
  }
  std::cout.rdbuf(old_buffer);

  EXPECT_THAT(buffer.str(), ::testing::HasSubstr("BN"));
}

// readIngridTopology
/// deliberately *not* using `WithQuietOutput` on `output_warn`, because the
struct ReadIngridTopologyTest : public ::testing::Test {
  ReadIngridTopologyTest() : old_cout_buffer(std::cout.rdbuf()) {
    std::cout.rdbuf(buffer.rdbuf());
  }
  ~ReadIngridTopologyTest() override { std::cout.rdbuf(old_cout_buffer); }
  static constexpr auto missing_topology_warning = "no 'topology' variable";

  /// Make a grid source containing `topology = some_value`
  static GridDataSource* gridWithTopology(const std::string& value) {
    Options values;
    values["topology"] = value;
    return new FakeGridDataSource{values};
  }

  std::stringstream buffer;
  std::streambuf* old_cout_buffer;
};

struct IngridTopologyParameters {
  std::string grid_value; // Value of `topology` in the grid file
  std::string expected;   // Expected value of `BoutMesh::IngridTopology`
  std::string test_name;
};

std::ostream& operator<<(std::ostream& out, const IngridTopologyParameters& value) {
  return out << "IngridTopologyParameters{grid_value='" << value.grid_value
             << "', expected='" << value.expected << "'}";
}

std::string IngridTopologyParametersToString(
    const ::testing::TestParamInfo<IngridTopologyParameters>& param) {
  return param.param.test_name;
}

struct ReadIngridTopologyParameterisedTest
    : public ReadIngridTopologyTest,
      public ::testing::WithParamInterface<IngridTopologyParameters> {};

INSTANTIATE_TEST_SUITE_P(
    KnownTopologies, ReadIngridTopologyParameterisedTest,
    ::testing::Values(
        IngridTopologyParameters{"LSN", "LSN", "LowerSingleNull"},
        IngridTopologyParameters{"UDN", "UDN", "UnconnectedDoubleNull"},
        IngridTopologyParameters{"SF75", "SF75", "SF75"},
        IngridTopologyParameters{"SF165", "SF165", "SF165"},
        IngridTopologyParameters{"SF+", "SF+", "SFplus"},
        IngridTopologyParameters{"SF-", "SF-", "SFminus"}),
    IngridTopologyParametersToString);


INSTANTIATE_TEST_SUITE_P(
    NormalizedTopologies, ReadIngridTopologyParameterisedTest,
    ::testing::Values(
        IngridTopologyParameters{"sf75", "SF75", "LowerCase"},
        IngridTopologyParameters{"Sf165", "SF165", "MixedCase"},
        IngridTopologyParameters{"  SF+  ", "SF+", "SurroundingSpaces"},
        IngridTopologyParameters{"\tsf-\n", "SF-", "SurroundingWhitespaceAndLowerCase"},
        IngridTopologyParameters{" udn ", "UDN", "UDNWithSpaces"}),
    IngridTopologyParametersToString);

TEST_P(ReadIngridTopologyParameterisedTest, ReadsTopologyFromGrid) {
  const auto params = GetParam();

  BoutMeshExposer mesh{gridWithTopology(params.grid_value), nullptr};

  EXPECT_EQ(mesh.readIngridTopology(), params.expected);
  EXPECT_EQ(mesh.IngridTopology, params.expected);

  // A grid that has no topology should not warn about a missing one
  EXPECT_THAT(buffer.str(),
              ::testing::Not(::testing::HasSubstr(missing_topology_warning)));
}

TEST_F(ReadIngridTopologyTest, NoTopologyInGridWarnsAndLeavesEmpty) {
  BoutMeshExposer mesh{new FakeGridDataSource, nullptr};

  EXPECT_EQ(mesh.readIngridTopology(), "");
  EXPECT_EQ(mesh.IngridTopology, "");
  EXPECT_THAT(buffer.str(), ::testing::HasSubstr(missing_topology_warning));
}

TEST_F(ReadIngridTopologyTest, EmptyTopologyStringWarnsAndLeavesEmpty) {
  // `topology` is present but empty
  BoutMeshExposer mesh{gridWithTopology(""), nullptr};

  EXPECT_EQ(mesh.readIngridTopology(), "");
  EXPECT_EQ(mesh.IngridTopology, "");
  EXPECT_THAT(buffer.str(), ::testing::HasSubstr(missing_topology_warning));
}

TEST_F(ReadIngridTopologyTest, WhitespaceOnlyTopologyWarnsAndLeavesEmpty) {
  // Trimming a whitespace-only string leaves nothing, so treat it as absent
  BoutMeshExposer mesh{gridWithTopology("  \t\n "), nullptr};

  EXPECT_EQ(mesh.readIngridTopology(), "");
  EXPECT_EQ(mesh.IngridTopology, "");
  EXPECT_THAT(buffer.str(), ::testing::HasSubstr(missing_topology_warning));
}

TEST_F(ReadIngridTopologyTest, NoTopologyIsNotFatal) {
  BoutMeshExposer mesh{new FakeGridDataSource, nullptr};

  EXPECT_NO_THROW(mesh.readIngridTopology());
}

TEST_F(ReadIngridTopologyTest, UnrecognisedTopologyIsPassedThroughUnchanged) {
  // A new INGRID configuration doesn't break loading the grid
  BoutMeshExposer mesh{gridWithTopology("SF-ideal"), nullptr};

  EXPECT_EQ(mesh.readIngridTopology(), "SF-IDEAL");
  EXPECT_THAT(buffer.str(),
              ::testing::Not(::testing::HasSubstr(missing_topology_warning)));
}

TEST_F(ReadIngridTopologyTest, RereadingOverwritesPreviousValue) {
  BoutMeshExposer mesh{gridWithTopology("SF165"), nullptr};

  EXPECT_EQ(mesh.readIngridTopology(), "SF165");
  EXPECT_EQ(mesh.readIngridTopology(), "SF165");
  EXPECT_EQ(mesh.IngridTopology, "SF165");
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

TEST_F(BoutMeshTest, YProc) {
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

TEST_F(BoutMeshTest, XProc) {
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

TEST_F(BoutMeshTest, GetGlobalXIndex) {
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

TEST_F(BoutMeshTest, GetGlobalXIndexNoBoundaries) {
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

TEST_F(BoutMeshTest, GlobalXIntSymmetricX) {
  BoutMeshExposer mesh01(4, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.GlobalX(0), -0.125);
  EXPECT_EQ(mesh01.GlobalX(1), 0.125);
  EXPECT_EQ(mesh01.GlobalX(2), 0.375);
  EXPECT_EQ(mesh01.GlobalX(3), 0.625);
  EXPECT_EQ(mesh01.GlobalX(4), 0.875);
}

TEST_F(BoutMeshTest, GlobalXIntAsymmetricX) {
  BoutMeshExposer mesh01(4, 3, 1, 2, 2, 0, 1, false, false);
  EXPECT_EQ(mesh01.GlobalX(0), 0.);
  EXPECT_EQ(mesh01.GlobalX(1), 0.25);
  EXPECT_EQ(mesh01.GlobalX(2), 0.5);
  EXPECT_EQ(mesh01.GlobalX(3), 0.75);
  EXPECT_EQ(mesh01.GlobalX(4), 1.0);
}

TEST_F(BoutMeshTest, GlobalXRealSymmetricX) {
  BoutMeshExposer mesh01(4, 3, 1, 2, 2, 0, 1);
  EXPECT_EQ(mesh01.GlobalX(0.5), 0.);
  EXPECT_EQ(mesh01.GlobalX(1.5), 0.25);
  EXPECT_EQ(mesh01.GlobalX(2.5), 0.5);
  EXPECT_EQ(mesh01.GlobalX(3.5), 0.75);
  EXPECT_EQ(mesh01.GlobalX(4.5), 1.0);
}

TEST_F(BoutMeshTest, GlobalXRealAsymmetricX) {
  BoutMeshExposer mesh01(4, 3, 1, 2, 2, 0, 1, false, false);
  EXPECT_EQ(mesh01.GlobalX(0.5), 0.125);
  EXPECT_EQ(mesh01.GlobalX(1.5), 0.375);
  EXPECT_EQ(mesh01.GlobalX(2.5), 0.625);
  EXPECT_EQ(mesh01.GlobalX(3.5), 0.875);
  EXPECT_EQ(mesh01.GlobalX(4.5), 1.125);
}

TEST_F(BoutMeshTest, GetLocalXIndex) {
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

TEST_F(BoutMeshTest, GetLocalXIndexNoBoundaries) {
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

TEST_F(BoutMeshTest, GetGlobalYIndexSingleNull) {
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

TEST_F(BoutMeshTest, GetGlobalYIndexDoubleNull) {
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

TEST_F(BoutMeshTest, GetGlobalYIndexNoBoundaries) {
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

TEST_F(BoutMeshTest, GetLocalYIndexSingleNull) {
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

TEST_F(BoutMeshTest, GetLocalYIndexDoubleNull) {
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

TEST_F(BoutMeshTest, GetLocalYIndexNoBoundaries) {
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

TEST_F(BoutMeshTest, GlobalYIntSymmetricY) {
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

TEST_F(BoutMeshTest, GlobalYIntAsymmetricY) {
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

TEST_F(BoutMeshTest, GlobalYRealSymmetricY) {
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

TEST_F(BoutMeshTest, GlobalYRealAsymmetricY) {
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

TEST_F(BoutMeshTest, GetGlobalZIndex) {
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

TEST_F(BoutMeshTest, GetGlobalZIndexNoBoundaries) {
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

TEST_F(BoutMeshTest, GetLocalZIndex) {
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

TEST_F(BoutMeshTest, GetLocalZIndexNoBoundaries) {
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

TEST_F(BoutMeshTest, FirstX) {
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

TEST_F(BoutMeshTest, LastX) {
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

TEST_F(BoutMeshTest, FirstY) {
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

TEST_F(BoutMeshTest, LastY) {
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

TEST_F(BoutMeshTest, DefaultConnectionsCore1x1) {
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

TEST_F(BoutMeshTest, TopologySOL2x2) {
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

TEST_F(BoutMeshTest, TopologySOLPeriodicX2x2) {
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

TEST_F(BoutMeshTest, TopologySingleNull2x3) {
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

TEST_F(BoutMeshTest, TopologyDisconnectedDoubleNull1x6) {
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

TEST_F(BoutMeshTest, SetDerivedGridSizes) {
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

TEST_F(BoutMeshTest, CreateXBoundariesPeriodicX) {
  // Periodic in X, so no boundaries
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 3, 6, 1, 0}));
  mesh.periodicX = true;
  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_TRUE(boundaries.empty());
}

TEST_F(BoutMeshTest, CreateXBoundariesNoGuards) {
  // No guards in X, so no boundaries
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 0, 1, 3, 6, 1, 0}));
  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_TRUE(boundaries.empty());
}

TEST_F(BoutMeshTest, CreateXBoundariesDoubleNullInsidePF) {
  // Three cores in X, inside core, one boundary
  BoutMeshExposer mesh_inside(createDisconnectedDoubleNull({12, 3, 1, 1, 3, 6, 0, 0}));
  mesh_inside.createXBoundaries();

  auto boundaries_inside = mesh_inside.getBoundaries();
  EXPECT_EQ(boundaries_inside.size(), 1);
  EXPECT_EQ(boundaries_inside[0]->label, "pf");
}

TEST_F(BoutMeshTest, CreateXBoundariesDoubleNullMiddlePF) {
  // Three cores in X, middle core, so no boundaries
  BoutMeshExposer mesh_middle(createDisconnectedDoubleNull({12, 3, 1, 1, 3, 6, 1, 0}));
  mesh_middle.createXBoundaries();

  auto boundaries_middle = mesh_middle.getBoundaries();
  EXPECT_TRUE(boundaries_middle.empty());
}

TEST_F(BoutMeshTest, CreateXBoundariesDoubleNullOutsidePF) {
  // Three cores in X, outside core, one boundary
  BoutMeshExposer mesh_inside(createDisconnectedDoubleNull({12, 3, 1, 1, 3, 6, 0, 0}));
  mesh_inside.createXBoundaries();

  auto boundaries_inside = mesh_inside.getBoundaries();
  EXPECT_EQ(boundaries_inside.size(), 1);
  EXPECT_EQ(boundaries_inside[0]->label, "pf");
}

TEST_F(BoutMeshTest, CreateXBoundariesDoubleNullInsideOutsideCore) {
  // One core in X, so we expect two boundaries
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 2);
  EXPECT_EQ(boundaries[0]->label, "core");
  EXPECT_EQ(boundaries[1]->label, "sol");
}

//New snowflake tests for X boundaries
TEST_F(BoutMeshTest, CreateXBoundariesSnowflakeOuterSOL) {
  WithQuietOutput info{output_info};
  // Snowflake topology: outer X boundary below ny_inner is SOL
  // PE_XIND = NXPE - 1, PE_YIND chosen so yg <= ny_inner
  BoutMeshExposer mesh(createSnowflake(
  {/* nx */ 12,
    /* ny */ 24,
    /* MXG */ 1,
    /* MYG */ 1,
    /* nxpe */ 2,
    /* nype */ 6,
    /* pe_xind */ 1,   // outer X
    /* pe_yind */ 1    // Y chosen so yg <= ny_inner
  }));

  // Ensure PE_YIND corresponds to yg <= ny_inner
  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  ASSERT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "sol");
}


TEST_F(BoutMeshTest, CreateXBoundariesSnowflakeSouthPFOuter) {
  WithQuietOutput info{output_info};
  // Snowflake topology: outer X boundary above ny_inner is PF outer region
  // PE_XIND = NXPE - 1, PE_YIND chosen so yg > ny_inner
  // ny >= nype is needed
  BoutMeshExposer mesh(createSnowflake(
      {12, 24, 1, 1,
       3, 6,
       2, 5}));   // PE_YIND high enough → yg > ny_inner

  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  ASSERT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "south_pf_outer");
}

TEST_F(BoutMeshTest, CreateXBoundariesSnowflakeInnerPF) {
  WithQuietOutput info{output_info};
  // Snowflake topology: inner X boundary inside PF region
  // PE_XIND = 0
  BoutMeshExposer mesh(createSnowflake(
      {12, 3, 1, 1,
       3, 6,
       0, 0}));

  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  ASSERT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "pf");
}

TEST_F(BoutMeshTest, CreateXBoundariesSnowflakeSingleCoreInX) {
  WithQuietOutput info{output_info};
  // Snowflake topology: one core in X → two boundaries
  BoutMeshExposer mesh(createSnowflake(
      {
        /* nx */ 12,
        /* ny */ 24,
        /* MXG */ 1,
        /* MYG */ 1,
        /* nxpe */ 1,
        /* nype */ 6,
        /* pe_xind */ 0,
        /* pe_yind */ 2
      }));

  mesh.createXBoundaries();

  auto boundaries = mesh.getBoundaries();
  ASSERT_EQ(boundaries.size(), 2);

  // Inner must be either core or PF
  EXPECT_TRUE(boundaries[0]->label == "core"
          || boundaries[0]->label == "pf");

  // Outer must be SOL or south PF outer
  EXPECT_TRUE(boundaries[1]->label == "sol"
          || boundaries[1]->label == "south_pf_outer");
}

////////////////////////////////////////////////////////////
// Mesh-level tests for each member of the snowflake family
//
// These build the mesh with `create_topology = false` so that
// `snowflake_type` can be set before `topology()` runs: the family member is
// not recoverable from the separatrix indices alone, it comes from the INGRID
// label.

bool isSnowflakeMinus(SnowflakeType snowflake_type) {
  return snowflake_type == SnowflakeType::SF_minus_low_field_side
         or snowflake_type == SnowflakeType::SF_minus_high_field_side;
}

/// Set up `mesh` (built with `nype` = 6 and MYSUB = `snowflake_mysub`) as the
/// given member of the snowflake family and build its topology.
///
/// SF+ has its second X-point in the private flux region, so its separatrix is
/// at smaller x than the primary (`ixseps2 < ixseps1`); SF- has it in the SOL,
/// so the other way round.
void buildSnowflakeTopology(BoutMeshExposer& mesh, SnowflakeType snowflake_type) {
  const bool is_minus = isSnowflakeMinus(snowflake_type);

  mesh.setXDecompositionIndices(is_minus ? BoutMeshExposer::XDecompositionIndices{2, 4}
                                         : BoutMeshExposer::XDecompositionIndices{4, 2});
  mesh.setYDecompositionIndices(is_minus ? snowflakeMinusIndices(snowflake_mysub)
                                         : snowflakePlusIndices(snowflake_mysub));
  mesh.mesh_topology = MeshTopology::snowflake;
  mesh.snowflake_type = snowflake_type;
  mesh.topology();
}

/// Global y index of the first cell of the core, for this family member.
int firstCoreCell(SnowflakeType snowflake_type) {
  const auto indices = isSnowflakeMinus(snowflake_type)
                           ? snowflakeMinusIndices(snowflake_mysub)
                           : snowflakePlusIndices(snowflake_mysub);
  // SF+ HFS is the one member whose core starts at jyseps2_1 instead of
  // jyseps1_1
  return (snowflake_type == SnowflakeType::SF_plus_high_field_side ? indices.jyseps2_1
                                                                   : indices.jyseps1_1)
         + 1;
}

struct SnowflakeFamilyMeshTest : public ::testing::TestWithParam<SnowflakeType> {
  WithQuietOutput debug{output_debug};
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  WithQuietOutput progress{output_progress};
};

std::string SnowflakeTypeToTestName(const ::testing::TestParamInfo<SnowflakeType>& param) {
  return toString(param.param);
}

INSTANTIATE_TEST_SUITE_P(SnowflakeFamily, SnowflakeFamilyMeshTest,
                         ::testing::Values(SnowflakeType::SF_plus_low_field_side,
                                           SnowflakeType::SF_plus_high_field_side,
                                           SnowflakeType::SF_minus_low_field_side,
                                           SnowflakeType::SF_minus_high_field_side,
                                           SnowflakeType::SF),
                         SnowflakeTypeToTestName);

/// A layout where every region is exactly one processor long has to be
/// accepted by `topology()` for every member of the family.
TEST_P(SnowflakeFamilyMeshTest, TopologyAcceptsAlignedBranchCuts) {
  BoutMeshExposer mesh(6, snowflake_mysub, 1, 1, snowflake_nype, 0, 0, false);

  EXPECT_NO_THROW(buildSnowflakeTopology(mesh, GetParam()));
}

/// The inner X face is always either core or private flux, never nothing.
TEST_P(SnowflakeFamilyMeshTest, CreateXBoundariesAlwaysLabelsTheInnerFace) {
  // nxpe = 2 and pe_xind = 0, so only the inner X face is a boundary
  BoutMeshExposer mesh(6, snowflake_mysub, 1, 2, snowflake_nype, 0, 2, false);
  buildSnowflakeTopology(mesh, GetParam());

  mesh.createXBoundaries();

  const auto boundaries = mesh.getBoundaries();
  ASSERT_EQ(boundaries.size(), 1);
  EXPECT_TRUE(boundaries[0]->label == "core" or boundaries[0]->label == "pf")
      << boundaries[0]->label;
}

/// The processor holding the first core cell must get a "core" inner boundary.
TEST_P(SnowflakeFamilyMeshTest, CreateXBoundariesLabelsTheCoreAsCore) {
  const int core_yproc = firstCoreCell(GetParam()) / snowflake_mysub;

  BoutMeshExposer mesh(6, snowflake_mysub, 1, 2, snowflake_nype, 0, core_yproc, false);
  buildSnowflakeTopology(mesh, GetParam());

  mesh.createXBoundaries();

  const auto boundaries = mesh.getBoundaries();
  ASSERT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "core");
}

/// The processor holding the west target is private flux, not core.
TEST_P(SnowflakeFamilyMeshTest, CreateXBoundariesLabelsTheWestTargetAsPF) {
  BoutMeshExposer mesh(6, snowflake_mysub, 1, 2, snowflake_nype, 0, 0, false);
  buildSnowflakeTopology(mesh, GetParam());

  mesh.createXBoundaries();

  const auto boundaries = mesh.getBoundaries();
  ASSERT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "pf");
}

/// `ySize` returns the number of y points along the field line through a given
/// x position, so it has to be positive and can never exceed `ny`. This is the
/// invariant that catches an x/y pair falling through every branch.
TEST_P(SnowflakeFamilyMeshTest, YSizeIsAlwaysAValidFieldLineLength) {
  for (int pe_yind = 0; pe_yind < snowflake_nype; ++pe_yind) {
    SCOPED_TRACE(fmt::format("pe_yind = {}", pe_yind));

    BoutMeshExposer mesh(6, snowflake_mysub, 1, 1, snowflake_nype, 0, pe_yind, false);
    buildSnowflakeTopology(mesh, GetParam());

    for (int xpos = 0; xpos < mesh.LocalNx; ++xpos) {
      SCOPED_TRACE(fmt::format("xpos = {}", xpos));
      const int y_size = mesh.ySize(xpos);
      EXPECT_GT(y_size, 0);
      EXPECT_LE(y_size, snowflake_ny);
    }
  }
}

/// Every member of the family has a single outer SOL, running from y = 0 up to
/// the east target at y = ny_inner - 1.
TEST_P(SnowflakeFamilyMeshTest, YSizeOfTheOuterSOL) {
  const auto indices = isSnowflakeMinus(GetParam())
                           ? snowflakeMinusIndices(snowflake_mysub)
                           : snowflakePlusIndices(snowflake_mysub);

  // pe_yind = 0 is below ny_inner for both layouts
  BoutMeshExposer mesh(6, snowflake_mysub, 1, 1, snowflake_nype, 0, 0, false);
  buildSnowflakeTopology(mesh, GetParam());

  // Outermost x point, well outside both separatrices
  EXPECT_EQ(mesh.ySize(mesh.LocalNx - 1), indices.ny_inner);
}

/// `GlobalY` normalises y over the core, so the first core cell sits half a
/// cell in: 0.5 / (number of core cells).
TEST_P(SnowflakeFamilyMeshTest, GlobalYIsNormalisedOverThisMembersCore) {
  const auto snowflake_type = GetParam();
  const auto indices = isSnowflakeMinus(snowflake_type)
                           ? snowflakeMinusIndices(snowflake_mysub)
                           : snowflakePlusIndices(snowflake_mysub);

  const int first_core = firstCoreCell(snowflake_type);
  const int core_length = snowflake_type == SnowflakeType::SF_plus_high_field_side
                              ? indices.jyseps1_2 - indices.jyseps2_1
                              : indices.jyseps2_1 - indices.jyseps1_1;

  BoutMeshExposer mesh(6, snowflake_mysub, 1, 1, snowflake_nype, 0,
                       first_core / snowflake_mysub, false);
  buildSnowflakeTopology(mesh, snowflake_type);

  // Local y = MYG is the first non-guard cell, which is the first core cell on
  // this processor
  ASSERT_EQ(mesh.getGlobalYIndexNoBoundaries(mesh.ystart), first_core);

  EXPECT_DOUBLE_EQ(mesh.GlobalY(mesh.ystart), 0.5 / core_length);
}

/// Regression: a closed field line has no X-points, so `topology()` used to
/// fall into the single-null branch via `jyseps2_1 == jyseps1_2`. It still has
/// to set up the same connections now that the branch is chosen by
/// `mesh_topology`.
TEST_F(BoutMeshTest, TopologyClosedFieldLineMatchesSingleNull) {
  const BoutMeshExposer::YDecompositionIndices core_indices{-1, 1, 1, 2, 1};
  const BoutMeshExposer::XDecompositionIndices no_separatrix{5, 5};

  BoutMeshExposer closed_field_line(5, 3, 1, 1, 1, 0, 0, false);
  closed_field_line.setXDecompositionIndices(no_separatrix);
  closed_field_line.setYDecompositionIndices(core_indices);
  closed_field_line.mesh_topology = MeshTopology::closed_field_line;
  closed_field_line.topology();

  BoutMeshExposer single_null(5, 3, 1, 1, 1, 0, 0, false);
  single_null.setXDecompositionIndices(no_separatrix);
  single_null.setYDecompositionIndices(core_indices);
  single_null.mesh_topology = MeshTopology::single_null;
  single_null.topology();

  EXPECT_EQ(closed_field_line.getConnectionInfo(), single_null.getConnectionInfo());
}


TEST_F(BoutMeshTest, CreateYBoundariesNoGuards) {
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 0, 1, 6, 0, 0}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_TRUE(boundaries.empty());
}

TEST_F(BoutMeshTest, CreateYBoundariesClosedFieldLines) {
  BoutMeshExposer mesh(createCore({4, 4, 2, 2, 4, 4}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_TRUE(boundaries.empty());
}

TEST_F(BoutMeshTest, CreateYBoundariesInnerLower) {
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 0}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "lower_target");
}

TEST_F(BoutMeshTest, CreateYBoundariesInnerUpper) {
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 2}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "upper_target");
}

TEST_F(BoutMeshTest, CreateYBoundariesOuterUpper) {
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 5}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "upper_target");
}

TEST_F(BoutMeshTest, CreateYBoundariesOuterLower) {
  BoutMeshExposer mesh(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 3}));
  mesh.createYBoundaries();

  auto boundaries = mesh.getBoundaries();
  EXPECT_EQ(boundaries.size(), 1);
  EXPECT_EQ(boundaries[0]->label, "lower_target");
}

TEST_F(BoutMeshTest, PeriodicY) {
  BoutMeshExposer mesh00(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 0}));
  EXPECT_FALSE(mesh00.periodicY(2));
  EXPECT_FALSE(mesh00.periodicY(10));

  BoutMeshExposer mesh01(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
  EXPECT_TRUE(mesh01.periodicY(2));
  EXPECT_FALSE(mesh01.periodicY(10));
}

TEST_F(BoutMeshTest, PeriodicYWithShiftAngle) {
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

TEST_F(BoutMeshTest, NumberOfYBoundaries) {
  BoutMeshExposer mesh_SOL(createSOL({3, 3, 1, 1, 2, 2, 1, 1}));
  EXPECT_EQ(mesh_SOL.numberOfYBoundaries(), 1);

  BoutMeshExposer mesh_DND(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 0}));
  EXPECT_EQ(mesh_DND.numberOfYBoundaries(), 2);
}

TEST_F(BoutMeshTest, HasBranchCutLower) {
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

TEST_F(BoutMeshTest, HasBranchCutUpper) {
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

TEST_F(BoutMeshTest, GetPossibleBoundariesCore) {
  BoutMeshExposer mesh_core_1x1(createCore({12, 3, 1, 1, 1, 1, 0, 0}));
  BoutMeshExposer mesh_core_32x64(createCore({12, 3, 1, 1, 32, 64, 7, 4}));

  std::set<std::string> boundaries{"core", "sol"};

  EXPECT_EQ(mesh_core_1x1.getPossibleBoundaries(), boundaries);
  EXPECT_EQ(mesh_core_32x64.getPossibleBoundaries(), boundaries);
}

TEST_F(BoutMeshTest, GetPossibleBoundariesCorePeriodicX) {
  BoutMeshExposer mesh_core_1x1(createCore({12, 3, 1, 1, 1, 1, 0, 0}), true);
  BoutMeshExposer mesh_core_32x64(createCore({12, 3, 1, 1, 32, 64, 7, 4}), true);

  EXPECT_TRUE(mesh_core_1x1.getPossibleBoundaries().empty());
  EXPECT_TRUE(mesh_core_32x64.getPossibleBoundaries().empty());
}

TEST_F(BoutMeshTest, GetPossibleBoundariesDND) {
  BoutMeshExposer mesh_DND_1x6(createDisconnectedDoubleNull({12, 3, 1, 1, 1, 6, 0, 1}));
  BoutMeshExposer mesh_DND_32x64(
      createDisconnectedDoubleNull({12, 3, 1, 1, 32, 64, 0, 4}));

  std::set<std::string> boundaries{"core", "pf", "sol", "upper_target", "lower_target"};

  EXPECT_EQ(mesh_DND_1x6.getPossibleBoundaries(), boundaries);
  EXPECT_EQ(mesh_DND_32x64.getPossibleBoundaries(), boundaries);
}

//New snowflake tests for possible boundaries: Needs modifying of more functions.
/*
TEST(BoutMeshTest, GetPossibleBoundariesSnowflake) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  // Minimal valid Snowflake (nype >= 6)
  BoutMeshExposer mesh_SF_1x6(createSnowflake({12, 3, 1, 1, 1, 6, 0, 2}));

  // Larger Snowflake decomposition
  BoutMeshExposer mesh_SF_32x64(createSnowflake({12, 3, 1, 1, 32, 64, 0, 4}));

  std::set<std::string> boundaries{
      "core",
      "pf",
      "sol",
      "upper_target",
      "lower_target",
      "south_pf_outer"
  };

  EXPECT_EQ(mesh_SF_1x6.getPossibleBoundaries(), boundaries);
  EXPECT_EQ(mesh_SF_32x64.getPossibleBoundaries(), boundaries);
}
*/
