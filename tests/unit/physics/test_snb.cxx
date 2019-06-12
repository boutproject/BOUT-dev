#include "gtest/gtest.h"
#include "test_extras.hxx"

#include "bout/mesh.hxx"

#include "bout/snb.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;


/// Test fixture to make sure the global mesh is our fake
/// one. Also initialize the global mesh_staggered for use in tests with
/// staggering. Multiple tests have exactly the same fixture, so use a type
/// alias to make a new test:
///
///     using MyTest = FakeParallelMeshFixture;
class FakeParallelMeshFixture : public ::testing::Test {
public:
  FakeParallelMeshFixture() {
    WithQuietOutput quiet{output_info};

    delete bout::globals::mesh;
    bout::globals::mesh = new FakeMesh(nx, ny, nz);
    bout::globals::mesh->createDefaultRegions();
    static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(nullptr);
    test_coords = std::make_shared<Coordinates>(
        bout::globals::mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
        false);
    static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(test_coords);
    static_cast<FakeMesh*>(bout::globals::mesh)->setGridDataSource(
        new FakeGridDataSource());
    // May need a ParallelTransform to create fields, because create3D calls
    // fromFieldAligned
    test_coords->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));

    delete mesh_staggered;
    mesh_staggered = new FakeMesh(nx, ny, nz);
    mesh_staggered->StaggerGrids = true;
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr);
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_XLOW);
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_YLOW);
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_ZLOW);
    mesh_staggered->createDefaultRegions();

    test_coords_staggered = std::make_shared<Coordinates>(
        mesh_staggered, Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        BoutReal{1.0}, Field2D{1.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered}, false);
    test_coords_staggered->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*mesh_staggered));
  }

  virtual ~FakeParallelMeshFixture() {
    delete bout::globals::mesh;
    bout::globals::mesh = nullptr;
    delete mesh_staggered;
    mesh_staggered = nullptr;
  }

  static constexpr int nx = 1;
  static constexpr int ny = 16;
  static constexpr int nz = 1;

  Mesh* mesh_staggered = nullptr;

  std::shared_ptr<Coordinates> test_coords{nullptr};
  std::shared_ptr<Coordinates> test_coords_staggered{nullptr};
};

using SNBTest = FakeParallelMeshFixture;

TEST_F(SNBTest, createInstance) {
  HeatFluxSNB snb;
}

