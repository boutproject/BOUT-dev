#include "gtest/gtest.h"
#include "test_extras.hxx"

#include "bout/mesh.hxx"
#include "field_factory.hxx"
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
    /// Note: dy = 1/ny so length in Y is 1
    test_coords = std::make_shared<Coordinates>(
        bout::globals::mesh, Field2D{1.0}, Field2D{1.0 / ny}, BoutReal{1.0}, Field2D{1.0},
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

// We can create an instance
TEST_F(SNBTest, CreateInstance) {
  HeatFluxSNB snb;
}

// When there is a temperature gradient the flux is nonzero
TEST_F(SNBTest, FluxNotZero) {
  
  FieldFactory factory;
  auto Te = factory.create3D("5 + cos(y)");
  auto Ne = factory.create3D("1e18 * (1 + 0.5*sin(y))");
  
  HeatFluxSNB snb;
  
  Field3D Div_q_SH;
  Field3D Div_q = snb.divHeatFlux(Te, Ne, &Div_q_SH);

  // Check that flux is not zero
  EXPECT_FALSE(IsFieldEqual(Div_q_SH, 0.0, "RGN_NOBNDRY"));
  EXPECT_FALSE(IsFieldEqual(Div_q, 0.0, "RGN_NOBNDRY"));
}

// When the temperature is constant there is no flux
TEST_F(SNBTest, ZeroTempGradient) {
  
  FieldFactory factory;
  auto Te = factory.create3D("1.5");
  auto Ne = factory.create3D("1e18 * (1 + 0.5*sin(y))");
  
  HeatFluxSNB snb;
  
  Field3D Div_q_SH;
  Field3D Div_q = snb.divHeatFlux(Te, Ne, &Div_q_SH);

  // Check that flux is zero
  EXPECT_TRUE(IsFieldEqual(Div_q_SH, 0.0, "RGN_NOBNDRY"));
  EXPECT_TRUE(IsFieldEqual(Div_q, 0.0, "RGN_NOBNDRY"));
}

// In the collisional limit the SH and SNB fluxes are close
TEST_F(SNBTest, CollisionalLimitClose) {
  
  FieldFactory factory;
  auto Te = factory.create3D("1 + 0.01*sin(y)");
  auto Ne = factory.create3D("1e20 * (1 + 0.5*sin(y))");
  
  HeatFluxSNB snb;
  
  Field3D Div_q_SH;
  Field3D Div_q = snb.divHeatFlux(Te, Ne, &Div_q_SH);

  // Check that flux is zero
  EXPECT_TRUE(IsFieldEqual(Div_q, Div_q_SH, "RGN_NOBNDRY"));
}

// In the collisionless limit the SH and SNB fluxes are different
TEST_F(SNBTest, CollisionlessLimitDifferent) {
  
  FieldFactory factory;
  auto Te = factory.create3D("1e3 + 0.01*sin(y)");
  auto Ne = factory.create3D("1e19 * (1 + 0.5*sin(y))");
  
  HeatFluxSNB snb;
  
  Field3D Div_q_SH;
  Field3D Div_q = snb.divHeatFlux(Te, Ne, &Div_q_SH);

  // Check that fluxes are not equal
  EXPECT_FALSE(IsFieldEqual(Div_q, Div_q_SH, "RGN_NOBNDRY"));
}

// Reversing the temperature gradient reverses the flux
TEST_F(SNBTest, FluxReverses) {
  
  FieldFactory factory;
  auto Te = factory.create3D("10 + 0.01*sin(y)");
  auto Ne = factory.create3D("1e19 * (1 + 0.5*sin(y))");
  
  HeatFluxSNB snb;
  
  Field3D Div_q_SH_1;
  Field3D Div_q_1 = snb.divHeatFlux(Te, Ne, &Div_q_SH_1);

  auto Te2 = factory.create3D("10 - 0.01*sin(y)");

  Field3D Div_q_SH_2;
  Field3D Div_q_2 = snb.divHeatFlux(Te2, Ne, &Div_q_SH_2);
  
  // Check that fluxes are reversed
  EXPECT_TRUE(IsFieldEqual(Div_q_2, -Div_q_1, "RGN_NOBNDRY"));
  EXPECT_TRUE(IsFieldEqual(Div_q_SH_2, -Div_q_SH_1, "RGN_NOBNDRY"));
}

// Spitzer-Harm is independent of density, but SNB is not
TEST_F(SNBTest, DensityDependence) {
  
  FieldFactory factory;
  auto Te = factory.create3D("10 + 0.01*sin(y)");
  auto Ne = factory.create3D("1e19 * (1 + 0.5*sin(y))");
  
  HeatFluxSNB snb;
  
  Field3D Div_q_SH_1;
  Field3D Div_q_1 = snb.divHeatFlux(Te, Ne, &Div_q_SH_1);

  auto Ne2 = factory.create3D("1e17 * (1 + 0.5*sin(y))");

  Field3D Div_q_SH_2;
  Field3D Div_q_2 = snb.divHeatFlux(Te, Ne2, &Div_q_SH_2);

  // SNB result different
  EXPECT_FALSE(IsFieldEqual(Div_q_2, Div_q_1, "RGN_NOBNDRY"));
  
  // Spitzer-Harm not changed
  EXPECT_TRUE(IsFieldEqual(Div_q_SH_2, Div_q_SH_1, "RGN_NOBNDRY"));
}
