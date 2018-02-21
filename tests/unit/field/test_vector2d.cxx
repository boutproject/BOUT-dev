#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "bout/boutexception.hxx"
#include "bout/vector2d.hxx"
#include "test_extras.hxx"
#include "bout/unused.hxx"

/// Global mesh
extern Mesh *mesh;

/// Test fixture to make sure the global mesh is our fake one
class Vector2DTest : public ::testing::Test {
protected:
  static void SetUpTestCase() {
    // Delete any existing mesh
    if (mesh != nullptr) {
      // Delete boundary regions
      for (auto &r : mesh->getBoundaries()) {
        delete r;
      }

      delete mesh;
      mesh = nullptr;
    }
    mesh = new FakeMesh(nx, ny, nz);

    mesh->addBoundary(new BoundaryRegionXIn("core", 1, ny - 2));
    mesh->addBoundary(new BoundaryRegionXOut("sol", 1, ny - 2));
    mesh->addBoundary(new BoundaryRegionYUp("upper_target", 1, nx - 2));
    mesh->addBoundary(new BoundaryRegionYDown("lower_target", 1, nx - 2));
  }

  static void TearDownTestCase() {
    if (mesh != nullptr) {
      // Delete boundary regions
      for (auto &r : mesh->getBoundaries()) {
        delete r;
      }
      delete mesh;
      mesh = nullptr;
    }
  }

public:
  static const int nx;
  static const int ny;
  static const int nz;
};

const int Vector2DTest::nx = 5;
const int Vector2DTest::ny = 5;
const int Vector2DTest::nz = 1;

TEST_F(Vector2DTest, ApplyBoundaryString) {
  Vector2D v;
  v = 0.0;
  v.applyBoundary("dirichlet(1.0)");

  // boundary cell in x
  EXPECT_DOUBLE_EQ(v.x(0,2), 2.0);
  EXPECT_DOUBLE_EQ(v.y(4,2), 2.0);
  
  // boundary cell in y
  EXPECT_DOUBLE_EQ(v.x(2,0), 2.0);
  EXPECT_DOUBLE_EQ(v.z(2,4), 2.0);

  // Middle cell not changed
  EXPECT_DOUBLE_EQ(v.x(2,2), 0.0);
}
