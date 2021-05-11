#include "gtest/gtest.h"

#include "boutexception.hxx"
#include "output.hxx"
#include "test_extras.hxx"
#include "unused.hxx"
#include "vector3d.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "bout/mpi_wrapper.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

/// Test fixture to make sure the global mesh is our fake one
class Vector3DTest : public ::testing::Test {
protected:
  Vector3DTest() {
    WithQuietOutput quiet{output_info};
    // Delete any existing mesh
    if (mesh != nullptr) {
      // Delete boundary regions
      for (auto& r : mesh->getBoundaries()) {
        delete r;
      }

      delete mesh;
      mesh = nullptr;
    }
    bout::globals::mpi = new MpiWrapper();
    mesh = new FakeMesh(nx, ny, nz);
    static_cast<FakeMesh*>(mesh)->setCoordinates(nullptr);
    mesh->createDefaultRegions();

    mesh->addBoundary(new BoundaryRegionXIn("core", 1, ny - 2, mesh));
    mesh->addBoundary(new BoundaryRegionXOut("sol", 1, ny - 2, mesh));
    mesh->addBoundary(new BoundaryRegionYUp("upper_target", 1, nx - 2, mesh));
    mesh->addBoundary(new BoundaryRegionYDown("lower_target", 1, nx - 2, mesh));

    static_cast<FakeMesh*>(mesh)->setCoordinates(std::make_shared<Coordinates>(
        mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{1.0}, Field2D{2.0}, Field2D{3.0}, Field2D{4.0}, Field2D{5.0},
        Field2D{6.0}, Field2D{1.0}, Field2D{2.0}, Field2D{3.0}, Field2D{4.0},
        Field2D{5.0}, Field2D{6.0}, Field2D{0.0}, Field2D{0.0}, false));

    delete mesh_staggered;
    mesh_staggered = new FakeMesh(nx, ny, nz);
    mesh_staggered->StaggerGrids = true;
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr);
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_XLOW);
    mesh_staggered->createDefaultRegions();
  }

  ~Vector3DTest() override {
    if (mesh != nullptr) {
      // Delete boundary regions
      for (auto &r : mesh->getBoundaries()) {
        delete r;
      }
    }
    delete mesh;
    mesh = nullptr;
    delete mesh_staggered;
    mesh_staggered = nullptr;
    delete bout::globals::mpi;
    bout::globals::mpi = nullptr;
  }

public:
  static const int nx;
  static const int ny;
  static const int nz;

  Mesh* mesh_staggered = nullptr;
};

const int Vector3DTest::nx = 5;
const int Vector3DTest::ny = 5;
const int Vector3DTest::nz = 3;

TEST_F(Vector3DTest, ApplyBoundaryString) {
  WithQuietOutput quiet{output_info};

  Vector3D v;
  v = 0.0;
  v.applyBoundary("dirichlet(1.0)");

  // boundary cell in x
  EXPECT_DOUBLE_EQ(v.x(0,2,0), 2.0);
  EXPECT_DOUBLE_EQ(v.y(4,2,1), 2.0);
  
  // boundary cell in y
  EXPECT_DOUBLE_EQ(v.x(2,0,2), 2.0);
  EXPECT_DOUBLE_EQ(v.z(2,4,0), 2.0);

  // Middle cell not changed
  EXPECT_DOUBLE_EQ(v.x(2,2,1), 0.0);
}

TEST_F(Vector3DTest, Is3D) {
  Vector3D vector;

  EXPECT_TRUE(vector.is3D());
}

TEST_F(Vector3DTest, BoutRealSize) {
  Vector3D vector;

  EXPECT_EQ(vector.elementSize(), 3);
}

TEST_F(Vector3DTest, TimeDeriv) {
  Vector3D vector;

  auto deriv = vector.timeDeriv();
  EXPECT_NE(&vector, deriv);

  auto deriv2 = vector.timeDeriv();
  EXPECT_EQ(deriv, deriv2);

  EXPECT_EQ(&(ddt(vector)), deriv);
}

TEST_F(Vector3DTest, TimeDerivComponents) {
  Vector3D vector;

  // Make time derivatives for components first, then check we rectify
  vector.x.timeDeriv();
  vector.y.timeDeriv();
  vector.z.timeDeriv();
  vector.timeDeriv();

  EXPECT_EQ(&(ddt(vector).x), &(ddt(vector.x)));
  EXPECT_EQ(&(ddt(vector).y), &(ddt(vector.y)));
  EXPECT_EQ(&(ddt(vector).z), &(ddt(vector.z)));
}

TEST_F(Vector3DTest, SetLocationNonStaggered) {
  Vector3D vector;
  EXPECT_EQ(vector.getLocation(), CELL_CENTRE);
  EXPECT_NO_THROW(vector.setLocation(CELL_CENTRE));
  EXPECT_EQ(vector.getLocation(), CELL_CENTRE);
#if CHECK > 0
  EXPECT_THROW(vector.setLocation(CELL_XLOW), BoutException);
#endif
}

TEST_F(Vector3DTest, SetLocationXLOW) {
  Vector3D vector(mesh_staggered);
  CELL_LOC targetLoc = CELL_XLOW;
  EXPECT_EQ(vector.getLocation(), CELL_CENTRE);
  EXPECT_NO_THROW(vector.setLocation(targetLoc));
  EXPECT_EQ(vector.getLocation(), targetLoc);
  EXPECT_EQ(vector.x.getLocation(), targetLoc);
  EXPECT_EQ(vector.y.getLocation(), targetLoc);
  EXPECT_EQ(vector.z.getLocation(), targetLoc);
}

TEST_F(Vector3DTest, SetLocationYLOW) {
  FakeMesh local_mesh{Vector3DTest::nx,Vector3DTest::ny,Vector3DTest::nz};
  local_mesh.setCoordinates(nullptr);
  local_mesh.StaggerGrids = true;
  local_mesh.setCoordinates(nullptr, CELL_YLOW);
  Vector3D vector(&local_mesh);
  CELL_LOC targetLoc = CELL_YLOW;
  EXPECT_EQ(vector.getLocation(), CELL_CENTRE);
  EXPECT_NO_THROW(vector.setLocation(targetLoc));
  EXPECT_EQ(vector.getLocation(), targetLoc);
  EXPECT_EQ(vector.x.getLocation(), targetLoc);
  EXPECT_EQ(vector.y.getLocation(), targetLoc);
  EXPECT_EQ(vector.z.getLocation(), targetLoc);
}

TEST_F(Vector3DTest, SetLocationZLOW) {
  FakeMesh local_mesh{Vector3DTest::nx,Vector3DTest::ny,Vector3DTest::nz};
  local_mesh.setCoordinates(nullptr);
  local_mesh.StaggerGrids = true;
  local_mesh.setCoordinates(nullptr, CELL_ZLOW);
  Vector3D vector(&local_mesh);
  CELL_LOC targetLoc = CELL_ZLOW;
  EXPECT_EQ(vector.getLocation(), CELL_CENTRE);
  EXPECT_NO_THROW(vector.setLocation(targetLoc));
  EXPECT_EQ(vector.getLocation(), targetLoc);
  EXPECT_EQ(vector.x.getLocation(), targetLoc);
  EXPECT_EQ(vector.y.getLocation(), targetLoc);
  EXPECT_EQ(vector.z.getLocation(), targetLoc);
}

TEST_F(Vector3DTest, SetLocationVSHIFT) {
  FakeMesh local_mesh{Vector3DTest::nx,Vector3DTest::ny,Vector3DTest::nz};
  local_mesh.setCoordinates(nullptr);
  local_mesh.StaggerGrids = true;
  local_mesh.setCoordinates(nullptr, CELL_XLOW);
  local_mesh.setCoordinates(nullptr, CELL_YLOW);
  local_mesh.setCoordinates(nullptr, CELL_ZLOW);
  Vector3D vector(&local_mesh);
  EXPECT_EQ(vector.getLocation(), CELL_CENTRE);
  EXPECT_NO_THROW(vector.setLocation(CELL_VSHIFT));
  EXPECT_EQ(vector.getLocation(), CELL_VSHIFT);
  EXPECT_EQ(vector.x.getLocation(), CELL_XLOW);
  EXPECT_EQ(vector.y.getLocation(), CELL_YLOW);
  EXPECT_EQ(vector.z.getLocation(), CELL_ZLOW);
}

TEST_F(Vector3DTest, SetLocationDEFAULT) {
  Vector3D vector;
  CELL_LOC targetLoc = CELL_CENTRE;
  vector.x.getMesh()->StaggerGrids = true;
  EXPECT_EQ(vector.getLocation(), CELL_CENTRE);
  EXPECT_NO_THROW(vector.setLocation(CELL_DEFAULT));
  EXPECT_EQ(vector.getLocation(), targetLoc);
  EXPECT_EQ(vector.x.getLocation(), targetLoc);
  EXPECT_EQ(vector.y.getLocation(), targetLoc);
  EXPECT_EQ(vector.z.getLocation(), targetLoc);
}

TEST_F(Vector3DTest, AssignFromBoutReal) {
  Vector3D vector;

  vector = 0.0;

  EXPECT_TRUE(IsFieldEqual(vector.x, 0.0));
  EXPECT_TRUE(IsFieldEqual(vector.y, 0.0));
  EXPECT_TRUE(IsFieldEqual(vector.z, 0.0));
}

TEST_F(Vector3DTest, AssignFromVector2D) {
  Vector2D vector1(mesh_staggered);
  Vector3D vector2(mesh_staggered);

  vector1.x = 1.0;
  vector1.y = 2.0;
  vector1.z = 3.0;
  vector1.setLocation(CELL_XLOW);

  vector2 = vector1;

  EXPECT_TRUE(IsFieldEqual(vector2.x, 1.0));
  EXPECT_TRUE(IsFieldEqual(vector2.y, 2.0));
  EXPECT_TRUE(IsFieldEqual(vector2.z, 3.0));
  EXPECT_EQ(vector1.getLocation(), vector2.getLocation());
}

TEST_F(Vector3DTest, AssignFromVector3D) {
  Vector3D vector1(mesh_staggered), vector2(mesh_staggered);

  vector1.x = 1.0;
  vector1.y = 2.0;
  vector1.z = 3.0;
  vector1.setLocation(CELL_XLOW);

  vector2 = vector1;

  EXPECT_TRUE(IsFieldEqual(vector2.x, 1.0));
  EXPECT_TRUE(IsFieldEqual(vector2.y, 2.0));
  EXPECT_TRUE(IsFieldEqual(vector2.z, 3.0));
  EXPECT_EQ(vector1.getLocation(), vector2.getLocation());
}

TEST_F(Vector3DTest, CreateFromVector3D) {
  Vector3D vector1(mesh_staggered);

  vector1.x = 4.0;
  vector1.y = 5.0;
  vector1.z = 6.0;
  vector1.setLocation(CELL_XLOW);

  Vector3D vector2{vector1};

  EXPECT_TRUE(IsFieldEqual(vector2.x, 4.0));
  EXPECT_TRUE(IsFieldEqual(vector2.y, 5.0));
  EXPECT_TRUE(IsFieldEqual(vector2.z, 6.0));
  EXPECT_EQ(vector1.getLocation(), vector2.getLocation());
}

TEST_F(Vector3DTest, UnaryMinus) {
  Vector3D vector1;
  vector1.x = 7.0;
  vector1.y = 8.0;
  vector1.z = 9.0;

  Vector3D vector2;

  vector2 = -vector1;

  EXPECT_TRUE(IsFieldEqual(vector2.x, -7.0));
  EXPECT_TRUE(IsFieldEqual(vector2.y, -8.0));
  EXPECT_TRUE(IsFieldEqual(vector2.z, -9.0));
}

TEST_F(Vector3DTest, AddEqualsVector3D) {
  Vector3D vector1;
  vector1.x = 10.0;
  vector1.y = 11.0;
  vector1.z = 12.0;

  Vector3D vector2;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  vector2 += vector1;

  EXPECT_TRUE(IsFieldEqual(vector2.x, 11.0));
  EXPECT_TRUE(IsFieldEqual(vector2.y, 13.0));
  EXPECT_TRUE(IsFieldEqual(vector2.z, 15.0));
}

TEST_F(Vector3DTest, AddVector3DVector2D) {
  Vector3D vector1;
  vector1.x = 13.0;
  vector1.y = 14.0;
  vector1.z = 15.0;

  Vector2D vector2;
  vector2.x = 4.0;
  vector2.y = 5.0;
  vector2.z = 6.0;

  Vector3D result = vector1 + vector2;

  EXPECT_TRUE(IsFieldEqual(result.x, 17.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 19.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 21.0));
}

TEST_F(Vector3DTest, AddVector3DVector3D) {
  Vector3D vector1;
  vector1.x = 3.0;
  vector1.y = 4.0;
  vector1.z = 5.0;

  Vector3D vector2;
  vector2.x = 4.0;
  vector2.y = 5.0;
  vector2.z = 6.0;

  Vector3D result = vector1 + vector2;

  EXPECT_TRUE(IsFieldEqual(result.x, 7.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 9.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 11.0));
}

TEST_F(Vector3DTest, MinusEqualsVector3D) {
  Vector3D vector1;
  vector1.x = 100.0;
  vector1.y = 101.0;
  vector1.z = 102.0;

  Vector3D vector2;
  vector2.x = 3.0;
  vector2.y = 2.0;
  vector2.z = 1.0;

  vector2 -= vector1;

  EXPECT_TRUE(IsFieldEqual(vector2.x, -97.0));
  EXPECT_TRUE(IsFieldEqual(vector2.y, -99.0));
  EXPECT_TRUE(IsFieldEqual(vector2.z, -101.0));
}

TEST_F(Vector3DTest, MinusVector3DVector2D) {
  Vector3D vector1;
  vector1.x = 0.0;
  vector1.y = 1.0;
  vector1.z = 2.0;

  Vector2D vector2;
  vector2.x = 3.0;
  vector2.y = 2.0;
  vector2.z = 1.0;

  Vector3D result = vector1 - vector2;

  EXPECT_TRUE(IsFieldEqual(result.x, -3.0));
  EXPECT_TRUE(IsFieldEqual(result.y, -1.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 1.0));
}

TEST_F(Vector3DTest, MinusVector3DVector3D) {
  Vector3D vector1;
  vector1.x = 10.0;
  vector1.y = 11.0;
  vector1.z = 12.0;

  Vector3D vector2;
  vector2.x = 3.0;
  vector2.y = 2.0;
  vector2.z = 1.0;

  Vector3D result = vector1 - vector2;

  EXPECT_TRUE(IsFieldEqual(result.x, 7.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 9.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 11.0));
}

TEST_F(Vector3DTest, MultiplyEqualsBoutReal) {
  Vector3D vector;
  vector.x = 4.0;
  vector.y = 5.0;
  vector.z = 6.0;

  BoutReal real {4.0};

  vector *= real;

  EXPECT_TRUE(IsFieldEqual(vector.x, 16.0));
  EXPECT_TRUE(IsFieldEqual(vector.y, 20.0));
  EXPECT_TRUE(IsFieldEqual(vector.z, 24.0));
}

TEST_F(Vector3DTest, MultiplyEqualsField2D) {
  Vector3D vector;
  vector.x = 4.0;
  vector.y = 5.0;
  vector.z = 6.0;

  Field2D field{40.0};

  vector *= field;

  EXPECT_TRUE(IsFieldEqual(vector.x, 160.0));
  EXPECT_TRUE(IsFieldEqual(vector.y, 200.0));
  EXPECT_TRUE(IsFieldEqual(vector.z, 240.0));
}

TEST_F(Vector3DTest, MultiplyVector3DBoutReal) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  BoutReal real {2.0};

  Vector3D result = vector * real;

  EXPECT_TRUE(IsFieldEqual(result.x, 2.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 4.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 6.0));
}

TEST_F(Vector3DTest, MultiplyVector3DField2D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field2D field{3.0};

  Vector3D result = vector * field;

  EXPECT_TRUE(IsFieldEqual(result.x, 3.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 6.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 9.0));
}

TEST_F(Vector3DTest, MultiplyVector3DField3D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field3D field{4.0};

  Vector3D result = vector * field;

  EXPECT_TRUE(IsFieldEqual(result.x, 4.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 8.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 12.0));
}

TEST_F(Vector3DTest, MultiplyBoutRealVector3D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  BoutReal real {2.0};

  Vector3D result = real * vector;

  EXPECT_TRUE(IsFieldEqual(result.x, 2.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 4.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 6.0));
}

TEST_F(Vector3DTest, MultiplyField2DVector3D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field2D field{3.0};

  Vector3D result = field * vector;

  EXPECT_TRUE(IsFieldEqual(result.x, 3.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 6.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 9.0));
}

TEST_F(Vector3DTest, MultiplyField3DVector3D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field3D field{4.0};

  Vector3D result = field * vector;

  EXPECT_TRUE(IsFieldEqual(result.x, 4.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 8.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 12.0));
}

TEST_F(Vector3DTest, DivideEqualsBoutReal) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  BoutReal real{10.0};

  vector /= real;

  EXPECT_TRUE(IsFieldEqual(vector.x, 0.1));
  EXPECT_TRUE(IsFieldEqual(vector.y, 0.2));
  EXPECT_TRUE(IsFieldEqual(vector.z, 0.3));
}

TEST_F(Vector3DTest, DivideEqualsConstField2D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field2D field{5.0};

  vector /= field;

  EXPECT_TRUE(IsFieldEqual(vector.x, 0.2));
  EXPECT_TRUE(IsFieldEqual(vector.y, 0.4));
  EXPECT_TRUE(IsFieldEqual(vector.z, 0.6));
}

TEST_F(Vector3DTest, DivideVector3DBoutReal) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  BoutReal real {2.0};

  Vector3D result = vector / real;

  EXPECT_TRUE(IsFieldEqual(result.x, 0.5));
  EXPECT_TRUE(IsFieldEqual(result.y, 1.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 1.5));
}

TEST_F(Vector3DTest, DivideVector3DField2D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field2D field{4.0};

  Vector3D result = vector / field;

  EXPECT_TRUE(IsFieldEqual(result.x, 0.25));
  EXPECT_TRUE(IsFieldEqual(result.y, 0.5));
  EXPECT_TRUE(IsFieldEqual(result.z, 0.75));
}

TEST_F(Vector3DTest, DivideVector3DField3D) {
  Vector3D vector;
  vector.x = 2.0;
  vector.y = 4.0;
  vector.z = 6.0;

  Field3D field{2.0};

  Vector3D result = vector / field;

  EXPECT_TRUE(IsFieldEqual(result.x, 1.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 2.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 3.0));
}

TEST_F(Vector3DTest, ToCovariant) {
  Vector3D vector;
  vector.covariant = false;
  vector.x = 2.0;
  vector.y = 4.0;
  vector.z = 6.0;

  vector.toCovariant();

  EXPECT_TRUE(IsFieldEqual(vector.x, 48.0));
  EXPECT_TRUE(IsFieldEqual(vector.y, 52.0));
  EXPECT_TRUE(IsFieldEqual(vector.z, 52.0));
}

TEST_F(Vector3DTest, ToContravariant) {
  Vector3D vector;
  vector.covariant = true;
  vector.x = 2.0;
  vector.y = 4.0;
  vector.z = 6.0;

  vector.toContravariant();

  EXPECT_TRUE(IsFieldEqual(vector.x, 48.0));
  EXPECT_TRUE(IsFieldEqual(vector.y, 52.0));
  EXPECT_TRUE(IsFieldEqual(vector.z, 52.0));
}

TEST_F(Vector3DTest, Cross3D3D) {
  Vector3D vector1;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  Vector3D vector2;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  auto result = cross(vector1, vector2);

  EXPECT_TRUE(IsFieldEqual(result.x, 0.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 0.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 0.0));
}

TEST_F(Vector3DTest, Cross3D2D) {
  Vector3D vector1;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  Vector2D vector2;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  auto result = cross(vector1, vector2);

  EXPECT_TRUE(IsFieldEqual(result.x, 0.0));
  EXPECT_TRUE(IsFieldEqual(result.y, 0.0));
  EXPECT_TRUE(IsFieldEqual(result.z, 0.0));
}

TEST_F(Vector3DTest, Dot3D3DCoContra) {
  Vector3D vector1;
  vector1.covariant = false;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  Vector3D vector2;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  auto result = vector1 * vector2;

  EXPECT_TRUE(IsFieldEqual(result, 28.0));
}

TEST_F(Vector3DTest, Dot3D2DCoContra) {
  Vector3D vector1;
  vector1.covariant = false;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  Vector2D vector2;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  auto result = vector1 * vector2;

  EXPECT_TRUE(IsFieldEqual(result, 28.0));
}

TEST_F(Vector3DTest, Dot3D3DCoCo) {
  Vector3D vector1;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  Vector3D vector2;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  auto result = vector1 * vector2;

  EXPECT_TRUE(IsFieldEqual(result, 308.0));
}

TEST_F(Vector3DTest, Dot3D2DCoCo) {
  Vector3D vector1;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  Vector2D vector2;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  auto result = vector1 * vector2;

  EXPECT_TRUE(IsFieldEqual(result, 308.0));
}

TEST_F(Vector3DTest, Dot3D3DContraContra) {
  Vector3D vector1;
  vector1.covariant = false;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  Vector3D vector2;
  vector2.covariant = false;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  auto result = vector1 * vector2;

  EXPECT_TRUE(IsFieldEqual(result, 308.0));
}

TEST_F(Vector3DTest, Dot3D2DContraContra) {
  Vector3D vector1;
  vector1.covariant = false;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  Vector2D vector2;
  vector2.covariant = false;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  auto result = vector1 * vector2;

  EXPECT_TRUE(IsFieldEqual(result, 308.0));
}

TEST_F(Vector3DTest, AbsCo) {
  Vector3D vector1;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  auto result = abs(vector1);

  EXPECT_TRUE(IsFieldEqual(result, 24.819347291981714));
}

TEST_F(Vector3DTest, AbsContra) {
  Vector3D vector1;
  vector1.covariant = false;
  vector1.x = 2.0;
  vector1.y = 4.0;
  vector1.z = 6.0;

  auto result = abs(vector1);

  EXPECT_TRUE(IsFieldEqual(result, 24.819347291981714));
}
