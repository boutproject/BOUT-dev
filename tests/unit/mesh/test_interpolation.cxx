#include "gtest/gtest.h"

#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "interpolation.hxx"
#include "output.hxx"
#include "test_extras.hxx"

////// delete these
#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "field3d.hxx"
#include "unused.hxx"
#include "utils.hxx"
#include <cmath>
#include <set>
#include <vector>
///////

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

/// Test fixture to make sure the global mesh is our fake one
class Field3DInterpToTest : public ::testing::Test {
protected:
  Field3DInterpToTest() : input(mesh) {
    input = 0.;
    input(2, 2, 2) = 2.;
    input(1, 2, 2) = 1.8;
    input(0, 2, 2) = 1.7;
    input(3, 2, 2) = 1.3;
    input(4, 2, 2) = 1.5;
    input(2, 1, 2) = 3.8;
    input(2, 0, 2) = 3.7;
    input(2, 3, 2) = 3.3;
    input(2, 4, 2) = 3.5;
    input(2, 2, 1) = 5.8;
    input(2, 2, 0) = 5.7;
    input(2, 2, 3) = 5.3;
    input(2, 2, 4) = 5.5;
  }

  static void SetUpTestCase() {
    WithQuietOutput quiet_info{output_info};
    WithQuietOutput quiet_warn{output_warn};
    delete mesh;
    mesh = new FakeMesh(nx, ny, nz);
    mesh->StaggerGrids = true;
    mesh->xstart = 2;
    mesh->ystart = 2;
    mesh->xend = nx - 3;
    mesh->yend = ny - 3;

    mesh->createDefaultRegions();

    // We need Coordinates so a parallel transform is available as
    // FieldFactory::create3D wants to un-field-align the result
    for (const auto& location
        : std::list<CELL_LOC>{CELL_CENTRE, CELL_XLOW, CELL_YLOW, CELL_ZLOW}) {

      static_cast<FakeMesh*>(mesh)->setCoordinates(nullptr, location);
      static_cast<FakeMesh*>(mesh)->setCoordinates(std::make_shared<Coordinates>(
          mesh, Field2D{1.0, mesh}, Field2D{1.0, mesh}, BoutReal{1.0}, Field2D{1.0, mesh},
          Field2D{0.0, mesh}, Field2D{1.0, mesh}, Field2D{1.0, mesh}, Field2D{1.0, mesh},
          Field2D{0.0, mesh}, Field2D{0.0, mesh}, Field2D{0.0, mesh}, Field2D{1.0, mesh},
          Field2D{1.0, mesh}, Field2D{1.0, mesh}, Field2D{0.0, mesh}, Field2D{0.0, mesh},
          Field2D{0.0, mesh}, Field2D{0.0, mesh}, Field2D{0.0, mesh}, false),
          location);
      mesh->getCoordinates(location)->setParallelTransform(
          bout::utils::make_unique<ParallelTransformIdentity>(*mesh));
    }
  }

  static void TearDownTestCase() {
    delete mesh;
    mesh = nullptr;
  }

  Field3D input;

public:
  static const int nx;
  static const int ny;
  static const int nz;
};

const int Field3DInterpToTest::nx = 5;
const int Field3DInterpToTest::ny = 5;
const int Field3DInterpToTest::nz = 7;

TEST_F(Field3DInterpToTest, CellCentreToXlow) {

  Field3D output = Field3D(mesh);

  // CELL_CENTRE -> CELL_XLOW
  input.setLocation(CELL_CENTRE);
  output = interp_to(input, CELL_XLOW);
  EXPECT_TRUE(output.getLocation() == CELL_XLOW);
  EXPECT_NEAR(output(2, 2, 2), 1.95, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellCentreToXlowNoBndry) {

  Field3D output = Field3D(mesh);

  // CELL_CENTRE -> CELL_XLOW
  input.setLocation(CELL_CENTRE);
  output = interp_to(input, CELL_XLOW, "RGN_NOBNDRY");
  EXPECT_TRUE(output.getLocation() == CELL_XLOW);
  EXPECT_NEAR(output(2, 2, 2), 1.95, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellXlowToCentre) {

  Field3D output = Field3D(mesh);

  // CELL_XLOW -> CELL_CENTRE
  input.setLocation(CELL_XLOW);
  output = interp_to(input, CELL_CENTRE);
  EXPECT_TRUE(output.getLocation() == CELL_CENTRE);
  EXPECT_NEAR(output(2, 2, 2), 1.65, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellXlowToCentreNoBndry) {

  Field3D output = Field3D(mesh);

  // CELL_XLOW -> CELL_CENTRE
  input.setLocation(CELL_XLOW);
  output = interp_to(input, CELL_CENTRE, "RGN_NOBNDRY");
  EXPECT_TRUE(output.getLocation() == CELL_CENTRE);
  EXPECT_NEAR(output(2, 2, 2), 1.65, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellCentreToYlow) {

  Field3D output = Field3D(mesh);

  // CELL_CENTRE -> CELL_YLOW
  input.setLocation(CELL_CENTRE);
  output = interp_to(input, CELL_YLOW);
  EXPECT_TRUE(output.getLocation() == CELL_YLOW);
  EXPECT_NEAR(output(2, 2, 2), 2.825, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellCentreToYlowNoBndry) {

  Field3D output = Field3D(mesh);

  // CELL_CENTRE -> CELL_YLOW
  input.setLocation(CELL_CENTRE);
  output = interp_to(input, CELL_YLOW, "RGN_NOBNDRY");
  EXPECT_TRUE(output.getLocation() == CELL_YLOW);
  EXPECT_NEAR(output(2, 2, 2), 2.825, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellYlowToCentre) {

  Field3D output = Field3D(mesh);

  // CELL_YLOW -> CELL_CENTRE
  input.setLocation(CELL_YLOW);
  output = interp_to(input, CELL_CENTRE);
  EXPECT_TRUE(output.getLocation() == CELL_CENTRE);
  EXPECT_NEAR(output(2, 2, 2), 2.525, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellYlowToCentreNoBndry) {

  Field3D output = Field3D(mesh);

  // CELL_YLOW -> CELL_CENTRE
  input.setLocation(CELL_YLOW);
  output = interp_to(input, CELL_CENTRE, "RGN_NOBNDRY");
  EXPECT_TRUE(output.getLocation() == CELL_CENTRE);
  EXPECT_NEAR(output(2, 2, 2), 2.525, 1.e-15);
}

TEST_F(Field3DInterpToTest, AlignedCellCentreToYlow) {

  Field3D output = Field3D(mesh);

  // CELL_CENTRE -> CELL_YLOW
  input.setLocation(CELL_CENTRE).setDirectionY(YDirectionType::Aligned);
  output = interp_to(input, CELL_YLOW);
  EXPECT_TRUE(output.getLocation() == CELL_YLOW);
  EXPECT_TRUE(output.getDirectionY() == YDirectionType::Aligned);
  EXPECT_NEAR(output(2, 2, 2), 2.825, 1.e-15);
}

TEST_F(Field3DInterpToTest, AlignedCellCentreToYlowNoBndry) {

  Field3D output = Field3D(mesh);

  // CELL_CENTRE -> CELL_YLOW
  input.setLocation(CELL_CENTRE).setDirectionY(YDirectionType::Aligned);
  output = interp_to(input, CELL_YLOW, "RGN_NOBNDRY");
  EXPECT_TRUE(output.getLocation() == CELL_YLOW);
  EXPECT_TRUE(output.getDirectionY() == YDirectionType::Aligned);
  EXPECT_NEAR(output(2, 2, 2), 2.825, 1.e-15);
}

TEST_F(Field3DInterpToTest, AlignedCellYlowToCentre) {

  Field3D output = Field3D(mesh);

  // CELL_YLOW -> CELL_CENTRE
  input.setLocation(CELL_YLOW).setDirectionY(YDirectionType::Aligned);
  output = interp_to(input, CELL_CENTRE);
  EXPECT_TRUE(output.getLocation() == CELL_CENTRE);
  EXPECT_TRUE(output.getDirectionY() == YDirectionType::Aligned);
  EXPECT_NEAR(output(2, 2, 2), 2.525, 1.e-15);
}

TEST_F(Field3DInterpToTest, AlignedCellYlowToCentreNoBndry) {

  Field3D output = Field3D(mesh);

  // CELL_YLOW -> CELL_CENTRE
  input.setLocation(CELL_YLOW).setDirectionY(YDirectionType::Aligned);
  output = interp_to(input, CELL_CENTRE, "RGN_NOBNDRY");
  EXPECT_TRUE(output.getLocation() == CELL_CENTRE);
  EXPECT_TRUE(output.getDirectionY() == YDirectionType::Aligned);
  EXPECT_NEAR(output(2, 2, 2), 2.525, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellCentreToZlow) {

  Field3D output = Field3D(mesh);

  // CELL_CENTRE -> CELL_ZLOW
  input.setLocation(CELL_CENTRE);
  output = interp_to(input, CELL_ZLOW);
  EXPECT_TRUE(output.getLocation() == CELL_ZLOW);
  EXPECT_NEAR(output(2, 2, 2), 3.7, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellCentreToZlowNoBndry) {

  Field3D output = Field3D(mesh);

  // CELL_CENTRE -> CELL_ZLOW
  input.setLocation(CELL_CENTRE);
  output = interp_to(input, CELL_ZLOW, "RGN_NOBNDRY");
  EXPECT_TRUE(output.getLocation() == CELL_ZLOW);
  EXPECT_NEAR(output(2, 2, 2), 3.7, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellZlowToCentre) {

  Field3D output = Field3D(mesh);

  // CELL_XLOW -> CELL_CENTRE
  input.setLocation(CELL_ZLOW);
  output = interp_to(input, CELL_CENTRE);
  EXPECT_TRUE(output.getLocation() == CELL_CENTRE);
  EXPECT_NEAR(output(2, 2, 2), 3.4, 1.e-15);
}

TEST_F(Field3DInterpToTest, CellZlowToCentreNoBndry) {

  Field3D output = Field3D(mesh);

  // CELL_XLOW -> CELL_CENTRE
  input.setLocation(CELL_ZLOW);
  output = interp_to(input, CELL_CENTRE, "RGN_NOBNDRY");
  EXPECT_TRUE(output.getLocation() == CELL_CENTRE);
  EXPECT_NEAR(output(2, 2, 2), 3.4, 1.e-15);
}
