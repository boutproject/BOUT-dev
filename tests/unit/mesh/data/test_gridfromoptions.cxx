#include "gtest/gtest.h"

#include "options.hxx"
#include "output.hxx"
#include "test_extras.hxx"
#include "bout/constants.hxx"
#include "bout/griddata.hxx"
#include "bout/mesh.hxx"

#include <numeric>
#include <string>
#include <vector>

// The unit tests use the global mesh
using namespace bout::globals;

class GridFromOptionsTest : public ::testing::Test {
public:
  GridFromOptionsTest() {

    mesh_from_options.StaggerGrids = true;
    mesh_from_options.xstart = 2;
    mesh_from_options.xend = nx - 3;
    mesh_from_options.ystart = 2;
    mesh_from_options.yend = ny - 3;

    mesh_from_options.createDefaultRegions();

    griddata = new GridFromOptions(&options);
    mesh_from_options.setGridDataSource(griddata);

    output_info.disable();
    output_progress.disable();
    output_warn.disable();
    options["f"] = expected_string;

    // modify mesh section in global options
    options["dx"] = "1.";
    options["dy"] = "1.";
    options["g11"] = expected_string + " + 5.";
    options["g22"] = expected_string + " + 4.";
    options["g33"] = expected_string + " + 3.";
    options["g12"] = expected_string + " + 2.";
    options["g13"] = expected_string + " + 1.";
    options["g23"] = expected_string;

    mesh_from_options.getCoordinates();

    // We need a parallel transform as FieldFactory::create3D wants to
    // un-field-align the result
    mesh_from_options.getCoordinates()->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(mesh_from_options));

    expected_2d = makeField<Field2D>(
        [](Field2D::ind_type& index) {
          return index.x() + (TWOPI * index.y()) + (TWOPI * index.z() / nz) + 3;
        },
        &mesh_from_options);

    expected_3d = makeField<Field3D>(
        [](Field3D::ind_type& index) {
          return index.x() + (TWOPI * index.y()) + (TWOPI * index.z() / nz) + 3;
        },
        &mesh_from_options);
  }

  ~GridFromOptionsTest() override {
    Options::cleanup();
    output_info.enable();
    output_progress.enable();
    output_warn.enable();
    // note GridFromOptions* griddata will be deleted by the ~Mesh() destructor
  }

  static const int nx{9};
  static const int ny{11};
  static const int nz{5};

  std::shared_ptr<Coordinates> test_coords;
  Options options;
  GridFromOptions* griddata{nullptr};
  std::string expected_string{"x + y + z + 3"};
  Field2D expected_2d;
  Field3D expected_3d;
  FakeMesh mesh_from_options{nx, ny, nz};
};

// higher tolerance used when field values are ~50
static constexpr BoutReal this_tolerance{1e-13};

TEST_F(GridFromOptionsTest, HasVar) {
  EXPECT_TRUE(griddata->hasVar("f"));
  EXPECT_FALSE(griddata->hasVar("non-existent"));
}

TEST_F(GridFromOptionsTest, GetString) {
  std::string result{"wrong"};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f"));
  EXPECT_EQ(result, expected_string);
}

TEST_F(GridFromOptionsTest, GetStringNone) {
  std::string result{"wrong"};

  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent"));
  EXPECT_EQ(result, std::string{});
}

TEST_F(GridFromOptionsTest, GetInt) {
  int result{-1};
  int expected{3};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f"));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetIntNone) {
  int result{-1};
  int expected{0};

  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent"));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetBoutReal) {
  BoutReal result{-1.};
  BoutReal expected{3.};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f"));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetBoutRealNone) {
  BoutReal result{-1.};
  BoutReal expected{0};

  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent"));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetField2D) {
  Field2D result{&mesh_from_options};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f"));
  EXPECT_TRUE(IsFieldEqual(result, expected_2d));
}

TEST_F(GridFromOptionsTest, GetField2DNone) {
  Field2D result{&mesh_from_options};
  BoutReal expected{0.};

  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent"));
  EXPECT_TRUE(IsFieldEqual(result, expected));
}

TEST_F(GridFromOptionsTest, GetField2DNoneWithDefault) {
  Field2D result{&mesh_from_options};
  BoutReal default_value{-32};

  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent", default_value));
  EXPECT_TRUE(IsFieldEqual(result, default_value));
}

TEST_F(GridFromOptionsTest, GetField3D) {
  Field3D result{&mesh_from_options};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f"));
  EXPECT_TRUE(IsFieldEqual(result, expected_3d));

  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent"));
  EXPECT_TRUE(IsFieldEqual(result, 0.));

  BoutReal default_value{-64};
  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent", default_value));
  EXPECT_TRUE(IsFieldEqual(result, default_value));
}

TEST_F(GridFromOptionsTest, GetField3DNone) {
  Field3D result{&mesh_from_options};
  BoutReal expected{0.};

  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent"));
  EXPECT_TRUE(IsFieldEqual(result, expected));
}

TEST_F(GridFromOptionsTest, GetField3DNoneWithDefault) {
  Field3D result{&mesh_from_options};
  BoutReal default_value{-64};

  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent", default_value));
  EXPECT_TRUE(IsFieldEqual(result, default_value));
}

TEST_F(GridFromOptionsTest, GetVectorInt) {
  // Getting a vector<int> from GridFromOptions is not currently implemented
  std::vector<int> result{};
  //std::vector<int> expected{3, 3, 3};

  EXPECT_THROW(griddata->get(&mesh_from_options, result, "f", 3), BoutException);
  //EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorIntNone) {
  std::vector<int> result{-1, -1, -1};
  std::vector<int> expected{};

  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent", 3));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealX) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3., 4., 5., 6., 7., 8., 9., 10., 11.};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f", nx));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealXOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{4., 5., 6., 7., 8., 9., 10., 11., 12.};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f", nx, 1, GridDataSource::Direction::X));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealXMeshOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{2., 3., 4., 5., 6., 7., 8., 9., 10.};

  mesh_from_options.OffsetX = 1;
  mesh_from_options.OffsetY = 100;
  mesh_from_options.OffsetZ = 100;

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f", nx, 0, GridDataSource::Direction::X));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealXNone) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> default_expected{};
  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent", nx));
  EXPECT_EQ(result, default_expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealY) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3., 3. + TWOPI, 3. + (2. * TWOPI), 3. + (3. * TWOPI),
                                 3. + (4. * TWOPI), 3. + (5. * TWOPI), 3. + (6. * TWOPI),
                                 3. + (7. * TWOPI), 3. + (8. * TWOPI), 3. + (9. * TWOPI),
                                 3. + (10. * TWOPI)};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f", ny, 0, GridDataSource::Direction::Y));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealYOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3. + TWOPI, 3. + (2. * TWOPI), 3. + (3. * TWOPI),
                                 3. + (4. * TWOPI), 3. + (5. * TWOPI), 3. + (6. * TWOPI),
                                 3. + (7. * TWOPI), 3. + (8. * TWOPI), 3. + (9. * TWOPI),
                                 3. + (10. * TWOPI), 3. + (11. * TWOPI)};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f", ny, 1, GridDataSource::Direction::Y));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealYMeshOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3. - TWOPI, 3., 3. + TWOPI, 3. + (2. * TWOPI),
                                 3. + (3. * TWOPI), 3. + (4. * TWOPI), 3. + (5. * TWOPI),
                                 3. + (6. * TWOPI), 3. + (7. * TWOPI), 3. + (8. * TWOPI),
                                 3. + (9. * TWOPI)};

  mesh_from_options.OffsetX = 100;
  mesh_from_options.OffsetY = 1;
  mesh_from_options.OffsetZ = 100;

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f", ny, 0, GridDataSource::Direction::Y));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealYNone) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> default_expected{};
  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent", ny));
  EXPECT_EQ(result, default_expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealZ) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3.,
                                 3. + (1. * TWOPI / nz),
                                 3. + (2. * TWOPI / nz),
                                 3. + (3. * TWOPI / nz),
                                 3. + (4. * TWOPI / nz)};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f", nz, 0, GridDataSource::Direction::Z));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealZOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3. + (1. * TWOPI / nz), 3. + (2. * TWOPI / nz),
                                 3. + (3. * TWOPI / nz), 3. + (4. * TWOPI / nz),
                                 3. + (5. * TWOPI / nz)};

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f", nz, 1, GridDataSource::Direction::Z));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealZMeshOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3. + (-1. * TWOPI / nz), 3.,
                                 3. + (1. * TWOPI / nz),  3. + (2. * TWOPI / nz),
                                 3. + (3. * TWOPI / nz)};

  mesh_from_options.OffsetX = 100;
  mesh_from_options.OffsetY = 100;
  mesh_from_options.OffsetZ = 1;

  EXPECT_TRUE(griddata->get(&mesh_from_options, result, "f", nz, 0, GridDataSource::Direction::Z));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealZNone) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> default_expected{};
  EXPECT_FALSE(griddata->get(&mesh_from_options, result, "non-existent", nz));
  EXPECT_EQ(result, default_expected);
}

TEST_F(GridFromOptionsTest, CoordinatesCentre) {
  auto coords = mesh_from_options.getCoordinates();

  mesh_from_options.communicate(expected_2d);

  EXPECT_TRUE(IsFieldEqual(coords->g11, expected_2d + 5.));
  EXPECT_TRUE(IsFieldEqual(coords->g22, expected_2d + 4.));
  EXPECT_TRUE(IsFieldEqual(coords->g33, expected_2d + 3.));
  EXPECT_TRUE(IsFieldEqual(coords->g12, expected_2d + 2.));
  EXPECT_TRUE(IsFieldEqual(coords->g13, expected_2d + 1.));
  EXPECT_TRUE(IsFieldEqual(coords->g23, expected_2d));
}

TEST_F(GridFromOptionsTest, CoordinatesZlow) {
  auto coords = mesh_from_options.getCoordinates(CELL_ZLOW);

  mesh_from_options.communicate(expected_2d);

  EXPECT_TRUE(IsFieldEqual(coords->g11, expected_2d + 5.));
  EXPECT_TRUE(IsFieldEqual(coords->g22, expected_2d + 4.));
  EXPECT_TRUE(IsFieldEqual(coords->g33, expected_2d + 3.));
  EXPECT_TRUE(IsFieldEqual(coords->g12, expected_2d + 2.));
  EXPECT_TRUE(IsFieldEqual(coords->g13, expected_2d + 1.));
  EXPECT_TRUE(IsFieldEqual(coords->g23, expected_2d));
}

TEST_F(GridFromOptionsTest, CoordinatesXlowInterp) {
  // *_xlow fields not present in options, Coordinates will be interpolated
  // from CELL_CENTRE

  // make the mesh have boundaries to avoid NaNs in guard cells after interpolating
  mesh_from_options.createBoundaries();

  auto coords = mesh_from_options.getCoordinates(CELL_XLOW);

  Field2D expected_xlow = makeField<Field2D>(
      [](Field2D::ind_type& index) {
        return index.x() - 0.5 + (TWOPI * index.y()) + (TWOPI * index.z() / nz) + 3;
      },
      &mesh_from_options);

  mesh_from_options.communicate(expected_xlow);

  EXPECT_TRUE(IsFieldEqual(coords->g11, expected_xlow + 5., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g22, expected_xlow + 4., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g33, expected_xlow + 3., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g12, expected_xlow + 2., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g13, expected_xlow + 1., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g23, expected_xlow, "RGN_NOBNDRY", this_tolerance));
}

TEST_F(GridFromOptionsTest, CoordinatesXlowRead) {
  // *_xlow fields added to options, will be read to initialise Coordinates

  // Note '(9 - x)' here because FakeMesh::GlobalX(int jx) returns jx, not a
  // global position between 0 and 1 (in grid cells, <0 or >1 in boundaries),
  // like a Mesh is supposed to.
  std::string expected_string_xlow{"(9 - x) + y + 3"};

  // modify mesh section in global options
  options["dx_xlow"] = "1.";
  options["dy_xlow"] = "1.";
  options["g11_xlow"] = expected_string_xlow + " + 5.";
  options["g22_xlow"] = expected_string_xlow + " + 4.";
  options["g33_xlow"] = expected_string_xlow + " + 3.";
  options["g12_xlow"] = expected_string_xlow + " + 2.";
  options["g13_xlow"] = expected_string_xlow + " + 1.";
  options["g23_xlow"] = expected_string_xlow;

  auto coords = mesh_from_options.getCoordinates(CELL_XLOW);

  Field2D expected_xlow = makeField<Field2D>(
      [](Field2D::ind_type& index) {
        return (nx - index.x()) + (TWOPI * index.y()) + (TWOPI * index.z() / nz) + 3;
      },
      &mesh_from_options);

  mesh_from_options.communicate(expected_xlow);

  EXPECT_TRUE(IsFieldEqual(coords->g11, expected_xlow + 5.));
  EXPECT_TRUE(IsFieldEqual(coords->g22, expected_xlow + 4.));
  EXPECT_TRUE(IsFieldEqual(coords->g33, expected_xlow + 3.));
  EXPECT_TRUE(IsFieldEqual(coords->g12, expected_xlow + 2.));
  EXPECT_TRUE(IsFieldEqual(coords->g13, expected_xlow + 1.));
  EXPECT_TRUE(IsFieldEqual(coords->g23, expected_xlow));
}

TEST_F(GridFromOptionsTest, CoordinatesYlowInterp) {
  // *_ylow fields not present in options, Coordinates will be interpolated
  // from CELL_CENTRE

  // make the mesh have boundaries to avoid NaNs in guard cells after interpolating
  mesh_from_options.createBoundaries();

  auto coords = mesh_from_options.getCoordinates(CELL_XLOW);

  Field2D expected_ylow = makeField<Field2D>(
      [](Field2D::ind_type& index) {
        return index.x() + (TWOPI * index.y() - 0.5) + (TWOPI * index.z() / nz) + 3;
      },
      &mesh_from_options);

  mesh_from_options.communicate(expected_ylow);

  EXPECT_TRUE(IsFieldEqual(coords->g11, expected_ylow + 5., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g22, expected_ylow + 4., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g33, expected_ylow + 3., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g12, expected_ylow + 2., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g13, expected_ylow + 1., "RGN_NOBNDRY", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g23, expected_ylow, "RGN_NOBNDRY", this_tolerance));
}

TEST_F(GridFromOptionsTest, CoordinatesYlowRead) {
  // *_ylow fields added to options, will be read to initialise Coordinates

  // Note '(2*pi*11 - y)' here because FakeMesh::GlobalY(int jy) returns jy, not a
  // global position between 0 and 1 (in grid cells, <0 or >1 in boundaries),
  // like a Mesh is supposed to. That means 'y' in input expressions varies
  // between 0 and 2*pi*ny.
  std::string expected_string_ylow{"x + (2*pi*11 - y) + 3"};

  // modify mesh section in global options
  options["dx_ylow"] = "1.";
  options["dy_ylow"] = "1.";
  options["g11_ylow"] = expected_string_ylow + " + 5.";
  options["g22_ylow"] = expected_string_ylow + " + 4.";
  options["g33_ylow"] = expected_string_ylow + " + 3.";
  options["g12_ylow"] = expected_string_ylow + " + 2.";
  options["g13_ylow"] = expected_string_ylow + " + 1.";
  options["g23_ylow"] = expected_string_ylow;

  auto coords = mesh_from_options.getCoordinates(CELL_YLOW);

  Field2D expected_ylow = makeField<Field2D>(
      [](Field2D::ind_type& index) {
        return index.x() + (TWOPI * (ny - index.y())) + (TWOPI * index.z() / nz) + 3;
      },
      &mesh_from_options);

  mesh_from_options.communicate(expected_ylow);

  EXPECT_TRUE(IsFieldEqual(coords->g11, expected_ylow + 5., "RGN_ALL", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g22, expected_ylow + 4., "RGN_ALL", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g33, expected_ylow + 3., "RGN_ALL", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g12, expected_ylow + 2., "RGN_ALL", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g13, expected_ylow + 1., "RGN_ALL", this_tolerance));
  EXPECT_TRUE(IsFieldEqual(coords->g23, expected_ylow, "RGN_ALL", this_tolerance));
}

TEST_F(GridFromOptionsTest, CoordinatesZlowRead) {
  // Grids are axisymmetric, so CELL_ZLOW Coordinates will be read from
  // CELL_CENTRE variables

  auto coords = mesh_from_options.getCoordinates(CELL_ZLOW);

  EXPECT_TRUE(IsFieldEqual(coords->g11, expected_2d + 5.));
  EXPECT_TRUE(IsFieldEqual(coords->g22, expected_2d + 4.));
  EXPECT_TRUE(IsFieldEqual(coords->g33, expected_2d + 3.));
  EXPECT_TRUE(IsFieldEqual(coords->g12, expected_2d + 2.));
  EXPECT_TRUE(IsFieldEqual(coords->g13, expected_2d + 1.));
  EXPECT_TRUE(IsFieldEqual(coords->g23, expected_2d));
}
