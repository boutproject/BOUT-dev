#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "field3d.hxx"
#include "output.hxx"
#include "test_extras.hxx"
#include "unused.hxx"
#include "utils.hxx"

#include <cmath>
#include <set>
#include <vector>

/// Global mesh
extern Mesh *mesh;

/// Test fixture to make sure the global mesh is our fake one
class Field3DTest : public ::testing::Test {
protected:
  static void SetUpTestCase() {
    // Delete any existing mesh
    if (mesh != nullptr) {
      delete mesh;
      mesh = nullptr;
    }
    mesh = new FakeMesh(nx, ny, nz);
    output_info.disable();
    mesh->createDefaultRegions();
    output_info.enable();
  }

  static void TearDownTestCase() {
    delete mesh;
    mesh = nullptr;
  }

public:
  static const int nx;
  static const int ny;
  static const int nz;
};

const int Field3DTest::nx = 3;
const int Field3DTest::ny = 5;
const int Field3DTest::nz = 7;

TEST_F(Field3DTest, IsReal) {
  Field3D field;

  EXPECT_TRUE(field.isReal());
}

TEST_F(Field3DTest, Is3D) {
  Field3D field;

  EXPECT_TRUE(field.is3D());
}

TEST_F(Field3DTest, ByteSize) {
  Field3D field;

  EXPECT_EQ(field.byteSize(), sizeof(BoutReal));
}

TEST_F(Field3DTest, BoutRealSize) {
  Field3D field;

  EXPECT_EQ(field.BoutRealSize(), 1);
}

TEST_F(Field3DTest, Allocate) {
  Field3D field;

  EXPECT_FALSE(field.isAllocated());

  field.allocate();

  EXPECT_TRUE(field.isAllocated());
}

TEST_F(Field3DTest, IsFinite) {
  Field3D field;

  EXPECT_FALSE(finite(field));

  field = 1.0;

  EXPECT_TRUE(finite(field));

  field(1, 1, 1) = std::nan("");

  EXPECT_FALSE(finite(field));
}

TEST_F(Field3DTest, GetGridSizes) {
  Field3D field;

  field.allocate();

  EXPECT_EQ(field.getNx(), nx);
  EXPECT_EQ(field.getNy(), ny);
  EXPECT_EQ(field.getNz(), nz);
}

TEST_F(Field3DTest, CreateOnGivenMesh) {
  int test_nx = Field3DTest::nx + 2;
  int test_ny = Field3DTest::ny + 2;
  int test_nz = Field3DTest::nz + 2;

  FakeMesh *fieldmesh = new FakeMesh(test_nx, test_ny, test_nz);

  Field3D field(fieldmesh);

  field.allocate();

  EXPECT_EQ(field.getNx(), test_nx);
  EXPECT_EQ(field.getNy(), test_ny);
  EXPECT_EQ(field.getNz(), test_nz);

  delete fieldmesh;
}

TEST_F(Field3DTest, CopyCheckFieldmesh) {
  int test_nx = Field3DTest::nx + 2;
  int test_ny = Field3DTest::ny + 2;
  int test_nz = Field3DTest::nz + 2;

  FakeMesh *fieldmesh = new FakeMesh(test_nx, test_ny, test_nz);

  Field3D field(0.0, fieldmesh);

  Field3D field2{field};

  EXPECT_EQ(field2.getNx(), test_nx);
  EXPECT_EQ(field2.getNy(), test_ny);
  EXPECT_EQ(field2.getNz(), test_nz);

  delete fieldmesh;
}

#if CHECK > 0
TEST_F(Field3DTest, CreateOnNullMesh) {
  auto old_mesh = mesh;
  mesh = nullptr;

  Field3D field;

  EXPECT_EQ(field.getNx(), -1);
  EXPECT_EQ(field.getNy(), -1);
  EXPECT_EQ(field.getNz(), -1);

  mesh = old_mesh;

  field.allocate();

  EXPECT_EQ(field.getNx(), Field3DTest::nx);
  EXPECT_EQ(field.getNy(), Field3DTest::ny);
  EXPECT_EQ(field.getNz(), Field3DTest::nz);
}
#endif

#if CHECK > 0 && CHECK <= 2
// We only want to run this test in a certain range of CHECK as we're
// checking some behaviour that is only enabled for CHECK above 0
// but there are checks that will throw before reaching these lines if
// check is greater than 2, so the test only makes sense in a certain range.
TEST_F(Field3DTest, CreateCopyOnNullMesh) {
  // Whilst the declaration of field below looks like it should create a Field2D
  // without a mesh, it in fact will result in a Field2D associated with the
  // global mesh as we end up calling the Field constructor that forces this.
  // Hence, to test the case of copying a field without a mesh we have to
  // temporarily hide the global mesh, before restoring it later.
  auto old_mesh = mesh;
  mesh = nullptr;

  Field3D field;
  // If CHECK > 2 then the following will throw due to the data
  // block in field not being allocated. We can't allocate as that
  // would force field to have a mesh associated with it.
  Field3D field2(field);

  EXPECT_EQ(field2.getNx(), -1);
  EXPECT_EQ(field2.getNy(), -1);
  EXPECT_EQ(field2.getNz(), -1);

  mesh = old_mesh;
  field2.allocate();

  EXPECT_EQ(field2.getNx(), Field3DTest::nx);
  EXPECT_EQ(field2.getNy(), Field3DTest::ny);
  EXPECT_EQ(field2.getNz(), Field3DTest::nz);
}
#endif

TEST_F(Field3DTest, TimeDeriv) {
  Field3D field;

  auto deriv = field.timeDeriv();
  EXPECT_NE(&field, deriv);

  auto deriv2 = field.timeDeriv();
  EXPECT_EQ(deriv, deriv2);

  EXPECT_EQ(&(ddt(field)), deriv);
}

TEST_F(Field3DTest, SplitYupYDown) {
  Field3D field;

  field = 0.;

  EXPECT_FALSE(field.hasYupYdown());

  field.createYupYdown();

  EXPECT_TRUE(field.hasYupYdown());

  auto& yup = field.yup();
  EXPECT_NE(&field, &yup);
  auto& ydown = field.ydown();
  EXPECT_NE(&field, &ydown);

  // Should be able to split again without any problems
  field.createYupYdown();

  // Would be nice to check yup2 != yup, but not sure this is possible
  // to do in general
  auto& yup2 = field.yup();
  EXPECT_NE(&field, &yup2);
  auto& ydown2 = field.ydown();
  EXPECT_NE(&field, &ydown2);
}

TEST_F(Field3DTest, Ynext) {
  Field3D field;

  field = 0.;
  field.createYupYdown();

  auto& yup = field.ynext(1);
  EXPECT_NE(&field, &yup);
  auto& ydown = field.ynext(-1);
  EXPECT_NE(&field, &ydown);
  EXPECT_NE(&yup, &ydown);

  EXPECT_THROW(field.ynext(99), BoutException);
}

TEST_F(Field3DTest, ConstYnext) {
  Field3D field(0.);

  field.createYupYdown();

  const Field3D& field2 = field;

  auto& yup = field2.ynext(1);
  EXPECT_NE(&field2, &yup);
  auto& ydown = field2.ynext(-1);
  EXPECT_NE(&field2, &ydown);
  EXPECT_NE(&yup, &ydown);

  EXPECT_THROW(field2.ynext(99), BoutException);
}

TEST_F(Field3DTest, GetGlobalMesh) {
  Field3D field;

  auto localmesh = field.getMesh();

  EXPECT_EQ(localmesh, mesh);
}

TEST_F(Field3DTest, GetLocalMesh) {
  FakeMesh myMesh{nx + 1, ny + 2, nz + 3};
  Field3D field(&myMesh);

  auto localmesh = field.getMesh();

  EXPECT_EQ(localmesh, &myMesh);
}

TEST_F(Field3DTest, SetGetLocation) {
  Field3D field;

  field.getMesh()->StaggerGrids = true;

  field.setLocation(CELL_XLOW);
  EXPECT_EQ(field.getLocation(), CELL_XLOW);

  field.setLocation(CELL_DEFAULT);
  EXPECT_EQ(field.getLocation(), CELL_CENTRE);

  EXPECT_THROW(field.setLocation(CELL_VSHIFT), BoutException);
}

TEST_F(Field3DTest, SetGetLocationNonStaggered) {
  Field3D field;

  field.getMesh()->StaggerGrids = false;

#if CHECK > 0
  EXPECT_THROW(field.setLocation(CELL_XLOW), BoutException);
  EXPECT_THROW(field.setLocation(CELL_VSHIFT), BoutException);

  field.setLocation(CELL_DEFAULT);
  EXPECT_EQ(field.getLocation(), CELL_CENTRE);
#else
  field.setLocation(CELL_XLOW);
  EXPECT_EQ(field.getLocation(), CELL_CENTRE);

  field.setLocation(CELL_DEFAULT);
  EXPECT_EQ(field.getLocation(), CELL_CENTRE);

  field.setLocation(CELL_VSHIFT);
  EXPECT_EQ(field.getLocation(), CELL_CENTRE);
#endif
}

//-------------------- Iteration tests --------------------

/// This test is split into two parts: a very basic sanity check first
/// (do we visit the right number of elements?), followed by a
/// slightly more complex check one which checks certain indices are
/// actually hit. The latter check works as follows: assign a constant
/// value to the whole field, apart from at certain points, which a
/// different "sentinel" value is used. When we iterate over the
/// field, check those sentinel values are in the correct locations,
/// and the sum over the whole field includes them.
///
/// A more rigorous test would assign a different power of two to each
/// grid point, and check the sum is correct.
TEST_F(Field3DTest, IterateOverWholeField) {
  Field3D field;

  field.allocate();

  // Basic test first: do we visit the correct number of elements?
  int count = 0;
  for (const auto &UNUSED(i) : field) {
    ++count;
  }

  // If this fails, no point doing second part
  ASSERT_EQ(count, nx * ny * nz);

  field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0, 0});
  test_indices.insert({0, 0, 1});
  test_indices.insert({0, 1, 0});
  test_indices.insert({1, 0, 0});
  test_indices.insert({0, 1, 1});
  test_indices.insert({1, 0, 1});
  test_indices.insert({1, 1, 0});
  test_indices.insert({1, 1, 1});
  const int num_sentinels = test_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto& index : test_indices) {
    field(index[0], index[1], index[2]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x, i.y, i.z});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * ny * nz) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(test_indices == result_indices);
}

TEST_F(Field3DTest, IterateOverRGN_ALL) {
  Field3D field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0, 0});
  test_indices.insert({0, 0, 1});
  test_indices.insert({0, 1, 0});
  test_indices.insert({1, 0, 0});
  test_indices.insert({0, 1, 1});
  test_indices.insert({1, 0, 1});
  test_indices.insert({1, 1, 0});
  test_indices.insert({1, 1, 1});
  const int num_sentinels = test_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto index : test_indices) {
    field(index[0], index[1], index[2]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field.region(RGN_ALL)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x, i.y, i.z});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * ny * nz) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(test_indices == result_indices);
}

TEST_F(Field3DTest, IterateOverRegionInd3D_RGN_ALL) {
  Field3D field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect
  std::set<std::vector<int>> test_indices{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0},
                                          {0, 1, 1}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};
  const int num_sentinels = test_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto &index : test_indices) {
    field(index[0], index[1], index[2]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field.getMesh()->getRegion("RGN_ALL")) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y(), i.z()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * ny * nz) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(test_indices == result_indices);
}

TEST_F(Field3DTest, IterateOverRGN_NOBNDRY) {
  Field3D field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0, 0});
  test_indices.insert({0, 0, 1});
  test_indices.insert({0, 1, 0});
  test_indices.insert({1, 0, 0});
  test_indices.insert({0, 1, 1});
  test_indices.insert({1, 0, 1});
  test_indices.insert({1, 1, 0});
  test_indices.insert({1, 1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({1, 1, 0});
  region_indices.insert({1, 1, 1});
  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto &index : test_indices) {
    field(index[0], index[1], index[2]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field.region(RGN_NOBNDRY)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x, i.y, i.z});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum,
            (((nx - 2) * (ny - 2) * nz) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field3DTest, IterateOverRGN_NOX) {
  Field3D field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0, 0});
  test_indices.insert({0, 0, 1});
  test_indices.insert({0, 1, 0});
  test_indices.insert({1, 0, 0});
  test_indices.insert({0, 1, 1});
  test_indices.insert({1, 0, 1});
  test_indices.insert({1, 1, 0});
  test_indices.insert({1, 1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({1, 0, 0});
  region_indices.insert({1, 1, 0});
  region_indices.insert({1, 0, 1});
  region_indices.insert({1, 1, 1});
  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto index : test_indices) {
    field(index[0], index[1], index[2]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field.region(RGN_NOX)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x, i.y, i.z});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, (((nx - 2) * ny * nz) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field3DTest, IterateOverRGN_NOY) {
  Field3D field;

  field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0, 0});
  test_indices.insert({0, 0, 1});
  test_indices.insert({0, 1, 0});
  test_indices.insert({1, 0, 0});
  test_indices.insert({0, 1, 1});
  test_indices.insert({1, 0, 1});
  test_indices.insert({1, 1, 0});
  test_indices.insert({1, 1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({0, 1, 0});
  region_indices.insert({1, 1, 0});
  region_indices.insert({0, 1, 1});
  region_indices.insert({1, 1, 1});
  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto index : test_indices) {
    field(index[0], index[1], index[2]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field.region(RGN_NOY)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x, i.y, i.z});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * (ny - 2) * nz) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field3DTest, IterateOverRGN_NOZ) {
  Field3D field = 1.0;

  // This is not a valid region for Field3D
  EXPECT_THROW(field.region(RGN_NOZ), BoutException);
}

TEST_F(Field3DTest, IterateOver2DRGN_ALL) {
  Field3D field;
  field.allocate();

  for (const auto &i : field) {
    field[i] = 1.0 + i.z;
  }

  BoutReal sum = 0.0;
  for (const auto &i : field.region2D(RGN_ALL)) {
    sum += field[i];
    EXPECT_EQ(field[i], 1.0);
    EXPECT_EQ(i.z, 0);
  }

  EXPECT_EQ(sum, nx * ny);
}

TEST_F(Field3DTest, IterateOver2DRGN_NOBNDRY) {
  Field3D field;
  field.allocate();

  for (const auto &i : field) {
    field[i] = 1.0 + i.z;
  }

  BoutReal sum = 0.0;
  for (const auto &i : field.region2D(RGN_NOBNDRY)) {
    sum += field[i];
    EXPECT_EQ(field[i], 1.0);
    EXPECT_EQ(i.z, 0);
  }

  EXPECT_EQ(sum, nx * ny - 2 * nx - 2 * (ny - 2));
}

TEST_F(Field3DTest, IterateOver2DRGN_NOX) {
  Field3D field;
  field.allocate();

  for (const auto &i : field) {
    field[i] = 1.0 + i.z;
  }

  BoutReal sum = 0.0;
  for (const auto &i : field.region2D(RGN_NOX)) {
    sum += field[i];
    EXPECT_EQ(field[i], 1.0);
    EXPECT_EQ(i.z, 0);
  }

  EXPECT_EQ(sum, nx * ny - 2 * ny);
}

TEST_F(Field3DTest, IterateOver2DRGN_NOY) {
  Field3D field;
  field.allocate();

  for (const auto &i : field) {
    field[i] = 1.0 + i.z;
  }

  BoutReal sum = 0.0;
  for (const auto &i : field.region2D(RGN_NOY)) {
    sum += field[i];
    EXPECT_EQ(field[i], 1.0);
    EXPECT_EQ(i.z, 0);
  }

  EXPECT_EQ(sum, nx * ny - 2 * nx);
}

TEST_F(Field3DTest, IterateOver2DRGN_NOZ) {
  Field3D field = 1.0;
  EXPECT_THROW(field.region2D(RGN_NOZ), BoutException);
}

TEST_F(Field3DTest, Indexing) {
  Field3D field;

  field.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        field(i, j, k) = i + j + k;
      }
    }
  }

  EXPECT_DOUBLE_EQ(field(2, 2, 2), 6);
}

TEST_F(Field3DTest, IndexingInd3D) {
  Field3D field;

  field.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        field(i, j, k) = i + j + k;
      }
    }
  }

  Ind3D ind{(2*ny + 2)*nz + 2};

  EXPECT_DOUBLE_EQ(field[ind], 6);
}

TEST_F(Field3DTest, ConstIndexingInd3D) {
  Field3D field1;

  field1.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        field1(i, j, k) = i + j + k;
      }
    }
  }

  const Field3D field2{field1};

  Ind3D ind{(2*ny + 2)*nz + 2};

  EXPECT_DOUBLE_EQ(field2[ind], 6);
}

TEST_F(Field3DTest, IndexingInd2D) {
  Field3D field(0.0);
  const BoutReal sentinel = 2.0;
  int ix = 1, iy = 2, iz = 3;
  field(ix, iy, iz) = sentinel;

  Ind2D ind{iy + ny * ix, ny, 1};
  EXPECT_DOUBLE_EQ(field(ind, iz), sentinel);
  field(ind, iz) = -sentinel;
  EXPECT_DOUBLE_EQ(field(ix, iy, iz), -sentinel);
}

TEST_F(Field3DTest, IndexingIndPerp) {
  Field3D field(0.0);
  const BoutReal sentinel = 2.0;
  int ix = 1, iy = 2, iz = 3;
  field(ix, iy, iz) = sentinel;

  IndPerp ind{iz + nz * ix, 1, nz};
  EXPECT_DOUBLE_EQ(field(ind, iy), sentinel);
  field(ind, iy) = -sentinel;
  EXPECT_DOUBLE_EQ(field(ix, iy, iz), -sentinel);
}

TEST_F(Field3DTest, IndexingToZPointer) {
  Field3D field;

  field.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        field(i, j, k) = i + j + k;
      }
    }
  }

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      auto tmp = field(i, j);
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(tmp[k], i + j + k);
        tmp[k] = -1.0;
      }
    }
  }

  for (const auto &i : field) {
    EXPECT_EQ(field[i], -1.0);
  }

#if CHECK > 2
  EXPECT_THROW(field(-1, 0), BoutException);
  EXPECT_THROW(field(0, -1), BoutException);
  EXPECT_THROW(field(nx, 0), BoutException);
  EXPECT_THROW(field(0, ny), BoutException);

  Field3D fieldUnassigned;
  EXPECT_THROW(fieldUnassigned(0, 0), BoutException);
#endif
}

TEST_F(Field3DTest, ConstIndexingToZPointer) {
  const Field3D field = 1.0;
  Field3D field2 = 0.0;

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      auto tmp = field(i, j);
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(tmp[k], 1.0);
        field2(i, j, k) = tmp[k];
      }
    }
  }

  for (const auto &i : field2) {
    EXPECT_EQ(field2[i], 1.0);
  }

#if CHECK > 2
  EXPECT_THROW(field(-1, 0), BoutException);
  EXPECT_THROW(field(0, -1), BoutException);
  EXPECT_THROW(field(nx, 0), BoutException);
  EXPECT_THROW(field(0, ny), BoutException);

  const Field3D fieldUnassigned;
  EXPECT_THROW(fieldUnassigned(0, 0), BoutException);

#endif
}

//-------------------- Checking tests --------------------

#if CHECK > 2
TEST_F(Field3DTest, CheckNotEmpty) {
  Field3D field;

  EXPECT_THROW(field(0, 0, 0), BoutException);
}

TEST_F(Field3DTest, BoundsCheck) {
  Field3D field;

  field = 1.0;
  EXPECT_THROW(field(-1, 0, 0), BoutException);
  EXPECT_THROW(field(0, -1, 0), BoutException);
  EXPECT_THROW(field(0, 0, -1), BoutException);
  EXPECT_THROW(field(nx, 0, 0), BoutException);
  EXPECT_THROW(field(0, ny, 0), BoutException);
  EXPECT_THROW(field(0, 0, nz), BoutException);
}

TEST_F(Field3DTest, ConstCheckNotEmpty) {
  const Field3D field;

  EXPECT_THROW(field(0, 0, 0), BoutException);
}

TEST_F(Field3DTest, ConstBoundsCheck) {
  const Field3D field = 1.0;

  EXPECT_THROW(field(-1, 0, 0), BoutException);
  EXPECT_THROW(field(0, -1, 0), BoutException);
  EXPECT_THROW(field(0, 0, -1), BoutException);
  EXPECT_THROW(field(nx, 0, 0), BoutException);
  EXPECT_THROW(field(0, ny, 0), BoutException);
  EXPECT_THROW(field(0, 0, nz), BoutException);
}

TEST_F(Field3DTest, CheckData) {
  Field3D field;

  EXPECT_THROW(checkData(field), BoutException);

  field = 1.0;

  EXPECT_NO_THROW(checkData(field));

  field(1, 1, 1) = std::nan("");

  EXPECT_THROW(checkData(field), BoutException);

  field = 1.0;
  field(0, 0, 0) = std::nan("");

  EXPECT_NO_THROW(checkData(field));
  EXPECT_NO_THROW(checkData(field, RGN_NOBNDRY));
  EXPECT_THROW(checkData(field, RGN_ALL), BoutException);
  
}

#if CHECK > 0
TEST_F(Field3DTest, BndryValid) {
  Field3D field = 1.0;
  field.bndry_xin = true;
  field.bndry_xout = true;
  field.bndry_yup = true;
  field.bndry_ydown = true;
  EXPECT_EQ(field.bndryValid(), true);

  field.bndry_xin = false;
  EXPECT_THROW(field.bndryValid(), BoutException);
  field.bndry_xin = true;

  field.bndry_xout = false;
  EXPECT_THROW(field.bndryValid(), BoutException);
  field.bndry_xout = true;

  field.bndry_yup = false;
  EXPECT_THROW(field.bndryValid(), BoutException);
  field.bndry_yup = true;

  field.bndry_ydown = false;
  EXPECT_THROW(field.bndryValid(), BoutException);
  field.bndry_ydown = true;
}

TEST_F(Field3DTest, DoneComms) {
  Field3D field = 1.0;
  field.bndry_xin = false;
  field.bndry_xout = false;
  field.bndry_yup = false;
  field.bndry_ydown = false;

  EXPECT_THROW(field.bndryValid(), BoutException);
  field.doneComms();
  EXPECT_EQ(field.bndryValid(), true);
}
#endif

TEST_F(Field3DTest, InvalidateGuards) {
  Field3D field;
  field.allocate(); // Calls invalidateGuards
  field = 1.0;      // Sets everywhere including boundaries

  const int nmesh = nx * ny * nz;

  int sum = 0;
  for (const auto &i : field.region(RGN_ALL)) {
    field[i] = 0.0; // Reset field value
    sum++;
  }
  EXPECT_EQ(sum, nmesh); // Field operator= hasn't been broken by invalidateGuards

  // Count the number of non-boundary points
  sum = 0;
  for (const auto &i : field.region(RGN_NOBNDRY)) {
    field[i] = 0.0; // Reset field value
    sum++;
  }
  const int nbndry = nmesh - sum;

  auto localmesh = field.getMesh();
  EXPECT_NO_THROW(checkData(field(0, 0, 0)));
  EXPECT_NO_THROW(checkData(field(localmesh->xstart, localmesh->ystart, 0)));

  invalidateGuards(field);

  EXPECT_THROW(checkData(field(0, 0, 0)), BoutException);
  EXPECT_NO_THROW(checkData(field(localmesh->xstart, localmesh->ystart, 0)));

  sum = 0;
  for (const auto &i : field.region(RGN_ALL)) {
    if (!finite(field[i]))
      sum++;
  }
  EXPECT_EQ(sum, nbndry);
}

#endif // CHECK > 2

//-------------------- Assignment tests --------------------

TEST_F(Field3DTest, CreateFromBoutReal) {
  Field3D field(1.0);

  EXPECT_TRUE(IsField3DEqualBoutReal(field, 1.0));
}

TEST_F(Field3DTest, CreateFromField3D) {
  Field3D field(99.0);
  Field3D result(field);

  EXPECT_TRUE(IsField3DEqualBoutReal(result, 99.0));
}

TEST_F(Field3DTest, CreateFromField2D) {
  Field2D field(99.0);
  Field3D result(field);

  EXPECT_TRUE(IsField3DEqualBoutReal(result, 99.0));
}

TEST_F(Field3DTest, AssignFromBoutReal) {
  Field3D field;

  field = 2.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(field, 2.0));
}

TEST_F(Field3DTest, AssignFromInvalid) {
  Field3D field;

#if CHECK > 0
  EXPECT_THROW(field = std::nan(""), BoutException);
#else
  EXPECT_NO_THROW(field = std::nan(""));
#endif
}

TEST_F(Field3DTest, AssignFromField2D) {
  Field3D field;
  Field2D field2(2.0);

  field = field2;

  EXPECT_TRUE(IsField3DEqualBoutReal(field, 2.0));

#if CHECK > 0
  Field2D field3;
  EXPECT_THROW(field = field3, BoutException);
#endif
}

TEST_F(Field3DTest, AssignFromFieldPerp) {
  Field3D field = 1.0;
  // Note we have to pass the mesh to the FieldPerp constructor
  // unlike other fields.
  FieldPerp field2(mesh);
  const int yindex = 2;
  field2.setIndex(yindex);
  field2 = 2.0;
  field = field2;

  for (const auto &i : field) {
    if (i.y == yindex) {
      EXPECT_EQ(field[i], 2.0);
    } else {
      EXPECT_EQ(field[i], 1.0);
    }
  }

#if CHECK > 0
  FieldPerp field3;
  EXPECT_THROW(field = field3, BoutException);
  FieldPerp field4(mesh);
  EXPECT_THROW(field = field4, BoutException);
#endif
}

TEST_F(Field3DTest, AssignFromField3D) {
  Field3D field, field2;

  field2 = -99.0;
  field = field2;

  EXPECT_TRUE(IsField3DEqualBoutReal(field, -99.0));

#if CHECK > 0
  Field3D field3;
  EXPECT_THROW(field = field3, BoutException);
#endif
}

//-------------------- Arithmetic tests --------------------

TEST_F(Field3DTest, UnaryMinus) {
  Field3D field;

  field = 2.0;
  field = -field;

  EXPECT_TRUE(IsField3DEqualBoutReal(field, -2.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(-field, 2.0));
}

TEST_F(Field3DTest, AddEqualsBoutReal) {
  Field3D a;

  a = 1.0;
  a += 5.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 6.0));

  // Check case where field is not unique
  auto c = a;
  c += 5.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 6.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, 11.0));
}

TEST_F(Field3DTest, AddEqualsField2D) {
  Field3D a;
  Field2D b;

  a = 2.0;
  b = 3.0;
  a += b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 5.0));

  // Check case where field is not unique
  auto c = a;
  c += b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 5.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, 8.0));
}

TEST_F(Field3DTest, AddEqualsField3D) {
  Field3D a, b;

  a = 2.0;
  b = 3.0;
  a += b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 5.0));

  // Check case where field is not unique
  auto c = a;
  c += b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 5.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, 8.0));
}

TEST_F(Field3DTest, AddEqualsField3DField3DStagger) {
  auto backup = mesh->StaggerGrids;
  mesh->StaggerGrids = true; // Force staggering

  Field3D a, b;

  a = 2.0;
  b = 3.0;

  a.setLocation(CELL_XLOW);
  b.setLocation(CELL_CENTRE);

// Throw as two rhs fields at different locations
#if CHECK > 0
  EXPECT_THROW(a += b, BoutException);
#else
  EXPECT_NO_THROW(a += b);
#endif

  mesh->StaggerGrids = backup;
}

TEST_F(Field3DTest, AddField3DBoutReal) {
  Field3D a, b;

  a = 1.0;
  b = a + 2.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(b, 3.0));
}

TEST_F(Field3DTest, AddBoutRealField3D) {
  Field3D a, b;

  a = 1.0;
  b = 3.0 + a;

  EXPECT_TRUE(IsField3DEqualBoutReal(b, 4.0));
}

TEST_F(Field3DTest, AddField2DField3D) {
  Field2D a;
  Field3D b, c;

  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 3.0));
}

TEST_F(Field3DTest, AddField3DField2D) {
  Field3D a, c;
  Field2D b;

  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 3.0));
}

TEST_F(Field3DTest, AddField3DField3D) {
  Field3D a, b, c;

  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 3.0));
}

TEST_F(Field3DTest, AddField3DField3DStagger) {
  auto backup = mesh->StaggerGrids;
  mesh->StaggerGrids = true; // Force staggering

  Field3D a, b, c;

  a = 1.0;
  b = 2.0;
  c = 3.0;

  a.setLocation(CELL_XLOW);
  b.setLocation(CELL_CENTRE);
  c.setLocation(CELL_CENTRE);

  // Throw as two rhs fields at different locations
#if CHECK > 0
  EXPECT_THROW(c = a + b, BoutException);
#else
  EXPECT_NO_THROW(c = a + b);
#endif

  // No throw as updates location of a
  EXPECT_NO_THROW(a = c + b);

  // Hence the first case should now not throw
  EXPECT_NO_THROW(c = a + b);

  mesh->StaggerGrids = backup;
}

TEST_F(Field3DTest, MultiplyEqualsBoutReal) {
  Field3D a;

  a = 2.0;
  a *= 1.5;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 3.0));

  // Check case where field is not unique
  auto c = a;
  c *= 1.5;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 3.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, 4.5));
}

TEST_F(Field3DTest, MultiplyEqualsField2D) {
  Field3D a;
  Field2D b;

  a = 2.5;
  b = 4.0;
  a *= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 10.0));

  // Check case where field is not unique
  auto c = a;
  c *= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 10.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, 40.0));
}

TEST_F(Field3DTest, MultiplyEqualsField3D) {
  Field3D a, b;

  a = 2.5;
  b = 4.0;
  a *= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 10.0));

  // Check case where field is not unique
  auto c = a;
  c *= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 10.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, 40.0));
}

TEST_F(Field3DTest, MultiplyEqualsField3DField3DStagger) {
  auto backup = mesh->StaggerGrids;
  mesh->StaggerGrids = true; // Force staggering

  Field3D a, b;

  a = 2.5;
  b = 4.0;

  a.setLocation(CELL_XLOW);
  b.setLocation(CELL_CENTRE);

// Throw as two rhs fields at different locations
#if CHECK > 0
  EXPECT_THROW(a *= b, BoutException);
#else
  EXPECT_NO_THROW(a *= b);
#endif

  mesh->StaggerGrids = backup;
}

TEST_F(Field3DTest, MultiplyField3DBoutReal) {
  Field3D a, b;

  a = 1.5;
  b = a * 2.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(b, 3.0));
}

TEST_F(Field3DTest, MultiplyBoutRealField3D) {
  Field3D a, b;

  a = 2.5;
  b = 3.0 * a;

  EXPECT_TRUE(IsField3DEqualBoutReal(b, 7.5));
}

TEST_F(Field3DTest, MultiplyField2DField3D) {
  Field2D a;
  Field3D b, c;

  a = 4.0;
  b = 4.0;
  c = a * b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 16.0));
}

TEST_F(Field3DTest, MultiplyField3DField2D) {
  Field3D a, c;
  Field2D b;

  a = 8.0;
  b = 8.0;
  c = a * b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 64.0));
}

TEST_F(Field3DTest, MultiplyField3DField3D) {
  Field3D a, b, c;

  a = 4.0;
  b = 8.0;
  c = a * b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 32.0));
}

TEST_F(Field3DTest, MultiplyField3DField3DStagger) {
  auto backup = mesh->StaggerGrids;
  mesh->StaggerGrids = true; // Force staggering

  Field3D a, b, c;

  a = 1.0;
  b = 2.0;
  c = 3.0;

  a.setLocation(CELL_XLOW);
  b.setLocation(CELL_CENTRE);
  c.setLocation(CELL_CENTRE);

  // Throw as two rhs fields at different locations
#if CHECK > 0
  EXPECT_THROW(c = a * b, BoutException);
#else
  EXPECT_NO_THROW(c = a * b);
#endif

  // No throw as updates location of a
  EXPECT_NO_THROW(a = c * b);

  // Hence the first case should now not throw
  EXPECT_NO_THROW(c = a * b);

  mesh->StaggerGrids = backup;
}

TEST_F(Field3DTest, SubtractEqualsBoutReal) {
  Field3D a;

  a = 1.0;
  a -= 5.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, -4.0));

  // Check case where field is not unique
  auto c = a;
  c -= 5.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, -4.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, -9.0));
}

TEST_F(Field3DTest, SubtractEqualsField2D) {
  Field3D a;
  Field2D b;

  a = 2.0;
  b = 7.0;
  a -= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, -5.0));

  // Check case where field is not unique
  auto c = a;
  c -= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, -5.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, -12.0));
}

TEST_F(Field3DTest, SubtractEqualsField3D) {
  Field3D a, b;

  a = 2.0;
  b = 7.0;
  a -= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, -5.0));

  // Check case where field is not unique
  auto c = a;
  c -= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, -5.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, -12.0));
}

TEST_F(Field3DTest, SubtractEqualsField3DField3DStagger) {
  auto backup = mesh->StaggerGrids;
  mesh->StaggerGrids = true; // Force staggering

  Field3D a, b;

  a = 2.0;
  b = 7.0;

  a.setLocation(CELL_XLOW);
  b.setLocation(CELL_CENTRE);

// Throw as two rhs fields at different locations
#if CHECK > 0
  EXPECT_THROW(a -= b, BoutException);
#else
  EXPECT_NO_THROW(a -= b);
#endif

  mesh->StaggerGrids = backup;
}

TEST_F(Field3DTest, SubtractField3DBoutReal) {
  Field3D a, b;

  a = 10.0;
  b = a - 2.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(b, 8.0));
}

TEST_F(Field3DTest, SubtractBoutRealField3D) {
  Field3D a, b;

  a = 10.0;
  b = 3.0 - a;

  EXPECT_TRUE(IsField3DEqualBoutReal(b, -7.0));
}

TEST_F(Field3DTest, SubtractField2DField3D) {
  Field2D a;
  Field3D b, c;

  a = 10.0;
  b = 20.0;
  c = a - b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, -10.0));
}

TEST_F(Field3DTest, SubtractField3DField2D) {
  Field3D a, c;
  Field2D b;

  a = 10.0;
  b = 20.0;
  c = a - b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, -10.0));
}

TEST_F(Field3DTest, SubtractField3DField3D) {
  Field3D a, b, c;

  a = 10.0;
  b = 20.0;
  c = a - b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, -10.0));
}

TEST_F(Field3DTest, SubtractField3DField3DStagger) {
  auto backup = mesh->StaggerGrids;
  mesh->StaggerGrids = true; // Force staggering

  Field3D a, b, c;

  a = 1.0;
  b = 2.0;
  c = 3.0;

  a.setLocation(CELL_XLOW);
  b.setLocation(CELL_CENTRE);
  c.setLocation(CELL_CENTRE);

  // Throw as two rhs fields at different locations
#if CHECK > 0
  EXPECT_THROW(c = a - b, BoutException);
#else
  EXPECT_NO_THROW(c = a - b);
#endif

  // No throw as updates location of a
  EXPECT_NO_THROW(a = c - b);

  // Hence the first case should now not throw
  EXPECT_NO_THROW(c = a - b);

  mesh->StaggerGrids = backup;
}

TEST_F(Field3DTest, DivideEqualsBoutReal) {
  Field3D a;

  a = 2.5;
  a /= 5.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 0.5));

  // Check case where field is not unique
  auto c = a;
  c /= 5.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 0.5));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, 0.1));
}

TEST_F(Field3DTest, DivideEqualsField2D) {
  Field3D a;
  Field2D b;

  a = 5.0;
  b = 2.5;
  a /= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 2.0));

  // Check case where field is not unique
  auto c = a;
  c /= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 2.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, 0.8));
}

TEST_F(Field3DTest, DivideEqualsField3D) {
  Field3D a, b;

  a = 5.0;
  b = 2.5;
  a /= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 2.0));

  // Check case where field is not unique
  auto c = a;
  c /= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 2.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(c, 0.8));
}

TEST_F(Field3DTest, DivideEqualsField3DField3DStagger) {
  auto backup = mesh->StaggerGrids;
  mesh->StaggerGrids = true; // Force staggering

  Field3D a, b;

  a = 5.0;
  b = 2.5;

  a.setLocation(CELL_XLOW);
  b.setLocation(CELL_CENTRE);

// Throw as two rhs fields at different locations
#if CHECK > 0
  EXPECT_THROW(a /= b, BoutException);
#else
  EXPECT_NO_THROW(a /= b);
#endif

  mesh->StaggerGrids = backup;
}

TEST_F(Field3DTest, DivideField3DBoutReal) {
  Field3D a, b;

  a = 3.0;
  b = a / 2.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(b, 1.5));
}

TEST_F(Field3DTest, DivideBoutRealField3D) {
  Field3D a, b;

  a = 2.5;
  b = 10.0 / a;

  EXPECT_TRUE(IsField3DEqualBoutReal(b, 4.0));
}

TEST_F(Field3DTest, DivideField2DField3D) {
  Field2D a;
  Field3D b, c;

  a = 32.0;
  b = 8.0;
  c = a / b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 4.0));
}

TEST_F(Field3DTest, DivideField3DField2D) {
  Field3D a, c;
  Field2D b;

  a = 32.0;
  b = 8.0;
  c = a / b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 4.0));
}

TEST_F(Field3DTest, DivideField3DField3D) {
  Field3D a, b, c;

  a = 32.0;
  b = 8.0;
  c = a / b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 4.0));
}

TEST_F(Field3DTest, DivideField3DField3DStagger) {
  auto backup = mesh->StaggerGrids;
  mesh->StaggerGrids = true; // Force staggering

  Field3D a, b, c;

  a = 1.0;
  b = 2.0;
  c = 3.0;

  a.setLocation(CELL_XLOW);
  b.setLocation(CELL_CENTRE);
  c.setLocation(CELL_CENTRE);

  // Throw as two rhs fields at different locations
#if CHECK > 0
  EXPECT_THROW(c = a / b, BoutException);
#else
  EXPECT_NO_THROW(c = a / b);
#endif

  // No throw as updates location of a
  EXPECT_NO_THROW(a = c / b);

  // Hence the first case should now not throw
  EXPECT_NO_THROW(c = a / b);

  mesh->StaggerGrids = backup;
}

TEST_F(Field3DTest, PowBoutRealField3D) {
  Field3D a, b;
  a = 5.0;
  b = pow(2.0, a);

  EXPECT_TRUE(IsField3DEqualBoutReal(b, 32.0));
}

TEST_F(Field3DTest, PowField3DBoutReal) {
  Field3D a, b;
  a = 5.0;
  b = pow(a, 2.0);

  EXPECT_TRUE(IsField3DEqualBoutReal(b, 25.0));
}

TEST_F(Field3DTest, PowField3DFieldPerp) {
  Field3D a{mesh};
  FieldPerp b{mesh}, c{mesh};
  const int yindex = 2;
  b.setIndex(yindex);

  a = 2.0;
  b = 6.0;
  c = pow(a, b);

  EXPECT_TRUE(IsFieldPerpEqualBoutReal(c, 64.0));
}

TEST_F(Field3DTest, PowField3DField2D) {
  Field3D a, c;
  Field2D b;

  a = 2.0;
  b = 6.0;
  c = pow(a, b);

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 64.0));
}

TEST_F(Field3DTest, PowField3DField3D) {
  Field3D a, b, c;
  a = 2.0;
  b = 6.0;
  c = pow(a, b);

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 64.0));
}

TEST_F(Field3DTest, Sqrt) {
  Field3D field;

  field = 16.0;
  EXPECT_TRUE(IsField3DEqualBoutReal(sqrt(field), 4.0));
}

TEST_F(Field3DTest, Abs) {
  Field3D field;

  field = -31.0;
  EXPECT_TRUE(IsField3DEqualBoutReal(abs(field), 31.0));
}

TEST_F(Field3DTest, Exp) {
  Field3D field;

  field = 2.5;
  const BoutReal expected = 12.182493960703473;
  EXPECT_TRUE(IsField3DEqualBoutReal(exp(field), expected));
}

TEST_F(Field3DTest, Log) {
  Field3D field;

  field = 12.182493960703473;
  const BoutReal expected = 2.5;
  EXPECT_TRUE(IsField3DEqualBoutReal(log(field), expected));
}

TEST_F(Field3DTest, LogExp) {
  Field3D field;

  field = 2.5;
  const BoutReal expected = 2.5;
  EXPECT_TRUE(IsField3DEqualBoutReal(log(exp(field)), expected));
}

TEST_F(Field3DTest, Sin) {
  Field3D field;

  field = PI / 2.0;
  EXPECT_TRUE(IsField3DEqualBoutReal(sin(field), 1.0));

  field = PI;
  EXPECT_TRUE(IsField3DEqualBoutReal(sin(field), 0.0));
}

TEST_F(Field3DTest, Cos) {
  Field3D field;

  field = PI / 2.0;
  EXPECT_TRUE(IsField3DEqualBoutReal(cos(field), 0.0));

  field = PI;
  EXPECT_TRUE(IsField3DEqualBoutReal(cos(field), -1.0));
}

TEST_F(Field3DTest, Tan) {
  Field3D field;

  field = PI / 4.0;
  EXPECT_TRUE(IsField3DEqualBoutReal(tan(field), 1.0));

  field = PI;
  EXPECT_TRUE(IsField3DEqualBoutReal(tan(field), 0.0));
}

TEST_F(Field3DTest, Sinh) {
  Field3D field;

  field = 1.0;
  const BoutReal expected = 1.1752011936438014;
  EXPECT_TRUE(IsField3DEqualBoutReal(sinh(field), expected));

  field = -1.0;
  EXPECT_TRUE(IsField3DEqualBoutReal(sinh(field), -expected));
}

TEST_F(Field3DTest, Cosh) {
  Field3D field;

  field = 1.0;
  const BoutReal expected = 1.5430806348152437;
  EXPECT_TRUE(IsField3DEqualBoutReal(cosh(field), expected));

  field = -1.0;
  EXPECT_TRUE(IsField3DEqualBoutReal(cosh(field), expected));
}

TEST_F(Field3DTest, Tanh) {
  Field3D field;

  field = 1.0;
  const BoutReal expected = 0.761594155955764;
  EXPECT_TRUE(IsField3DEqualBoutReal(tanh(field), expected));

  field = -1.0;
  EXPECT_TRUE(IsField3DEqualBoutReal(tanh(field), -expected));
}

TEST_F(Field3DTest, Floor) {
  Field3D field;

  field = 50.0;
  field(1, 1, 1) = 49.9;
  field(2, 3, 4) = -20;

  const BoutReal floor_value = 50.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(floor(field, floor_value), floor_value));
}

TEST_F(Field3DTest, Min) {
  Field3D field;

  field = 50.0;
  field(0, 0, 0) = -99.0;
  field(1, 1, 1) = 60.0;
  field(1, 2, 2) = 40.0;
  field(2, 4, 3) = 99.0;

  // min doesn't include guard cells
  const BoutReal min_value = 40.0;

  EXPECT_EQ(min(field, false), min_value);
  EXPECT_EQ(min(field, false, RGN_ALL), -99.0);
  EXPECT_EQ(min(field, true, RGN_ALL), -99.0);
}

TEST_F(Field3DTest, Max) {
  Field3D field;

  field = 50.0;
  field(0, 0, 0) = -99.0;
  field(1, 1, 1) = 40.0;
  field(1, 2, 2) = 60.0;
  field(2, 4, 3) = 99.0;

  // max doesn't include guard cells
  const BoutReal max_value = 60.0;

  EXPECT_EQ(max(field, false), max_value);
  EXPECT_EQ(max(field, false, RGN_ALL), 99.0);
  EXPECT_EQ(max(field, true, RGN_ALL), 99.0);
}

TEST_F(Field3DTest, Mean) {
  Field3D field;

  field = 50.0;
  field(0, 0, 0) = 1.0;
  field(1, 1, 1) = 40.0;
  field(1, 2, 2) = 60.0;
  field(2, 4, 3) = 109.0;

  // mean doesn't include guard cells by default
  const int npoints_all = nx*ny*nz;
  const BoutReal mean_value_nobndry = 50.0;
  const BoutReal mean_value_all = 50.0 + 10.0/npoints_all;

  EXPECT_EQ(mean(field, false), mean_value_nobndry);
  EXPECT_EQ(mean(field, false, RGN_ALL), mean_value_all);
  EXPECT_EQ(mean(field, true, RGN_ALL), mean_value_all);
}

TEST_F(Field3DTest, DC) {
  Field3D field;

  field = 1.0;
  for (const auto& i : field) {
    field[i] = i.z;
  }

  EXPECT_TRUE(IsField2DEqualBoutReal(DC(field), 3.0));
}

#ifdef _OPENMP
// This test may be more of a DataIterator test so should perhaps
// be migrated to a separate test file.
TEST_F(Field3DTest, OpenMPIterator) {
  const int fields = 10;
  Field3D *d3 = new Field3D[fields];
  for (int i = 0; i < fields; ++i) {
    d3[i] = 1;
  }

  BoutReal expected = 1;

  BOUT_OMP(parallel for)
  for (int j = 0; j < fields; ++j) {
    BOUT_OMP(parallel)
    for (auto i : d3[j]) {
      EXPECT_DOUBLE_EQ(d3[j][i], expected);
      d3[j][i] = 2;
    }
  }

  expected = 2;

  for (int i = 0; i < fields; ++i) {
    for (int jx = 0; jx < mesh->LocalNx; ++jx) {
      for (int jy = 0; jy < mesh->LocalNy; ++jy) {
        for (int jz = 0; jz < mesh->LocalNz; ++jz) {
          EXPECT_DOUBLE_EQ(d3[i](jx, jy, jz), expected);
        }
      }
    }
  }

  delete[] d3;
}
#endif
