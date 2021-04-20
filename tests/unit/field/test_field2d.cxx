// We know stuff might be deprecated, but we still want to test it
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "field2d.hxx"
#include "output.hxx"
#include "test_extras.hxx"
#include "unused.hxx"
#include "utils.hxx"

#include <cmath>
#include <set>
#include <vector>

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using Field2DTest = FakeMeshFixture;

TEST_F(Field2DTest, IsReal) {
  Field2D field;

  EXPECT_TRUE(field.isReal());
}

TEST_F(Field2DTest, Is3D) {
  Field2D field;

  EXPECT_FALSE(field.is3D());
}

TEST_F(Field2DTest, ByteSize) {
  Field2D field;

  EXPECT_EQ(field.byteSize(), sizeof(BoutReal));
}

TEST_F(Field2DTest, BoutRealSize) {
  Field2D field;

  EXPECT_EQ(field.BoutRealSize(), 1);
}

TEST_F(Field2DTest, Allocate) {
  Field2D field;

  EXPECT_FALSE(field.isAllocated());

  field.allocate();

  EXPECT_TRUE(field.isAllocated());
}

TEST_F(Field2DTest, IsFinite) {
  Field2D field;

  EXPECT_FALSE(finite(field));

  field = 1.0;

  EXPECT_TRUE(finite(field));

  field(1, 1) = std::nan("");

  EXPECT_FALSE(finite(field));
}

TEST_F(Field2DTest, GetGridSizes) {
  Field2D field;

  EXPECT_EQ(field.getNx(), nx);
  EXPECT_EQ(field.getNy(), ny);
  EXPECT_EQ(field.getNz(), 1);
}

TEST_F(Field2DTest, CreateOnGivenMesh) {
  int test_nx = Field2DTest::nx + 2;
  int test_ny = Field2DTest::ny + 2;
  int test_nz = Field2DTest::nz + 2;

  FakeMesh fieldmesh{test_nx, test_ny, test_nz};
  fieldmesh.setCoordinates(nullptr);

  Field2D field{&fieldmesh};

  EXPECT_EQ(field.getNx(), test_nx);
  EXPECT_EQ(field.getNy(), test_ny);
  EXPECT_EQ(field.getNz(), 1);
}

TEST_F(Field2DTest, CopyCheckFieldmesh) {
  int test_nx = Field2DTest::nx + 2;
  int test_ny = Field2DTest::ny + 2;
  int test_nz = Field2DTest::nz + 2;

  FakeMesh fieldmesh{test_nx, test_ny, test_nz};
  fieldmesh.setCoordinates(nullptr);

  // createDefaultRegions is noisy
  WithQuietOutput quiet{output_info};
  fieldmesh.createDefaultRegions();

  Field2D field{1.0, &fieldmesh};
  Field2D field2{field};

  EXPECT_EQ(field2.getNx(), test_nx);
  EXPECT_EQ(field2.getNy(), test_ny);
  EXPECT_EQ(field2.getNz(), 1);
  EXPECT_TRUE(areFieldsCompatible(field, field2));
}

#if CHECK > 0
TEST_F(Field2DTest, CreateOnNullMesh) {
  auto old_mesh = mesh;
  mesh = nullptr;

  Field2D field;

  EXPECT_EQ(field.getNx(), -1);
  EXPECT_EQ(field.getNy(), -1);
  EXPECT_EQ(field.getNz(), 1);

  mesh = old_mesh;

  field.allocate();

  EXPECT_EQ(field.getNx(), Field2DTest::nx);
  EXPECT_EQ(field.getNy(), Field2DTest::ny);
  EXPECT_EQ(field.getNz(), 1);
}
#endif

#if CHECK > 0 && CHECK <= 2
// We only want to run this test in a certain range of CHECK as we're
// checking some behaviour that is only enabled for CHECK above 0
// but there are checks that will throw before reaching these lines if
// check is greater than 2, so the test only makes sense in a certain range.
TEST_F(Field2DTest, CreateCopyOnNullMesh) {
  // Whilst the declaration of field below looks like it should create a Field2D
  // without a mesh, it in fact will result in a Field2D associated with the
  // global mesh as we end up calling the Field constructor that forces this.
  // Hence, to test the case of copying a field without a mesh we have to
  // temporarily hide the global mesh, before restoring it later.
  auto old_mesh = mesh;
  mesh = nullptr;

  Field2D field;
  // If CHECK > 2 then the following will throw due to the data
  // block in field not being allocated. We can't allocate as that
  // would force field to have a mesh associated with it.
  Field2D field2(field);

  EXPECT_EQ(field2.getNx(), -1);
  EXPECT_EQ(field2.getNy(), -1);
  EXPECT_EQ(field2.getNz(), 1);

  mesh = old_mesh;
  field2.allocate();

  EXPECT_EQ(field2.getNx(), Field2DTest::nx);
  EXPECT_EQ(field2.getNy(), Field2DTest::ny);
  EXPECT_EQ(field2.getNz(), 1);
}
#endif

TEST_F(Field2DTest, TimeDeriv) {
  Field2D field;

  auto deriv = field.timeDeriv();
  EXPECT_NE(&field, deriv);

  auto deriv2 = field.timeDeriv();
  EXPECT_EQ(deriv, deriv2);

  EXPECT_EQ(&(ddt(field)), deriv);
}

TEST_F(Field2DTest, SetGetLocation) {
  Field2D field(mesh_staggered);

  field.setLocation(CELL_XLOW);
  EXPECT_EQ(field.getLocation(), CELL_XLOW);

  field.setLocation(CELL_DEFAULT);
  EXPECT_EQ(field.getLocation(), CELL_CENTRE);

  EXPECT_THROW(field.setLocation(CELL_VSHIFT), BoutException);
}

TEST_F(Field2DTest, SetGetLocationNonStaggered) {
  Field2D field;

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
TEST_F(Field2DTest, IterateOverWholeField) {
  Field2D field;

  field.allocate();

  // Basic test first: do we visit the correct number of elements?
  int count = 0;
  for (auto &UNUSED(i) : field) {
    ++count;
  }

  // If this fails, no point doing second part
  ASSERT_EQ(count, nx * ny);

  field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0});
  test_indices.insert({0, 1});
  test_indices.insert({1, 0});
  test_indices.insert({1, 1});
  const int num_sentinels = test_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (auto index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (auto &i : field) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * ny) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(test_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRegionInd2D_RGN_ALL) {
  Field2D field{1.0};

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect
  std::set<std::vector<int>> test_indices{{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  const int num_sentinels = test_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto &index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (auto &i : field.getRegion("RGN_ALL")) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * ny) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(test_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRegionInd3D_RGN_ALL) {
  Field2D field{1.0};

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect
  std::set<std::vector<int>> test_indices{{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  // nz more sentinels than expected, as we're looping over a 3D region
  const int num_sentinels = test_indices.size() * nz;

  // Assign sentinel value to watch out for to our chosen points
  for (const auto &index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (auto &i : field.getMesh()->getRegion3D("RGN_ALL")) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nz * nx * ny) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(test_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRGN_NOBNDRY) {
  Field2D field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0});
  test_indices.insert({0, 1});
  test_indices.insert({1, 0});
  test_indices.insert({1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({1, 1});
  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (auto index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (auto &i : field.getRegion(RGN_NOBNDRY)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, (((nx - 2) * (ny - 2)) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRGN_NOX) {
  Field2D field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0});
  test_indices.insert({0, 1});
  test_indices.insert({1, 0});
  test_indices.insert({1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({1, 0});
  region_indices.insert({1, 1});
  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (auto index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (auto &i : field.getRegion(RGN_NOX)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, (((nx - 2) * ny) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRGN_NOY) {
  Field2D field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0});
  test_indices.insert({0, 1});
  test_indices.insert({1, 0});
  test_indices.insert({1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({0, 1});
  region_indices.insert({1, 1});
  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (auto index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (auto &i : field.getRegion(RGN_NOY)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * (ny - 2)) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRGN_NOZ) {
  Field2D field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0});
  test_indices.insert({0, 1});
  test_indices.insert({1, 0});
  test_indices.insert({1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({0, 0});
  region_indices.insert({0, 1});
  region_indices.insert({1, 0});
  region_indices.insert({1, 1});
  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (auto index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (auto &i : field.getRegion(RGN_NOZ)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * ny) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRGN_XGUARDS) {
  Field2D field;

  field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0});
  test_indices.insert({0, 1});
  test_indices.insert({1, 0});
  test_indices.insert({1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({0, 1});

  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field.getRegion("RGN_XGUARDS")) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((2 * (ny - 2)) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRGN_YGUARDS) {
  Field2D field;

  field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0});
  test_indices.insert({0, 1});
  test_indices.insert({1, 0});
  test_indices.insert({1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({1, 0});

  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field.getRegion("RGN_YGUARDS")) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, (((nx - 2) * 2) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRGN_ZGUARDS) {
  Field2D field;

  field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0});
  test_indices.insert({0, 1});
  test_indices.insert({1, 0});
  test_indices.insert({1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;

  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field.getRegion("RGN_ZGUARDS")) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx - 2) * (ny - 2) * 0 - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(Field2DTest, IterateOverRGN_NOCORNERS) {
  Field2D field;

  field = 1.0;

  const BoutReal sentinel = -99.0;

  // We use a set in case for some reason the iterator doesn't visit
  // each point in the order we expect.
  std::set<std::vector<int>> test_indices;
  test_indices.insert({0, 0});
  test_indices.insert({0, 1});
  test_indices.insert({1, 0});
  test_indices.insert({1, 1});

  // This is the set of indices actually inside the region we want
  std::set<std::vector<int>> region_indices;
  region_indices.insert({0, 1});
  region_indices.insert({1, 0});
  region_indices.insert({1, 1});

  const int num_sentinels = region_indices.size();

  // Assign sentinel value to watch out for to our chosen points
  for (const auto index : test_indices) {
    field(index[0], index[1]) = sentinel;
  }

  int found_sentinels = 0;
  BoutReal sum = 0.0;
  std::set<std::vector<int>> result_indices;

  for (const auto &i : field.getRegion("RGN_NOCORNERS")) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * ny - 4) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}


TEST_F(Field2DTest, Indexing) {
  Field2D field;

  field.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      field(i, j) = i + j;
    }
  }

  EXPECT_DOUBLE_EQ(field(2, 2), 4);
}

TEST_F(Field2DTest, IndexingAs3D) {
  Field2D field;

  field.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        field(i, j, k) = i + j + k;
      }
    }
  }

  EXPECT_DOUBLE_EQ(field(2, 2), 4 + nz -1);
}

TEST_F(Field2DTest, ConstIndexingAs3D) {
  const Field2D field = 3.0;
  Field2D field2;
  field2.allocate();
  
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        field2(i, j, k) = field(i, j, k) + i + j + k;
      }
    }
  }

  EXPECT_DOUBLE_EQ(field2(2, 2), 3 + 4 + nz - 1);
}

TEST_F(Field2DTest, IndexingInd3D) {
  Field2D field;

  field.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        field(i, j, k) = i + j + k;
      }
    }
  }

  Ind3D ind{(2*ny + 2)*nz + 2};

  EXPECT_DOUBLE_EQ(field[ind], 10);
}

TEST_F(Field2DTest, ConstIndexingInd3D) {
  Field2D field1;

  field1.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        field1(i, j, k) = i + j + k;
      }
    }
  }

  const Field2D field2{field1};

  Ind3D ind{(2*ny + 2)*nz + 2};

  EXPECT_DOUBLE_EQ(field2[ind], 10);
}

#if CHECK > 2
TEST_F(Field2DTest, CheckNotEmpty) {
  Field2D field;

  EXPECT_THROW(field(0, 0), BoutException);
}

TEST_F(Field2DTest, BoundsCheck) {
  Field2D field;

  field = 1.0;
  EXPECT_THROW(field(-1, 0), BoutException);
  EXPECT_THROW(field(0, -1), BoutException);
  EXPECT_THROW(field(nx, 0), BoutException);
  EXPECT_THROW(field(0, ny), BoutException);
}

TEST_F(Field2DTest, ConstCheckNotEmpty) {
  const Field2D field;

  EXPECT_THROW(field(0, 0), BoutException);
}

TEST_F(Field2DTest, ConstBoundsCheck) {
  const Field2D field = 1.0;

  EXPECT_THROW(field(-1, 0), BoutException);
  EXPECT_THROW(field(0, -1), BoutException);
  EXPECT_THROW(field(nx, 0), BoutException);
  EXPECT_THROW(field(0, ny), BoutException);
}

TEST_F(Field2DTest, CheckData) {
  Field2D field;

  EXPECT_THROW(checkData(field), BoutException);

  field = 1.0;

  EXPECT_NO_THROW(checkData(field));

  field(1, 1) = std::nan("");

  EXPECT_THROW(checkData(field), BoutException);
  
  field = 1.0;
  field(0, 0) = std::nan("");

  EXPECT_NO_THROW(checkData(field));
  EXPECT_NO_THROW(checkData(field, RGN_NOBNDRY));
  EXPECT_THROW(checkData(field, RGN_ALL), BoutException);

}

#if CHECK > 0
TEST_F(Field2DTest, DoneComms) {
  Field2D field = 1.0;
  field.bndry_xin = false;
  field.bndry_xout = false;
  field.bndry_yup = false;
  field.bndry_ydown = false;

  EXPECT_THROW(field.bndryValid(), BoutException);
  field.doneComms();
  EXPECT_EQ(field.bndryValid(), true);
}
#endif

TEST_F(Field2DTest, InvalidateGuards) {
  Field2D field;
  field.allocate(); // Calls invalidateGuards
  field = 1.0;      // Sets everywhere including boundaries

  const int nmesh = nx * ny;

  int sum = 0;
  for (const auto &i : field) {
    field[i] = 0.0; // Reset field value
    sum++;
  }
  EXPECT_EQ(sum, nmesh); // Field operator= hasn't been broken by invalidateGuards

  // Count the number of non-boundary points
  sum = 0;
  for (const auto &i : field.getRegion(RGN_NOBNDRY)) {
    field[i] = 0.0; // Reset field value
    sum++;
  }
  const int nbndry = nmesh - sum;

  auto localmesh = field.getMesh();
  EXPECT_NO_THROW(checkData(field(0, 0)));
  EXPECT_NO_THROW(checkData(field(localmesh->xstart, localmesh->ystart)));

  invalidateGuards(field);

  EXPECT_THROW(checkData(field(0, 0)), BoutException);
  EXPECT_NO_THROW(checkData(field(localmesh->xstart, localmesh->ystart)));

  sum = 0;
  for (const auto &i : field) {
    if (!finite(field[i]))
      sum++;
  }
  EXPECT_EQ(sum, nbndry);
}

#endif // CHECK > 2

TEST_F(Field2DTest, CreateFromBoutReal) {
  Field2D field(1.0);

  EXPECT_TRUE(IsFieldEqual(field, 1.0));
}

TEST_F(Field2DTest, CreateFromField2D) {
  Field2D field(99.0);
  Field2D result(field);

  EXPECT_TRUE(IsFieldEqual(result, 99.0));
}

TEST_F(Field2DTest, AssignFromBoutReal) {
  Field2D field;

  field = 2.0;

  EXPECT_TRUE(IsFieldEqual(field, 2.0));
}

TEST_F(Field2DTest, AssignFromInvalid) {
  Field2D field;

  EXPECT_NO_THROW(field = std::nan(""));
  EXPECT_TRUE(IsFieldEqual(field, std::nan("")));
}

TEST_F(Field2DTest, UnaryMinus) {
  Field2D field;

  field = 2.0;
  field = -field;

  EXPECT_TRUE(IsFieldEqual(field, -2.0));
}

TEST_F(Field2DTest, AddEqualsBoutReal) {
  Field2D a;

  a = 1.0;
  a += 5.0;

  EXPECT_TRUE(IsFieldEqual(a, 6.0));

  // Check case where field is not unique
  auto c = a;
  c += 5.0;

  EXPECT_TRUE(IsFieldEqual(a, 6.0));
  EXPECT_TRUE(IsFieldEqual(c, 11.0));
}

TEST_F(Field2DTest, AddEqualsField2D) {
  Field2D a, b;

  a = 2.0;
  b = 3.0;
  a += b;

  EXPECT_TRUE(IsFieldEqual(a, 5.0));

  // Check case where field is not unique
  auto c = a;
  c += b;

  EXPECT_TRUE(IsFieldEqual(a, 5.0));
  EXPECT_TRUE(IsFieldEqual(c, 8.0));
}

TEST_F(Field2DTest, AddField2DBoutReal) {
  Field2D a, b;

  a = 1.0;
  b = a + 2.0;

  EXPECT_TRUE(IsFieldEqual(b, 3.0));
}

TEST_F(Field2DTest, AddBoutRealField2D) {
  Field2D a, b;

  a = 1.0;
  b = 3.0 + a;

  EXPECT_TRUE(IsFieldEqual(b, 4.0));
}

TEST_F(Field2DTest, AddField2DField2D) {
  Field2D a, b, c;

  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsFieldEqual(c, 3.0));
}

TEST_F(Field2DTest, MultiplyEqualsBoutReal) {
  Field2D a;

  a = 2.0;
  a *= 1.5;

  EXPECT_TRUE(IsFieldEqual(a, 3.0));

  // Check case where field is not unique
  auto c = a;
  c *= 1.5;

  EXPECT_TRUE(IsFieldEqual(a, 3.0));
  EXPECT_TRUE(IsFieldEqual(c, 4.5));
}

TEST_F(Field2DTest, MultiplyEqualsField2D) {
  Field2D a, b;

  a = 2.5;
  b = 4.0;
  a *= b;

  EXPECT_TRUE(IsFieldEqual(a, 10.0));

  // Check case where field is not unique
  auto c = a;
  c *= b;

  EXPECT_TRUE(IsFieldEqual(a, 10.0));
  EXPECT_TRUE(IsFieldEqual(c, 40.0));
}

TEST_F(Field2DTest, MultiplyField2DBoutReal) {
  Field2D a, b;

  a = 1.5;
  b = a * 2.0;

  EXPECT_TRUE(IsFieldEqual(b, 3.0));
}

TEST_F(Field2DTest, MultiplyBoutRealField2D) {
  Field2D a, b;

  a = 2.5;
  b = 3.0 * a;

  EXPECT_TRUE(IsFieldEqual(b, 7.5));
}

TEST_F(Field2DTest, MultiplyField2DField2D) {
  Field2D a, b, c;

  a = 4.0;
  b = 8.0;
  c = a * b;

  EXPECT_TRUE(IsFieldEqual(c, 32.0));
}

TEST_F(Field2DTest, SubtractEqualsBoutReal) {
  Field2D a;

  a = 1.0;
  a -= 5.0;

  EXPECT_TRUE(IsFieldEqual(a, -4.0));

  // Check case where field is not unique
  auto c = a;
  c -= 5.0;

  EXPECT_TRUE(IsFieldEqual(a, -4.0));
  EXPECT_TRUE(IsFieldEqual(c, -9.0));
}

TEST_F(Field2DTest, SubtractEqualsField2D) {
  Field2D a, b;

  a = 2.0;
  b = 7.0;
  a -= b;

  EXPECT_TRUE(IsFieldEqual(a, -5.0));

  // Check case where field is not unique
  auto c = a;
  c -= b;

  EXPECT_TRUE(IsFieldEqual(a, -5.0));
  EXPECT_TRUE(IsFieldEqual(c, -12.0));
}

TEST_F(Field2DTest, SubtractField2DBoutReal) {
  Field2D a, b;

  a = 10.0;
  b = a - 2.0;

  EXPECT_TRUE(IsFieldEqual(b, 8.0));
}

TEST_F(Field2DTest, SubtractBoutRealField2D) {
  Field2D a, b;

  a = 10.0;
  b = 3.0 - a;

  EXPECT_TRUE(IsFieldEqual(b, -7.0));
}

TEST_F(Field2DTest, SubtractField2DField2D) {
  Field2D a, b, c;

  a = 10.0;
  b = 20.0;
  c = a - b;

  EXPECT_TRUE(IsFieldEqual(c, -10.0));
}

TEST_F(Field2DTest, DivideEqualsBoutReal) {
  Field2D a;

  a = 2.5;
  a /= 5.0;

  EXPECT_TRUE(IsFieldEqual(a, 0.5));

  // Check case where field is not unique
  auto c = a;
  c /= 5.0;

  EXPECT_TRUE(IsFieldEqual(a, 0.5));
  EXPECT_TRUE(IsFieldEqual(c, 0.1));
}

TEST_F(Field2DTest, DivideEqualsField2D) {
  Field2D a, b;

  a = 5.0;
  b = 2.5;
  a /= b;

  EXPECT_TRUE(IsFieldEqual(a, 2.0));

  // Check case where field is not unique
  auto c = a;
  c /= b;

  EXPECT_TRUE(IsFieldEqual(a, 2.0));
  EXPECT_TRUE(IsFieldEqual(c, 0.8));
}

TEST_F(Field2DTest, DivideField2DBoutReal) {
  Field2D a, b;

  a = 3.0;
  b = a / 2.0;

  EXPECT_TRUE(IsFieldEqual(b, 1.5));
}

TEST_F(Field2DTest, DivideBoutRealField2D) {
  Field2D a, b;

  a = 2.5;
  b = 10.0 / a;

  EXPECT_TRUE(IsFieldEqual(b, 4.0));
}

TEST_F(Field2DTest, DivideField2DField2D) {
  Field2D a, b, c;

  a = 32.0;
  b = 8.0;
  c = a / b;

  EXPECT_TRUE(IsFieldEqual(c, 4.0));
}

TEST_F(Field2DTest, PowBoutRealField2D) {
  Field2D a, b;
  a = 5.0;
  b = pow(2.0, a);

  EXPECT_TRUE(IsFieldEqual(b, 32.0));
}

TEST_F(Field2DTest, PowField2DBoutReal) {
  Field2D a, b;
  a = 5.0;
  b = pow(a, 2.0);

  EXPECT_TRUE(IsFieldEqual(b, 25.0));
}

TEST_F(Field2DTest, PowField2DField2D) {
  Field2D a, b, c;
  a = 2.0;
  b = 6.0;
  c = pow(a, b);

  EXPECT_TRUE(IsFieldEqual(c, 64.0));
}

TEST_F(Field2DTest, Sqrt) {
  Field2D field;

  field = 16.0;
  EXPECT_TRUE(IsFieldEqual(sqrt(field), 4.0));
}

TEST_F(Field2DTest, Abs) {
  Field2D field;

  field = -31.0;
  EXPECT_TRUE(IsFieldEqual(abs(field), 31.0));
}

TEST_F(Field2DTest, Exp) {
  Field2D field;

  field = 2.5;
  const BoutReal expected = 12.182493960703473;
  EXPECT_TRUE(IsFieldEqual(exp(field), expected));
}

TEST_F(Field2DTest, Log) {
  Field2D field;

  field = 12.182493960703473;
  const BoutReal expected = 2.5;
  EXPECT_TRUE(IsFieldEqual(log(field), expected));
}

TEST_F(Field2DTest, LogExp) {
  Field2D field;

  field = 2.5;
  const BoutReal expected = 2.5;
  EXPECT_TRUE(IsFieldEqual(log(exp(field)), expected));
}

TEST_F(Field2DTest, Sin) {
  Field2D field;

  field = PI / 2.0;
  EXPECT_TRUE(IsFieldEqual(sin(field), 1.0));

  field = PI;
  EXPECT_TRUE(IsFieldEqual(sin(field), 0.0));
}

TEST_F(Field2DTest, Cos) {
  Field2D field;

  field = PI / 2.0;
  EXPECT_TRUE(IsFieldEqual(cos(field), 0.0));

  field = PI;
  EXPECT_TRUE(IsFieldEqual(cos(field), -1.0));
}

TEST_F(Field2DTest, Tan) {
  Field2D field;

  field = PI / 4.0;
  EXPECT_TRUE(IsFieldEqual(tan(field), 1.0));

  field = PI;
  EXPECT_TRUE(IsFieldEqual(tan(field), 0.0));
}

TEST_F(Field2DTest, Sinh) {
  Field2D field;

  field = 1.0;
  const BoutReal expected = 1.1752011936438014;
  EXPECT_TRUE(IsFieldEqual(sinh(field), expected));

  field = -1.0;
  EXPECT_TRUE(IsFieldEqual(sinh(field), -expected));
}

TEST_F(Field2DTest, Cosh) {
  Field2D field;

  field = 1.0;
  const BoutReal expected = 1.5430806348152437;
  EXPECT_TRUE(IsFieldEqual(cosh(field), expected));

  field = -1.0;
  EXPECT_TRUE(IsFieldEqual(cosh(field), expected));
}

TEST_F(Field2DTest, Tanh) {
  Field2D field;

  field = 1.0;
  const BoutReal expected = 0.761594155955764;
  EXPECT_TRUE(IsFieldEqual(tanh(field), expected));

  field = -1.0;
  EXPECT_TRUE(IsFieldEqual(tanh(field), -expected));
}

TEST_F(Field2DTest, Floor) {
  Field2D field;

  field = 50.0;
  field(1, 1) = 49.9;
  field(2, 3) = -20;

  const BoutReal floor_value = 50.0;

  EXPECT_TRUE(IsFieldEqual(floor(field, floor_value), floor_value));
}

TEST_F(Field2DTest, Min) {
  Field2D field;

  field = 50.0;
  field(0, 0) = -99.0;
  field(1, 1) = 60.0;
  field(1, 2) = 40.0;
  field(2, 4) = 99.0;

  // min doesn't include guard cells
  const BoutReal min_value = 40.0;

  EXPECT_EQ(min(field, false), min_value);
  EXPECT_EQ(min(field, false, "RGN_ALL"), -99.0);
  EXPECT_EQ(min(field, true, "RGN_ALL"), -99.0);
}

TEST_F(Field2DTest, Max) {
  Field2D field;

  field = 50.0;
  field(0, 0) = -99.0;
  field(1, 1) = 40.0;
  field(1, 2) = 60.0;
  field(2, 4) = 99.0;

  // max doesn't include guard cells
  const BoutReal max_value = 60.0;

  EXPECT_EQ(max(field, false), max_value);
  EXPECT_EQ(max(field, false, "RGN_ALL"), 99.0);
  EXPECT_EQ(max(field, true, "RGN_ALL"), 99.0);
}

TEST_F(Field2DTest, Swap) {
  WithQuietOutput quiet{output_info};

  Field2D first(1., mesh_staggered);

  first.setLocation(CELL_XLOW);

  ddt(first) = 1.1;

  // Mesh for second field
  constexpr int second_nx = Field2DTest::nx + 2;
  constexpr int second_ny = Field2DTest::ny + 2;
  constexpr int second_nz = Field2DTest::nz + 2;

  FakeMesh second_mesh{second_nx, second_ny, second_nz};
  second_mesh.setCoordinates(nullptr);
  second_mesh.StaggerGrids = false;
  second_mesh.createDefaultRegions();

  // Second field
  Field2D second(2., &second_mesh);

  ddt(second) = 2.4;

  // Basic sanity check
  EXPECT_TRUE(IsFieldEqual(first, 1.0));
  EXPECT_TRUE(IsFieldEqual(second, 2.0));

  // swap is marked noexcept, so absolutely should not throw!
  ASSERT_NO_THROW(swap(first, second));

  // Values
  EXPECT_TRUE(IsFieldEqual(first, 2.0));
  EXPECT_TRUE(IsFieldEqual(second, 1.0));

  EXPECT_TRUE(IsFieldEqual(ddt(first), 2.4));
  EXPECT_TRUE(IsFieldEqual(ddt(second), 1.1));

  // Mesh properties
  EXPECT_EQ(first.getMesh(), &second_mesh);
  EXPECT_EQ(second.getMesh(), mesh_staggered);

  EXPECT_EQ(first.getNx(), second_nx);
  EXPECT_EQ(first.getNy(), second_ny);
  EXPECT_EQ(first.getNz(), 1);

  EXPECT_EQ(second.getNx(), Field2DTest::nx);
  EXPECT_EQ(second.getNy(), Field2DTest::ny);
  EXPECT_EQ(second.getNz(), 1);

  EXPECT_EQ(first.getLocation(), CELL_CENTRE);
  EXPECT_EQ(second.getLocation(), CELL_XLOW);

  // We don't check the boundaries, but the data is protected and
  // there are no inquiry functions
}

TEST_F(Field2DTest, MoveCtor) {
  // First field
  Field2D first(1., mesh_staggered);

  first.setLocation(CELL_XLOW);

  ddt(first) = 1.1;

  // Second field
  Field2D second{std::move(first)};

  // Values
  EXPECT_TRUE(IsFieldEqual(second, 1.0));

  EXPECT_TRUE(IsFieldEqual(ddt(second), 1.1));

  // Mesh properties
  EXPECT_EQ(second.getMesh(), mesh_staggered);

  EXPECT_EQ(second.getNx(), Field2DTest::nx);
  EXPECT_EQ(second.getNy(), Field2DTest::ny);
  EXPECT_EQ(second.getNz(), 1);

  EXPECT_EQ(second.getLocation(), CELL_XLOW);

  // We don't check the boundaries, but the data is protected and
  // there are no inquiry functions
}

TEST_F(Field2DTest, FillField) {
  Field2D f{mesh};

  fillField(f, {{1., 1., 1., 1., 1.}, {1., 1., 1., 1., 1.}, {1., 1., 1., 1., 1.}});

  EXPECT_TRUE(IsFieldEqual(f, 1.));

  fillField(f, {{0., 1., 2., 3., 4.}, {0., 1., 2., 3., 4.}, {0., 1., 2., 3., 4.}});

  Field2D g{mesh};
  g.allocate();
  BOUT_FOR_SERIAL(i, g.getRegion("RGN_ALL")) { g[i] = i.y(); }

  EXPECT_TRUE(IsFieldEqual(f, g));
}

TEST_F(Field2DTest, OperatorEqualsField2D) {
  Field2D field;

  // Create field with non-default arguments so we can check they get copied
  // to 'field'.
  // Note that Aligned y-direction type is not really allowed for Field2D, but
  // we don't check anywhere at the moment.
  Field2D field2{mesh_staggered, CELL_XLOW, {YDirectionType::Aligned, ZDirectionType::Average}};

  field = field2;

  EXPECT_TRUE(areFieldsCompatible(field, field2));
  EXPECT_EQ(field.getMesh(), field2.getMesh());
  EXPECT_EQ(field.getLocation(), field2.getLocation());
  EXPECT_EQ(field.getDirectionY(), field2.getDirectionY());
  EXPECT_EQ(field.getDirectionZ(), field2.getDirectionZ());
}

TEST_F(Field2DTest, EmptyFrom) {
  // Create field with non-default arguments so we can check they get copied
  // to 'field2'.
  // Note that Aligned y-direction type is not really allowed for Field2D, but
  // we don't check anywhere at the moment.
  Field2D field{mesh_staggered, CELL_XLOW, {YDirectionType::Aligned, ZDirectionType::Average}};
  field = 5.;

  Field2D field2{emptyFrom(field)};
  EXPECT_EQ(field2.getMesh(), mesh_staggered);
  EXPECT_EQ(field2.getLocation(), CELL_XLOW);
  EXPECT_EQ(field2.getDirectionY(), YDirectionType::Aligned);
  EXPECT_EQ(field2.getDirectionZ(), ZDirectionType::Average);
  EXPECT_TRUE(field2.isAllocated());
}

TEST_F(Field2DTest, ZeroFrom) {
  // Create field with non-default arguments so we can check they get copied
  // to 'field2'.
  // Note that Aligned y-direction type is not really allowed for Field2D, but
  // we don't check anywhere at the moment.
  Field2D field{mesh_staggered, CELL_XLOW, {YDirectionType::Aligned, ZDirectionType::Average}};
  field = 5.;

  Field2D field2{zeroFrom(field)};
  EXPECT_EQ(field2.getMesh(), mesh_staggered);
  EXPECT_EQ(field2.getLocation(), CELL_XLOW);
  EXPECT_EQ(field2.getDirectionY(), YDirectionType::Aligned);
  EXPECT_EQ(field2.getDirectionZ(), ZDirectionType::Average);
  EXPECT_TRUE(field2.isAllocated());
  EXPECT_TRUE(IsFieldEqual(field2, 0.));
}
#pragma GCC diagnostic pop
