#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"
#include "unused.hxx"

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

    mesh->createDefaultRegions();
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
  int test_nx = 4;
  int test_ny = 8;
  int test_nz = 9;

  FakeMesh *fieldmesh = new FakeMesh(test_nx, test_ny, test_nz);

  Field3D field(fieldmesh);

  field.allocate();

  EXPECT_EQ(field.getNx(), test_nx);
  EXPECT_EQ(field.getNy(), test_ny);
  EXPECT_EQ(field.getNz(), test_nz);

  delete fieldmesh;
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

TEST_F(Field3DTest, IterateOverRGN_ALL_SDI) {
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

  for (const auto &i : field.sdi_region(RGN_ALL)) {
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

TEST_F(Field3DTest, IterateOverRGN_ALL_SDI_NonRangeBased) {
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

  for (auto i = field.sdi_region(RGN_ALL).begin(), end = field.sdi_region(RGN_ALL).end();
       i < end; ++i) {
    sum += field(*i);
    if (field(*i) == sentinel) {
      result_indices.insert({i->x(), i->y(), i->z()});
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

TEST_F(Field3DTest, IterateOverRGN_NOBNDRY_SDI) {
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

  for (const auto &i : field.sdi_region(RGN_NOBNDRY)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y(), i.z()});
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

TEST_F(Field3DTest, IterateOverRGN_NOX_SDI) {
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

  for (const auto &i : field.sdi_region(RGN_NOX)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y(), i.z()});
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

TEST_F(Field3DTest, IterateOverRGN_NOY_SDI) {
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

  for (const auto &i : field.sdi_region(RGN_NOY)) {
    sum += field[i];
    if (field[i] == sentinel) {
      result_indices.insert({i.x(), i.y(), i.z()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * (ny - 2) * nz) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
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

TEST_F(Field3DTest, AssignFromBoutReal) {
  Field3D field;

  field = 2.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(field, 2.0));
}

TEST_F(Field3DTest, AssignFromField2D) {
  Field3D field;
  Field2D field2(2.0);

  field = field2;

  EXPECT_TRUE(IsField3DEqualBoutReal(field, 2.0));
}

TEST_F(Field3DTest, AssignFromField3D) {
  Field3D field, field2;

  field2 = -99.0;
  field = field2;

  EXPECT_TRUE(IsField3DEqualBoutReal(field, -99.0));
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
}

TEST_F(Field3DTest, AddEqualsField2D) {
  Field3D a;
  Field2D b;

  a = 2.0;
  b = 3.0;
  a += b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 5.0));
}

TEST_F(Field3DTest, AddEqualsField3D) {
  Field3D a, b;

  a = 2.0;
  b = 3.0;
  a += b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 5.0));
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

TEST_F(Field3DTest, MultiplyEqualsBoutReal) {
  Field3D a;

  a = 2.0;
  a *= 1.5;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 3.0));
}

TEST_F(Field3DTest, MultiplyEqualsField2D) {
  Field3D a;
  Field2D b;

  a = 2.5;
  b = 4.0;
  a *= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 10.0));
}

TEST_F(Field3DTest, MultiplyEqualsField3D) {
  Field3D a, b;

  a = 2.5;
  b = 4.0;
  a *= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 10.0));
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

TEST_F(Field3DTest, SubtractEqualsBoutReal) {
  Field3D a;

  a = 1.0;
  a -= 5.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, -4.0));
}

TEST_F(Field3DTest, SubtractEqualsField2D) {
  Field3D a;
  Field2D b;

  a = 2.0;
  b = 7.0;
  a -= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, -5.0));
}

TEST_F(Field3DTest, SubtractEqualsField3D) {
  Field3D a, b;

  a = 2.0;
  b = 7.0;
  a -= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, -5.0));
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

TEST_F(Field3DTest, DivideEqualsBoutReal) {
  Field3D a;

  a = 2.5;
  a /= 5.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 0.5));
}

TEST_F(Field3DTest, DivideEqualsField2D) {
  Field3D a;
  Field2D b;

  a = 5.0;
  b = 2.5;
  a /= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 2.0));
}

TEST_F(Field3DTest, DivideEqualsField3D) {
  Field3D a, b;

  a = 5.0;
  b = 2.5;
  a /= b;

  EXPECT_TRUE(IsField3DEqualBoutReal(a, 2.0));
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

  EXPECT_TRUE(IsField3DEqualBoutReal(min(field, false), min_value));
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

  EXPECT_TRUE(IsField3DEqualBoutReal(max(field, false), max_value));
}

TEST_F(Field3DTest, DC) {
  Field3D field;

  field = 1.0;
  for (const auto& i : field) {
    field[i] = i.z;
  }

  EXPECT_TRUE(IsField2DEqualBoutReal(DC(field), 3.0));
}
