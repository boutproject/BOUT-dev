// We know stuff might be deprecated, but we still want to test it
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "fieldperp.hxx"
#include "test_extras.hxx"
#include "output.hxx"
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

/// Test fixture to make sure the global mesh is our fake one
using FieldPerpTest = FakeMeshFixture;

TEST_F(FieldPerpTest, Allocate) {
  FieldPerp field;

  EXPECT_FALSE(field.isAllocated());

  field.allocate();

  EXPECT_TRUE(field.isAllocated());

  int counter = 0;
  for (const auto &i : field) {
    field[i] = 1.; // Hits Array bounds checking
    counter++;
  }
  EXPECT_EQ(counter, FieldPerpTest::nx * FieldPerpTest::nz);
}

TEST_F(FieldPerpTest, IsFinite) {
  FieldPerp field;

  EXPECT_FALSE(finite(field));

  field = 1.0;

  EXPECT_TRUE(finite(field));

  field(1, 1) = std::nan("");

  EXPECT_FALSE(finite(field));
}

TEST_F(FieldPerpTest, SliceXZ) {
  Field3D masterField;
#if CHECK > 0
  EXPECT_THROW(sliceXZ(masterField, 0), BoutException);
#endif
  masterField = 1.0;
  EXPECT_NO_THROW(sliceXZ(masterField, 0));

  const int yindex = 2;
  auto result = sliceXZ(masterField, yindex);

  EXPECT_EQ(result.getIndex(), yindex);
  EXPECT_TRUE(areFieldsCompatible(masterField, result));

  for (const auto &i : result) {
    EXPECT_EQ(result[i], 1.0);
  }
}

TEST_F(FieldPerpTest, GetGridSizes) {
  FieldPerp field;

  field.allocate();

  EXPECT_EQ(field.getNx(), FieldPerpTest::nx);
  EXPECT_EQ(field.getNy(), 1);
  EXPECT_EQ(field.getNz(), FieldPerpTest::nz);
}

TEST_F(FieldPerpTest, GetSetIndex) {
  FieldPerp field;
  EXPECT_EQ(field.getIndex(), -1);
  field.setIndex(1);
  EXPECT_EQ(field.getIndex(), 1);
  FieldPerp field2(field);
  field.setIndex(2);
  EXPECT_EQ(field2.getIndex(), 1);
  EXPECT_EQ(field.getIndex(), 2);

  // Currently don't bounds check setIndex so the following is
  // allowed -- make this explicit by testing
  EXPECT_NO_THROW(field.setIndex(-10));
  EXPECT_EQ(field.getIndex(), -10);
}

TEST_F(FieldPerpTest, CreateOnGivenMesh) {
  WithQuietOutput quiet{output_info};

  int test_nx = FieldPerpTest::nx + 2;
  int test_ny = FieldPerpTest::ny + 2;
  int test_nz = FieldPerpTest::nz + 2;

  FakeMesh fieldmesh{test_nx, test_ny, test_nz};
  fieldmesh.setCoordinates(nullptr);
  fieldmesh.createDefaultRegions();

  FieldPerp field{&fieldmesh};

  field.allocate();

  EXPECT_EQ(field.getNx(), test_nx);
  EXPECT_EQ(field.getNy(), 1);
  EXPECT_EQ(field.getNz(), test_nz);
}

TEST_F(FieldPerpTest, CopyCheckFieldmesh) {
  WithQuietOutput quiet{output_info};

  int test_nx = FieldPerpTest::nx + 2;
  int test_ny = FieldPerpTest::ny + 2;
  int test_nz = FieldPerpTest::nz + 2;

  FakeMesh fieldmesh{test_nx, test_ny, test_nz};
  fieldmesh.setCoordinates(nullptr);
  fieldmesh.createDefaultRegions();

  FieldPerp field{&fieldmesh};
  field = 1.0;

  FieldPerp field2{field};

  EXPECT_EQ(field2.getNx(), test_nx);
  EXPECT_EQ(field2.getNy(), 1);
  EXPECT_EQ(field2.getNz(), test_nz);
  EXPECT_TRUE(areFieldsCompatible(field, field2));
}

#if CHECK > 0
TEST_F(FieldPerpTest, CreateOnNullMesh) {
  auto old_mesh = mesh;
  mesh = nullptr;

  FieldPerp field;

  EXPECT_EQ(field.getNx(), -1);
  EXPECT_EQ(field.getNy(), 1);
  EXPECT_EQ(field.getNz(), -1);

  mesh = old_mesh;

  field.allocate();

  EXPECT_EQ(field.getNx(), FieldPerpTest::nx);
  EXPECT_EQ(field.getNy(), 1);
  EXPECT_EQ(field.getNz(), FieldPerpTest::nz);
}
#endif

#if CHECK > 0 && CHECK <= 2
// We only want to run this test in a certain range of CHECK as we're
// checking some behaviour that is only enabled for CHECK above 0
// but there are checks that will throw before reaching these lines if
// check is greater than 2, so the test only makes sense in a certain range.
TEST_F(FieldPerpTest, CreateCopyOnNullMesh) {
  // Whilst the declaration of field below looks like it should create a FieldPerp
  // without a mesh, it in fact will result in a FieldPerp associated with the
  // global mesh as we end up calling the Field constructor that forces this.
  // Hence, to test the case of copying a field without a mesh we have to
  // temporarily hide the global mesh, before restoring it later.
  auto old_mesh = mesh;
  mesh = nullptr;

  FieldPerp field;
  // If CHECK > 2 then the following will throw due to the data
  // block in field not being allocated. We can't allocate as that
  // would force field to have a mesh associated with it.
  FieldPerp field2(field);

  EXPECT_EQ(field2.getNx(), -1);
  EXPECT_EQ(field2.getNy(), 1);
  EXPECT_EQ(field2.getNz(), -1);

  mesh = old_mesh;
  field2.allocate();

  EXPECT_EQ(field2.getNx(), FieldPerpTest::nx);
  EXPECT_EQ(field2.getNy(), 1);
  EXPECT_EQ(field2.getNz(), FieldPerpTest::nz);
}
#endif

TEST_F(FieldPerpTest, SetGetLocation) {
  FieldPerp field(mesh_staggered);

  field.setLocation(CELL_XLOW);
  EXPECT_EQ(field.getLocation(), CELL_XLOW);

  field.setLocation(CELL_DEFAULT);
  EXPECT_EQ(field.getLocation(), CELL_CENTRE);

  EXPECT_THROW(field.setLocation(CELL_VSHIFT), BoutException);
}

TEST_F(FieldPerpTest, SetGetLocationNonStaggered) {
  FieldPerp field;

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
TEST_F(FieldPerpTest, IterateOverWholeField) {
  FieldPerp field(mesh);

  field.allocate();

  // Basic test first: do we visit the correct number of elements?
  int count = 0;
  for (auto &UNUSED(i) : field) {
    ++count;
  }

  // If this fails, no point doing second part
  ASSERT_EQ(count, nx * nz);

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
      result_indices.insert({i.x(), i.z()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * nz) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(test_indices == result_indices);
}

TEST_F(FieldPerpTest, IterateOverRGN_ALL) {
  FieldPerp field(mesh);
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
      result_indices.insert({i.x(), i.z()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nx * nz) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(test_indices == result_indices);
}

TEST_F(FieldPerpTest, IterateOverRGN_NOZ) {
  FieldPerp field(mesh);
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
      result_indices.insert({i.x(), i.z()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, (nx * nz - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(FieldPerpTest, IterateOverRGN_NOX) {
  FieldPerp field(mesh);
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
      result_indices.insert({i.x(), i.z()});
      ++found_sentinels;
    }
  }

  EXPECT_EQ(found_sentinels, num_sentinels);
  EXPECT_EQ(sum, ((nz * (nx - 2)) - num_sentinels) + (num_sentinels * sentinel));
  EXPECT_TRUE(region_indices == result_indices);
}

TEST_F(FieldPerpTest, Indexing) {
  FieldPerp field;

  field.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      field(i, j) = i + j;
    }
  }

  EXPECT_DOUBLE_EQ(field(2, 2), 4);
}

TEST_F(FieldPerpTest, IndexingAs3D) {
  FieldPerp field;

  field.allocate();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        field(i, j, k) = i + j + k;
      }
    }
  }

  EXPECT_DOUBLE_EQ(field(2, 2), 4 + ny - 1);
}

TEST_F(FieldPerpTest, IndexingWithIndPerp) {
  FieldPerp field(0.0);
  const BoutReal sentinel = 2.0;
  field(1, 5) = sentinel;
  IndPerp firstPoint(1 * nz + 5, 1, nz);
  EXPECT_DOUBLE_EQ(field[firstPoint], sentinel);
  IndPerp secondPoint((nx - 1) * nz + 1, 1, nz);
  field[secondPoint] = -sentinel;
  EXPECT_DOUBLE_EQ(field(nx - 1, 1), -sentinel);
}

TEST_F(FieldPerpTest, ConstIndexingWithIndPerp) {
  FieldPerp field(0.0);
  const BoutReal sentinel = 2.0;
  field(1, 5) = sentinel;
  const IndPerp firstPoint(1 * nz + 5, 1, nz);
  EXPECT_DOUBLE_EQ(field[firstPoint], sentinel);
  const IndPerp secondPoint((nx - 1) * nz + 1, 1, nz);
  field[secondPoint] = -sentinel;
  EXPECT_DOUBLE_EQ(field(nx - 1, 1), -sentinel);
}

TEST_F(FieldPerpTest, IndexingWithInd3D) {
  FieldPerp field(0.0);
  field.setIndex(2);
  const BoutReal sentinel = 2.0;
  int ix = 1, iy = field.getIndex(), iz = 4;
  field(ix, iz) = sentinel;
  Ind3D firstPoint(iz + nz * (iy + ny * ix), ny, nz);
  EXPECT_DOUBLE_EQ(field[firstPoint], sentinel);
  field[firstPoint] = -sentinel;
  EXPECT_DOUBLE_EQ(field(ix, iz), -sentinel);
#if CHECK > 2
  iy++;
  Ind3D secondPoint(iz + nz * (iy + ny * ix), ny, nz);
  EXPECT_THROW(field[secondPoint], BoutException);
#endif
}

TEST_F(FieldPerpTest, ConstIndexingWithInd3D) {
  FieldPerp field(0.0);
  field.setIndex(2);
  const BoutReal sentinel = 2.0;
  int ix = 1, iy = field.getIndex(), iz = 4;
  field(ix, iz) = sentinel;
  const Ind3D firstPoint(iz + nz * (iy + ny * ix), ny, nz);
  EXPECT_DOUBLE_EQ(field[firstPoint], sentinel);
  field[firstPoint] = -sentinel;
  EXPECT_DOUBLE_EQ(field(ix, iz), -sentinel);
#if CHECK > 2
  iy++;
  const Ind3D secondPoint(iz + nz * (iy + ny * ix), ny, nz);
  EXPECT_THROW(field[secondPoint], BoutException);
#endif
}

TEST_F(FieldPerpTest, IndexingToPointer) {
  FieldPerp field;

#if CHECK > 2
  // Data is empty
  EXPECT_THROW(field[0], BoutException);
#endif

  field.allocate();

#if CHECK > 2
  // Out of bounds
  EXPECT_THROW(field[-1], BoutException);
  EXPECT_NO_THROW(field[0]);
  EXPECT_THROW(field[nx], BoutException);
#endif

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      field[i][j] = i + j;
    }
  }

  EXPECT_DOUBLE_EQ(field(2, 2), 4);
}

TEST_F(FieldPerpTest, ConstIndexingToPointer) {
#if CHECK > 2
  const FieldPerp empty;
  EXPECT_THROW(empty[0], BoutException);
#endif

  const FieldPerp field = 2.0;

#if CHECK > 2
  // Out of bounds
  EXPECT_THROW(field[-1], BoutException);
  EXPECT_NO_THROW(field[0]);
  EXPECT_THROW(field[nx], BoutException);
#endif

  FieldPerp field2 = 0.0;

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      field2[i][j] = field[i][j];
    }
  }

  EXPECT_DOUBLE_EQ(field2(2, 2), 2.0);
}

#if CHECK > 2
TEST_F(FieldPerpTest, CheckNotEmpty) {
  FieldPerp field;

  EXPECT_THROW(field(0, 0), BoutException);
}

TEST_F(FieldPerpTest, BoundsCheck) {
  FieldPerp field;
  field = 1.0;

  EXPECT_THROW(field(-1, 0), BoutException);
  EXPECT_THROW(field(0, -1), BoutException);
  EXPECT_THROW(field(nx, 0), BoutException);
  EXPECT_THROW(field(0, nz), BoutException);
}

TEST_F(FieldPerpTest, ConstBoundsCheck) {
  const FieldPerp field = 1.0;

  EXPECT_THROW(field(-1, 0), BoutException);
  EXPECT_THROW(field(0, -1), BoutException);
  EXPECT_THROW(field(nx, 0), BoutException);
  EXPECT_THROW(field(0, nz), BoutException);
}

TEST_F(FieldPerpTest, ConstCheckNotEmpty) {
  const FieldPerp field;

  EXPECT_THROW(field(0, 0), BoutException);
}

TEST_F(FieldPerpTest, CheckData) {
  FieldPerp field;

  EXPECT_THROW(checkData(field), BoutException);

  field = 1.0;

  EXPECT_THROW(checkData(field), BoutException);

  field.setIndex(-10);

  EXPECT_THROW(checkData(field), BoutException);

  field.setIndex(100);

  EXPECT_THROW(checkData(field), BoutException);

  field.setIndex(0);

  EXPECT_NO_THROW(checkData(field));

  field(1, 1) = BoutNaN;

  EXPECT_THROW(checkData(field), BoutException);

  field = 1.0;
  field(0, 0) = BoutNaN;

  EXPECT_NO_THROW(checkData(field));
  EXPECT_NO_THROW(checkData(field, RGN_NOX));
  EXPECT_THROW(checkData(field, RGN_ALL), BoutException);
}

TEST_F(FieldPerpTest, InvalidateGuards) {
  FieldPerp field;
  field.allocate(); // Calls invalidateGuards
  field = 1.0;      // Sets everywhere including boundaries

  const int nmesh = nx * nz;

  int sum = 0;
  for (const auto &i : field) {
    field[i] = 0.0; // Reset field value
    sum++;
  }
  EXPECT_EQ(sum, nmesh); // Field operator= hasn't been broken by invalidateGuards

  // Count the number of non-boundary points
  sum = 0;
  for (const auto &i : field.getRegion(RGN_NOX)) {
    field[i] = 0.0; // Reset field value
    sum++;
  }
  const int nbndry = nmesh - sum;

  auto localmesh = field.getMesh();
  EXPECT_NO_THROW(checkData(field(0, 0)));
  EXPECT_NO_THROW(checkData(field(localmesh->xstart, 0)));

  invalidateGuards(field);

  EXPECT_THROW(checkData(field(0, 0)), BoutException);
  EXPECT_NO_THROW(checkData(field(localmesh->xstart, 0)));

  sum = 0;
  for (const auto &i : field) {
    if (!finite(field[i]))
      sum++;
  }
  EXPECT_EQ(sum, nbndry);
}

#endif // CHECK > 2

TEST_F(FieldPerpTest, CreateFromBoutReal) {
  FieldPerp field(1.0);

  EXPECT_TRUE(IsFieldEqual(field, 1.0));
}

TEST_F(FieldPerpTest, CreateFromFieldPerp) {
  FieldPerp field(99.0);
  FieldPerp result(field);

  EXPECT_TRUE(IsFieldEqual(result, 99.0));
}

TEST_F(FieldPerpTest, AssignFromBoutReal) {
  FieldPerp field;

  field = 2.0;

  EXPECT_TRUE(IsFieldEqual(field, 2.0));
}

TEST_F(FieldPerpTest, AssignFromInvalid) {
  FieldPerp field;

  EXPECT_NO_THROW(field = std::nan(""));
  EXPECT_TRUE(IsFieldEqual(field, std::nan("")));
}

TEST_F(FieldPerpTest, UnaryMinus) {
  FieldPerp field;
  field.setIndex(0);

  field = 2.0;
  field = -field;

  EXPECT_TRUE(IsFieldEqual(field, -2.0));
}

TEST_F(FieldPerpTest, AddEqualsBoutReal) {
  FieldPerp a;
  a.setIndex(0);

  a = 1.0;
  a += 5.0;

  EXPECT_TRUE(IsFieldEqual(a, 6.0));

  // Check case where field is not unique
  auto c = a;
  c += 5.0;

  EXPECT_TRUE(IsFieldEqual(a, 6.0));
  EXPECT_TRUE(IsFieldEqual(c, 11.0));
}

TEST_F(FieldPerpTest, AddEqualsFieldPerp) {
  FieldPerp a, b;
  a.setIndex(0);
  b.setIndex(0);

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

TEST_F(FieldPerpTest, AddEqualsField2D) {
  FieldPerp a;
  Field2D b;

  a = 2.0;
  b = 3.0;
#if CHECK > 2
  EXPECT_THROW(a += b, BoutException);
#endif
  a.setIndex(1);
  EXPECT_NO_THROW(a += b);

  EXPECT_TRUE(IsFieldEqual(a, 5.0));

  // Check case where field is not unique
  auto c = a;
  c += b;

  EXPECT_TRUE(IsFieldEqual(a, 5.0));
  EXPECT_TRUE(IsFieldEqual(c, 8.0));
}

TEST_F(FieldPerpTest, AddEqualsField3D) {
  FieldPerp a;
  Field3D b;

  a = 2.0;
  b = 3.0;
#if CHECK > 2
  EXPECT_THROW(a += b, BoutException);
#endif
  a.setIndex(1);
  EXPECT_NO_THROW(a += b);

  EXPECT_TRUE(IsFieldEqual(a, 5.0));

  // Check case where field is not unique
  auto c = a;
  c += b;

  EXPECT_TRUE(IsFieldEqual(a, 5.0));
  EXPECT_TRUE(IsFieldEqual(c, 8.0));
}

TEST_F(FieldPerpTest, AddFieldPerpBoutReal) {
  FieldPerp a, b;
  a.setIndex(0);

  a = 1.0;
  b = a + 2.0;

  EXPECT_TRUE(IsFieldEqual(b, 3.0));
}

TEST_F(FieldPerpTest, AddBoutRealFieldPerp) {
  FieldPerp a, b;
  a.setIndex(0);

  a = 1.0;
  b = 3.0 + a;

  EXPECT_TRUE(IsFieldEqual(b, 4.0));
}

TEST_F(FieldPerpTest, AddFieldPerpFieldPerp) {
  FieldPerp a, b, c;
  a.setIndex(0);
  b.setIndex(0);

  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsFieldEqual(c, 3.0));
}

TEST_F(FieldPerpTest, AddFieldPerpField2D) {
  FieldPerp a, c;
  Field2D b;

  a.setIndex(1);
  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsFieldEqual(c, 3.0));
}

TEST_F(FieldPerpTest, AddField2DFieldPerp) {
  FieldPerp b, c;
  Field2D a;

  b.setIndex(1);
  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsFieldEqual(c, 3.0));
}

TEST_F(FieldPerpTest, AddFieldPerpField3D) {
  FieldPerp a, c;
  Field3D b;

  a.setIndex(1);
  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsFieldEqual(c, 3.0));
}

TEST_F(FieldPerpTest, AddField3DFieldPerp) {
  FieldPerp b, c;
  Field3D a;

  b.setIndex(1);
  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsFieldEqual(c, 3.0));
}

TEST_F(FieldPerpTest, MultiplyEqualsBoutReal) {
  FieldPerp a;
  a.setIndex(0);

  a = 2.0;
  a *= 1.5;

  EXPECT_TRUE(IsFieldEqual(a, 3.0));

  // Check case where field is not unique
  auto c = a;
  c *= 1.5;

  EXPECT_TRUE(IsFieldEqual(a, 3.0));
  EXPECT_TRUE(IsFieldEqual(c, 4.5));
}

TEST_F(FieldPerpTest, MultiplyEqualsFieldPerp) {
  FieldPerp a, b;
  a.setIndex(0);
  b.setIndex(0);

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

TEST_F(FieldPerpTest, MultiplyEqualsField2D) {
  FieldPerp a;
  Field2D b;

  a = 2.5;
  b = 4.0;
#if CHECK > 2
  EXPECT_THROW(a *= b, BoutException);
#endif
  a.setIndex(1);
  EXPECT_NO_THROW(a *= b);

  EXPECT_TRUE(IsFieldEqual(a, 10.0));

  // Check case where field is not unique
  auto c = a;
  c *= b;

  EXPECT_TRUE(IsFieldEqual(a, 10.0));
  EXPECT_TRUE(IsFieldEqual(c, 40.0));
}

TEST_F(FieldPerpTest, MultiplyEqualsField3D) {
  FieldPerp a;
  Field3D b;

  a = 2.5;
  b = 4.0;
#if CHECK > 2
  EXPECT_THROW(a *= b, BoutException);
#endif
  a.setIndex(1);
  EXPECT_NO_THROW(a *= b);

  EXPECT_TRUE(IsFieldEqual(a, 10.0));

  // Check case where field is not unique
  auto c = a;
  c *= b;

  EXPECT_TRUE(IsFieldEqual(a, 10.0));
  EXPECT_TRUE(IsFieldEqual(c, 40.0));
}

TEST_F(FieldPerpTest, MultiplyFieldPerpBoutReal) {
  FieldPerp a, b;
  a.setIndex(0);

  a = 1.5;
  b = a * 2.0;

  EXPECT_TRUE(IsFieldEqual(b, 3.0));
}

TEST_F(FieldPerpTest, MultiplyBoutRealFieldPerp) {
  FieldPerp a, b;
  a.setIndex(0);

  a = 2.5;
  b = 3.0 * a;

  EXPECT_TRUE(IsFieldEqual(b, 7.5));
}

TEST_F(FieldPerpTest, MultiplyFieldPerpFieldPerp) {
  FieldPerp a, b, c;
  a.setIndex(0);
  b.setIndex(0);

  a = 4.0;
  b = 8.0;
  c = a * b;

  EXPECT_TRUE(IsFieldEqual(c, 32.0));
}

TEST_F(FieldPerpTest, MultiplyFieldPerpField2D) {
  FieldPerp a, c;
  Field2D b;

  a.setIndex(1);
  a = 4.0;
  b = 8.0;
  c = a * b;

  EXPECT_TRUE(IsFieldEqual(c, 32.0));
}

TEST_F(FieldPerpTest, MultiplyField2DFieldPerp) {
  FieldPerp b, c;
  Field2D a;

  b.setIndex(1);
  a = 4.0;
  b = 8.0;
  c = a * b;

  EXPECT_TRUE(IsFieldEqual(c, 32.0));
}

TEST_F(FieldPerpTest, MultiplyFieldPerpField3D) {
  FieldPerp a, c;
  Field3D b;

  a.setIndex(1);
  a = 4.0;
  b = 8.0;
  c = a * b;

  EXPECT_TRUE(IsFieldEqual(c, 32.0));
}

TEST_F(FieldPerpTest, MultiplyField3DFieldPerp) {
  FieldPerp b, c;
  Field3D a;

  b.setIndex(1);
  a = 4.0;
  b = 8.0;
  c = a * b;

  EXPECT_TRUE(IsFieldEqual(c, 32.0));
}

TEST_F(FieldPerpTest, SubtractEqualsBoutReal) {
  FieldPerp a;
  a.setIndex(0);

  a = 1.0;
  a -= 5.0;

  EXPECT_TRUE(IsFieldEqual(a, -4.0));

  // Check case where field is not unique
  auto c = a;
  c -= 5.0;

  EXPECT_TRUE(IsFieldEqual(a, -4.0));
  EXPECT_TRUE(IsFieldEqual(c, -9.0));
}

TEST_F(FieldPerpTest, SubtractEqualsFieldPerp) {
  FieldPerp a, b;
  a.setIndex(0);
  b.setIndex(0);

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

TEST_F(FieldPerpTest, SubtractEqualsField2D) {
  FieldPerp a;
  Field2D b;

  a = 2.0;
  b = 7.0;
#if CHECK > 2
  EXPECT_THROW(a -= b, BoutException);
#endif
  a.setIndex(1);
  EXPECT_NO_THROW(a -= b);

  EXPECT_TRUE(IsFieldEqual(a, -5.0));

  // Check case where field is not unique
  auto c = a;
  c -= b;

  EXPECT_TRUE(IsFieldEqual(a, -5.0));
  EXPECT_TRUE(IsFieldEqual(c, -12.0));
}

TEST_F(FieldPerpTest, SubtractEqualsField3D) {
  FieldPerp a;
  Field3D b;

  a = 2.0;
  b = 7.0;
#if CHECK > 2
  EXPECT_THROW(a -= b, BoutException);
#endif
  a.setIndex(1);
  EXPECT_NO_THROW(a -= b);

  EXPECT_TRUE(IsFieldEqual(a, -5.0));

  // Check case where field is not unique
  auto c = a;
  c -= b;

  EXPECT_TRUE(IsFieldEqual(a, -5.0));
  EXPECT_TRUE(IsFieldEqual(c, -12.0));
}

TEST_F(FieldPerpTest, SubtractFieldPerpBoutReal) {
  FieldPerp a, b;
  a.setIndex(0);

  a = 10.0;
  b = a - 2.0;

  EXPECT_TRUE(IsFieldEqual(b, 8.0));
}

TEST_F(FieldPerpTest, SubtractBoutRealFieldPerp) {
  FieldPerp a, b;
  a.setIndex(0);

  a = 10.0;
  b = 3.0 - a;

  EXPECT_TRUE(IsFieldEqual(b, -7.0));
}

TEST_F(FieldPerpTest, SubtractFieldPerpFieldPerp) {
  FieldPerp a, b, c;
  a.setIndex(0);
  b.setIndex(0);

  a = 10.0;
  b = 20.0;
  c = a - b;

  EXPECT_TRUE(IsFieldEqual(c, -10.0));
}

TEST_F(FieldPerpTest, SubtractFieldPerpField2D) {
  FieldPerp a, c;
  Field2D b;

  a.setIndex(1);
  a = 10.0;
  b = 20.0;
  c = a - b;

  EXPECT_TRUE(IsFieldEqual(c, -10.0));
}

TEST_F(FieldPerpTest, SubtractField2DFieldPerp) {
  FieldPerp b, c;
  Field2D a;

  b.setIndex(1);
  a = 10.0;
  b = 20.0;
  c = a - b;

  EXPECT_TRUE(IsFieldEqual(c, -10.0));
}

TEST_F(FieldPerpTest, SubtractFieldPerpField3D) {
  FieldPerp a, c;
  Field3D b;

  a.setIndex(1);
  a = 10.0;
  b = 20.0;
  c = a - b;

  EXPECT_TRUE(IsFieldEqual(c, -10.0));
}

TEST_F(FieldPerpTest, SubtractField3DFieldPerp) {
  FieldPerp b, c;
  Field3D a;

  b.setIndex(1);
  a = 10.0;
  b = 20.0;
  c = a - b;

  EXPECT_TRUE(IsFieldEqual(c, -10.0));
}

TEST_F(FieldPerpTest, DivideEqualsBoutReal) {
  FieldPerp a;
  a.setIndex(0);

  a = 2.5;
  a /= 5.0;

  EXPECT_TRUE(IsFieldEqual(a, 0.5));

  // Check case where field is not unique
  auto c = a;
  c /= 5.0;

  EXPECT_TRUE(IsFieldEqual(a, 0.5));
  EXPECT_TRUE(IsFieldEqual(c, 0.1));
}

TEST_F(FieldPerpTest, DivideEqualsFieldPerp) {
  FieldPerp a, b;
  a.setIndex(0);
  b.setIndex(0);

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

TEST_F(FieldPerpTest, DivideEqualsField2D) {
  FieldPerp a;
  Field2D b;

  a = 5.0;
  b = 2.5;
#if CHECK > 2
  EXPECT_THROW(a /= b, BoutException);
#endif
  a.setIndex(1);
  EXPECT_NO_THROW(a /= b);

  EXPECT_TRUE(IsFieldEqual(a, 2.0));

  // Check case where field is not unique
  auto c = a;
  c /= b;

  EXPECT_TRUE(IsFieldEqual(a, 2.0));
  EXPECT_TRUE(IsFieldEqual(c, 0.8));
}

TEST_F(FieldPerpTest, DivideEqualsField3D) {
  FieldPerp a;
  Field3D b;

  a = 5.0;
  b = 2.5;
#if CHECK > 2
  EXPECT_THROW(a /= b, BoutException);
#endif
  a.setIndex(1);
  EXPECT_NO_THROW(a /= b);

  EXPECT_TRUE(IsFieldEqual(a, 2.0));

  // Check case where field is not unique
  auto c = a;
  c /= b;

  EXPECT_TRUE(IsFieldEqual(a, 2.0));
  EXPECT_TRUE(IsFieldEqual(c, 0.8));
}

TEST_F(FieldPerpTest, DivideFieldPerpBoutReal) {
  FieldPerp a, b;
  a.setIndex(0);

  a = 3.0;
  b = a / 2.0;

  EXPECT_TRUE(IsFieldEqual(b, 1.5));
}

TEST_F(FieldPerpTest, DivideBoutRealFieldPerp) {
  FieldPerp a, b;
  a.setIndex(0);

  a = 2.5;
  b = 10.0 / a;

  EXPECT_TRUE(IsFieldEqual(b, 4.0));
}

TEST_F(FieldPerpTest, DivideFieldPerpFieldPerp) {
  FieldPerp a, b, c;
  a.setIndex(0);
  b.setIndex(0);

  a = 32.0;
  b = 8.0;
  c = a / b;

  EXPECT_TRUE(IsFieldEqual(c, 4.0));
}

TEST_F(FieldPerpTest, DivideFieldPerpField2D) {
  FieldPerp a, c;
  Field2D b;

  a.setIndex(1);
  a = 32.0;
  b = 8.0;
  c = a / b;

  EXPECT_TRUE(IsFieldEqual(c, 4.0));
}

TEST_F(FieldPerpTest, DivideField2DFieldPerp) {
  FieldPerp b, c;
  Field2D a;

  b.setIndex(1);
  a = 32.0;
  b = 8.0;
  c = a / b;

  EXPECT_TRUE(IsFieldEqual(c, 4.0));
}

TEST_F(FieldPerpTest, DivideFieldPerpField3D) {
  FieldPerp a, c;
  Field3D b;

  a.setIndex(1);
  a = 32.0;
  b = 8.0;
  c = a / b;

  EXPECT_TRUE(IsFieldEqual(c, 4.0));
}

TEST_F(FieldPerpTest, DivideField3DFieldPerp) {
  FieldPerp b, c;
  Field3D a;

  b.setIndex(1);
  a = 32.0;
  b = 8.0;
  c = a / b;

  EXPECT_TRUE(IsFieldEqual(c, 4.0));
}

TEST_F(FieldPerpTest, PowBoutRealFieldPerp) {
  FieldPerp a, b;
  a.setIndex(0);
  b.setIndex(0);

  a = 5.0;
  b = pow(2.0, a);

  EXPECT_TRUE(IsFieldEqual(b, 32.0));
}

TEST_F(FieldPerpTest, PowFieldPerpBoutReal) {
  FieldPerp a, b;
  a.setIndex(0);

  a = 5.0;
  b = pow(a, 2.0);

  EXPECT_TRUE(IsFieldEqual(b, 25.0));
}

TEST_F(FieldPerpTest, PowFieldPerpFieldPerp) {
  FieldPerp a, b, c;
  a.setIndex(0);
  b.setIndex(0);

  a = 2.0;
  b = 6.0;
  c = pow(a, b);

  EXPECT_TRUE(IsFieldEqual(c, 64.0));
}

TEST_F(FieldPerpTest, Sqrt) {
  FieldPerp field;
  field.setIndex(0);

  field = 16.0;
  EXPECT_TRUE(IsFieldEqual(sqrt(field), 4.0));
}

TEST_F(FieldPerpTest, Abs) {
  FieldPerp field;
  field.setIndex(0);

  field = -31.0;
  EXPECT_TRUE(IsFieldEqual(abs(field), 31.0));
}

TEST_F(FieldPerpTest, Exp) {
  FieldPerp field;
  field.setIndex(0);

  field = 2.5;
  const BoutReal expected = 12.182493960703473;
  EXPECT_TRUE(IsFieldEqual(exp(field), expected));
}

TEST_F(FieldPerpTest, Log) {
  FieldPerp field;
  field.setIndex(0);

  field = 12.182493960703473;
  const BoutReal expected = 2.5;
  EXPECT_TRUE(IsFieldEqual(log(field), expected));
}

TEST_F(FieldPerpTest, LogExp) {
  FieldPerp field;
  field.setIndex(0);

  field = 2.5;
  const BoutReal expected = 2.5;
  EXPECT_TRUE(IsFieldEqual(log(exp(field)), expected));
}

TEST_F(FieldPerpTest, Sin) {
  FieldPerp field;
  field.setIndex(0);

  field = PI / 2.0;
  EXPECT_TRUE(IsFieldEqual(sin(field), 1.0));

  field = PI;
  EXPECT_TRUE(IsFieldEqual(sin(field), 0.0));
}

TEST_F(FieldPerpTest, Cos) {
  FieldPerp field;
  field.setIndex(0);

  field = PI / 2.0;
  EXPECT_TRUE(IsFieldEqual(cos(field), 0.0));

  field = PI;
  EXPECT_TRUE(IsFieldEqual(cos(field), -1.0));
}

TEST_F(FieldPerpTest, Tan) {
  FieldPerp field;
  field.setIndex(0);

  field = PI / 4.0;
  EXPECT_TRUE(IsFieldEqual(tan(field), 1.0));

  field = PI;
  EXPECT_TRUE(IsFieldEqual(tan(field), 0.0));
}

TEST_F(FieldPerpTest, Sinh) {
  FieldPerp field;
  field.setIndex(0);

  field = 1.0;
  const BoutReal expected = 1.1752011936438014;
  EXPECT_TRUE(IsFieldEqual(sinh(field), expected));

  field = -1.0;
  EXPECT_TRUE(IsFieldEqual(sinh(field), -expected));
}

TEST_F(FieldPerpTest, Cosh) {
  FieldPerp field;
  field.setIndex(0);

  field = 1.0;
  const BoutReal expected = 1.5430806348152437;
  EXPECT_TRUE(IsFieldEqual(cosh(field), expected));

  field = -1.0;
  EXPECT_TRUE(IsFieldEqual(cosh(field), expected));
}

TEST_F(FieldPerpTest, Tanh) {
  FieldPerp field;
  field.setIndex(0);

  field = 1.0;
  const BoutReal expected = 0.761594155955764;
  EXPECT_TRUE(IsFieldEqual(tanh(field), expected));

  field = -1.0;
  EXPECT_TRUE(IsFieldEqual(tanh(field), -expected));
}

TEST_F(FieldPerpTest, Floor) {
  FieldPerp field;
  field.setIndex(0);

  field = 50.0;
  field(1, 1) = 49.9;
  field(2, 3) = -20;

  const BoutReal floor_value = 50.0;

  EXPECT_TRUE(IsFieldEqual(floor(field, floor_value), floor_value));
}

TEST_F(FieldPerpTest, Min) {
  FieldPerp field;
  field.setIndex(0);

  field = 50.0;
  field(0, 0) = -99.0;
  field(1, 1) = 60.0;
  field(1, 2) = 40.0;
  field(2, 4) = 99.0;

  // min doesn't include guard cells
  const BoutReal min_value = 40.0;

  EXPECT_EQ(min(field, false), min_value);
  EXPECT_EQ(min(field, false, RGN_ALL), -99.0);
  EXPECT_EQ(min(field, true, RGN_ALL), -99.0);
}

TEST_F(FieldPerpTest, Max) {
  FieldPerp field;
  field.setIndex(0);

  field = 50.0;
  field(0, 0) = -99.0;
  field(1, 1) = 40.0;
  field(1, 2) = 60.0;
  field(2, 4) = 99.0;

  // max doesn't include guard cells
  const BoutReal max_value = 60.0;

  EXPECT_EQ(max(field, false), max_value);
  EXPECT_EQ(max(field, false, RGN_ALL), 99.0);
  EXPECT_EQ(max(field, true, RGN_ALL), 99.0);
}

TEST_F(FieldPerpTest, OperatorEqualsFieldPerp) {
  FieldPerp field;

  // Create field with non-default arguments so we can check they get copied
  // to 'field'.
  // Note that Average z-direction type is not really allowed for FieldPerp, but
  // we don't check anywhere at the moment.
  FieldPerp field2{mesh_staggered, CELL_XLOW, 2, {YDirectionType::Aligned, ZDirectionType::Average}};

  field = field2;

  EXPECT_TRUE(areFieldsCompatible(field, field2));
  EXPECT_EQ(field.getMesh(), field2.getMesh());
  EXPECT_EQ(field.getLocation(), field2.getLocation());
  EXPECT_EQ(field.getDirectionY(), field2.getDirectionY());
  EXPECT_EQ(field.getDirectionZ(), field2.getDirectionZ());
}

TEST_F(FieldPerpTest, EmptyFrom) {
  // Create field with non-default arguments so we can check they get copied
  // to 'field2'.
  // Note that Average z-direction type is not really allowed for FieldPerp, but
  // we don't check anywhere at the moment.
  FieldPerp field{mesh_staggered, CELL_XLOW, 3, {YDirectionType::Aligned, ZDirectionType::Average}};
  field = 5.;

  FieldPerp field2{emptyFrom(field)};
  EXPECT_EQ(field2.getMesh(), mesh_staggered);
  EXPECT_EQ(field2.getLocation(), CELL_XLOW);
  EXPECT_EQ(field2.getIndex(), 3);
  EXPECT_EQ(field2.getDirectionY(), YDirectionType::Aligned);
  EXPECT_EQ(field2.getDirectionZ(), ZDirectionType::Average);
  EXPECT_TRUE(field2.isAllocated());
}

TEST_F(FieldPerpTest, ZeroFrom) {
  // Create field with non-default arguments so we can check they get copied
  // to 'field2'.
  // Note that Average z-direction type is not really allowed for FieldPerp, but
  // we don't check anywhere at the moment.
  FieldPerp field{mesh_staggered, CELL_XLOW, 3, {YDirectionType::Aligned, ZDirectionType::Average}};
  field = 5.;

  FieldPerp field2{zeroFrom(field)};
  EXPECT_EQ(field2.getMesh(), mesh_staggered);
  EXPECT_EQ(field2.getLocation(), CELL_XLOW);
  EXPECT_EQ(field2.getIndex(), 3);
  EXPECT_EQ(field2.getDirectionY(), YDirectionType::Aligned);
  EXPECT_EQ(field2.getDirectionZ(), ZDirectionType::Average);
  EXPECT_TRUE(field2.isAllocated());
  EXPECT_TRUE(IsFieldEqual(field2, 0.));
}
#pragma GCC diagnostic pop
