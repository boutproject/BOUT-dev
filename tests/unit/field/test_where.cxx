#include "gtest/gtest.h"

#include "field.hxx"
#include "where.hxx"
#include "test_extras.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
class WhereTest : public FakeMeshFixture {
public:
  WhereTest() : FakeMeshFixture(), expected_2d(mesh), expected_3d(mesh) {
    fillField(expected_2d,
              {{4., 4., 4., 2., 2.}, {4., 4., 4., 2., 2.}, {4., 4., 4., 2., 2.}});

    fillField(expected_3d, {{{4., 4., 4., 4., 4., 4., 4.},
                             {4., 4., 4., 4., 4., 4., 4.},
                             {4., 4., 4., 4., 4., 4., 4.},
                             {2., 2., 2., 2., 2., 2., 2.},
                             {2., 2., 2., 2., 2., 2., 2.}},

                            {{4., 4., 4., 4., 4., 4., 4.},
                             {4., 4., 4., 4., 4., 4., 4.},
                             {4., 4., 4., 4., 4., 4., 4.},
                             {2., 2., 2., 2., 2., 2., 2.},
                             {2., 2., 2., 2., 2., 2., 2.}},

                            {{4., 4., 4., 4., 4., 4., 4.},
                             {4., 4., 4., 4., 4., 4., 4.},
                             {4., 4., 4., 4., 4., 4., 4.},
                             {2., 2., 2., 2., 2., 2., 2.},
                             {2., 2., 2., 2., 2., 2., 2.}}});

    test_2d = makeField<Field2D>(
        [](Field2D::ind_type i) { return i.y() - (WhereTest::ny / 2); });
    test_3d = makeField<Field3D>(
        [](Field3D::ind_type i) { return i.y() - (WhereTest::ny / 2); });
  }

  virtual ~WhereTest() = default;

  Field2D expected_2d;
  Field3D expected_3d;

  Field2D test_2d;
  Field3D test_3d;
};

TEST_F(WhereTest, Field2DField3DField3D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field3D{2.0}, Field3D{4.0}), expected_3d));
}

TEST_F(WhereTest, Field2DField3DBoutReal) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field3D{2.0}, 4.0), expected_3d));
}

TEST_F(WhereTest, Field2DBoutRealField3D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, 2.0, Field3D{4.0}), expected_3d));
}

TEST_F(WhereTest, Field2DField3DField2D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field3D{2.0}, Field2D{4.0}), expected_3d));
}

TEST_F(WhereTest, Field2DField2DField3D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field2D{2.0}, Field3D{4.0}), expected_3d));
}

TEST_F(WhereTest, Field2DField2DField2D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field2D{2.0}, Field2D{4.0}), expected_2d));
}

TEST_F(WhereTest, Field2DField2DBoutReal) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field2D{2.0}, 4.0), expected_2d));
}

TEST_F(WhereTest, Field2DBoutRealField2D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, 2.0, Field2D{4.0}), expected_2d));
}

TEST_F(WhereTest, Field2DBoutRealBoutReal) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, 2.0, 4.0), expected_2d));
}

TEST_F(WhereTest, Field3DBoutRealField3D) {
  EXPECT_TRUE(IsFieldEqual(where(test_3d, 2.0, Field3D{4.0}), expected_3d));
}
