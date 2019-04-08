#include "gtest/gtest.h"

#include "field.hxx"
#include "test_extras.hxx"
#include "where.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

namespace {
constexpr auto le0 = 4.0;
constexpr auto gt0 = 2.0;
} // namespace

// Reuse the "standard" fixture for FakeMesh
class WhereTest : public FakeMeshFixture {
public:
  WhereTest()
      : FakeMeshFixture(), expected_2d2d(mesh), expected_2d3d(mesh), test_2d(mesh) {
    fillField(expected_2d2d, {{le0, le0, le0, gt0, gt0},
                              {gt0, le0, le0, le0, le0},
                              {le0, le0, gt0, gt0, gt0}});

    fillField(expected_2d3d, {{{le0, le0, le0, le0, le0, le0, le0},
                               {le0, le0, le0, le0, le0, le0, le0},
                               {le0, le0, le0, le0, le0, le0, le0},
                               {gt0, gt0, gt0, gt0, gt0, gt0, gt0},
                               {gt0, gt0, gt0, gt0, gt0, gt0, gt0}},

                              {{gt0, gt0, gt0, gt0, gt0, gt0, gt0},
                               {le0, le0, le0, le0, le0, le0, le0},
                               {le0, le0, le0, le0, le0, le0, le0},
                               {le0, le0, le0, le0, le0, le0, le0},
                               {le0, le0, le0, le0, le0, le0, le0}},

                              {{le0, le0, le0, le0, le0, le0, le0},
                               {le0, le0, le0, le0, le0, le0, le0},
                               {gt0, gt0, gt0, gt0, gt0, gt0, gt0},
                               {gt0, gt0, gt0, gt0, gt0, gt0, gt0},
                               {gt0, gt0, gt0, gt0, gt0, gt0, gt0}}});

    fillField(test_2d,
              {{-2., -1., 0., 1., 2.}, {1., 0., -1., -2., -3.}, {-1., 0., 1., 2., 3.}});
  }

  virtual ~WhereTest() = default;

  Field2D expected_2d2d;
  Field3D expected_2d3d;

  Field2D test_2d;
};

TEST_F(WhereTest, Field2DField3DField3D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field3D{gt0}, Field3D{le0}), expected_2d3d));
}

TEST_F(WhereTest, Field2DField3DBoutReal) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field3D{gt0}, le0), expected_2d3d));
}

TEST_F(WhereTest, Field2DBoutRealField3D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, gt0, Field3D{le0}), expected_2d3d));
}

TEST_F(WhereTest, Field2DField3DField2D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field3D{gt0}, Field2D{le0}), expected_2d3d));
}

TEST_F(WhereTest, Field2DField2DField3D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field2D{gt0}, Field3D{le0}), expected_2d3d));
}

TEST_F(WhereTest, Field2DField2DField2D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field2D{gt0}, Field2D{le0}), expected_2d2d));
}

TEST_F(WhereTest, Field2DField2DBoutReal) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, Field2D{gt0}, le0), expected_2d2d));
}

TEST_F(WhereTest, Field2DBoutRealField2D) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, gt0, Field2D{le0}), expected_2d2d));
}

TEST_F(WhereTest, Field2DBoutRealBoutReal) {
  EXPECT_TRUE(IsFieldEqual(where(test_2d, gt0, le0), expected_2d2d));
}

TEST_F(WhereTest, Field3DBoutRealField3D) {
  auto test_3d = makeField<Field3D>([](Field3D::ind_type i) {
    return (i.y() - (ny / 2)) + (i.x() - (nx / 2)) + (i.z() - (nz / 2));
  });

  Field3D expected_3d3d;
  fillField(expected_3d3d, {{{le0, le0, le0, le0, le0, le0, le0},
                             {le0, le0, le0, le0, le0, le0, gt0},
                             {le0, le0, le0, le0, le0, gt0, gt0},
                             {le0, le0, le0, le0, gt0, gt0, gt0},
                             {le0, le0, le0, gt0, gt0, gt0, gt0}},
                            {{le0, le0, le0, le0, le0, le0, gt0},
                             {le0, le0, le0, le0, le0, gt0, gt0},
                             {le0, le0, le0, le0, gt0, gt0, gt0},
                             {le0, le0, le0, gt0, gt0, gt0, gt0},
                             {le0, le0, gt0, gt0, gt0, gt0, gt0}},
                            {{le0, le0, le0, le0, le0, gt0, gt0},
                             {le0, le0, le0, le0, gt0, gt0, gt0},
                             {le0, le0, le0, gt0, gt0, gt0, gt0},
                             {le0, le0, gt0, gt0, gt0, gt0, gt0},
                             {le0, gt0, gt0, gt0, gt0, gt0, gt0}}});

  EXPECT_TRUE(IsFieldEqual(where(test_3d, gt0, Field3D{le0}), expected_3d3d));
}
