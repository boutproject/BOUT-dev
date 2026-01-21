#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "bout/facefield3d.hxx"
#include "bout/mesh.hxx"
#include "../fake_mesh_fixture.hxx"

/// Test fixture for FaceField3D
class FaceField3DTest : public FakeMeshFixture {
public:
  FaceField3DTest() : FakeMeshFixture() {
    // Enable staggered grids for face fields
    bout::globals::mesh->StaggerGrids = true;
  }
};

TEST_F(FaceField3DTest, Constructor) {
  FaceField3D field(bout::globals::mesh);
  
  // Check that components have correct locations
  EXPECT_EQ(field.x().getLocation(), CELL_XLOW);
  EXPECT_EQ(field.y().getLocation(), CELL_YLOW);
  EXPECT_EQ(field.z().getLocation(), CELL_ZLOW);
  
  // Check mesh is set correctly
  EXPECT_EQ(field.getMesh(), bout::globals::mesh);
}

TEST_F(FaceField3DTest, AssignmentScalar) {
  FaceField3D field(bout::globals::mesh);
  field = 3.14;
  
  // Check all components are set
  for (int i = 0; i < field.x().getNx(); ++i) {
    for (int j = 0; j < field.x().getNy(); ++j) {
      for (int k = 0; k < field.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field.x()(i, j, k), 3.14);
        EXPECT_DOUBLE_EQ(field.y()(i, j, k), 3.14);
        EXPECT_DOUBLE_EQ(field.z()(i, j, k), 3.14);
      }
    }
  }
}

TEST_F(FaceField3DTest, CopyConstructor) {
  FaceField3D field1(bout::globals::mesh);
  field1 = 2.0;
  
  FaceField3D field2(field1);
  
  // Check values are copied
  for (int i = 0; i < field2.x().getNx(); ++i) {
    for (int j = 0; j < field2.x().getNy(); ++j) {
      for (int k = 0; k < field2.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field2.x()(i, j, k), 2.0);
        EXPECT_DOUBLE_EQ(field2.y()(i, j, k), 2.0);
        EXPECT_DOUBLE_EQ(field2.z()(i, j, k), 2.0);
      }
    }
  }
}

TEST_F(FaceField3DTest, MoveConstructor) {
  FaceField3D field1(bout::globals::mesh);
  field1 = 5.0;
  
  FaceField3D field2(std::move(field1));
  
  // Check values are moved
  for (int i = 0; i < field2.x().getNx(); ++i) {
    for (int j = 0; j < field2.x().getNy(); ++j) {
      for (int k = 0; k < field2.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field2.x()(i, j, k), 5.0);
        EXPECT_DOUBLE_EQ(field2.y()(i, j, k), 5.0);
        EXPECT_DOUBLE_EQ(field2.z()(i, j, k), 5.0);
      }
    }
  }
}

TEST_F(FaceField3DTest, Addition) {
  FaceField3D field1(bout::globals::mesh);
  FaceField3D field2(bout::globals::mesh);
  
  field1 = 1.0;
  field2 = 2.0;
  
  FaceField3D field3 = field1 + field2;
  
  // Check addition
  for (int i = 0; i < field3.x().getNx(); ++i) {
    for (int j = 0; j < field3.x().getNy(); ++j) {
      for (int k = 0; k < field3.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field3.x()(i, j, k), 3.0);
        EXPECT_DOUBLE_EQ(field3.y()(i, j, k), 3.0);
        EXPECT_DOUBLE_EQ(field3.z()(i, j, k), 3.0);
      }
    }
  }
}

TEST_F(FaceField3DTest, AdditionAssignment) {
  FaceField3D field1(bout::globals::mesh);
  FaceField3D field2(bout::globals::mesh);
  
  field1 = 1.0;
  field2 = 2.0;
  
  field1 += field2;
  
  // Check addition assignment
  for (int i = 0; i < field1.x().getNx(); ++i) {
    for (int j = 0; j < field1.x().getNy(); ++j) {
      for (int k = 0; k < field1.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field1.x()(i, j, k), 3.0);
        EXPECT_DOUBLE_EQ(field1.y()(i, j, k), 3.0);
        EXPECT_DOUBLE_EQ(field1.z()(i, j, k), 3.0);
      }
    }
  }
}

TEST_F(FaceField3DTest, Subtraction) {
  FaceField3D field1(bout::globals::mesh);
  FaceField3D field2(bout::globals::mesh);
  
  field1 = 5.0;
  field2 = 2.0;
  
  FaceField3D field3 = field1 - field2;
  
  // Check subtraction
  for (int i = 0; i < field3.x().getNx(); ++i) {
    for (int j = 0; j < field3.x().getNy(); ++j) {
      for (int k = 0; k < field3.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field3.x()(i, j, k), 3.0);
        EXPECT_DOUBLE_EQ(field3.y()(i, j, k), 3.0);
        EXPECT_DOUBLE_EQ(field3.z()(i, j, k), 3.0);
      }
    }
  }
}

TEST_F(FaceField3DTest, MultiplicationScalar) {
  FaceField3D field(bout::globals::mesh);
  field = 2.0;
  
  FaceField3D field2 = field * 3.0;
  
  // Check scalar multiplication
  for (int i = 0; i < field2.x().getNx(); ++i) {
    for (int j = 0; j < field2.x().getNy(); ++j) {
      for (int k = 0; k < field2.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field2.x()(i, j, k), 6.0);
        EXPECT_DOUBLE_EQ(field2.y()(i, j, k), 6.0);
        EXPECT_DOUBLE_EQ(field2.z()(i, j, k), 6.0);
      }
    }
  }
}

TEST_F(FaceField3DTest, DivisionScalar) {
  FaceField3D field(bout::globals::mesh);
  field = 6.0;
  
  FaceField3D field2 = field / 2.0;
  
  // Check scalar division
  for (int i = 0; i < field2.x().getNx(); ++i) {
    for (int j = 0; j < field2.x().getNy(); ++j) {
      for (int k = 0; k < field2.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field2.x()(i, j, k), 3.0);
        EXPECT_DOUBLE_EQ(field2.y()(i, j, k), 3.0);
        EXPECT_DOUBLE_EQ(field2.z()(i, j, k), 3.0);
      }
    }
  }
}

TEST_F(FaceField3DTest, UnaryNegation) {
  FaceField3D field(bout::globals::mesh);
  field = 2.0;
  
  FaceField3D field2 = -field;
  
  // Check negation
  for (int i = 0; i < field2.x().getNx(); ++i) {
    for (int j = 0; j < field2.x().getNy(); ++j) {
      for (int k = 0; k < field2.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field2.x()(i, j, k), -2.0);
        EXPECT_DOUBLE_EQ(field2.y()(i, j, k), -2.0);
        EXPECT_DOUBLE_EQ(field2.z()(i, j, k), -2.0);
      }
    }
  }
}

TEST_F(FaceField3DTest, MultiplicationField3D) {
  FaceField3D field1(bout::globals::mesh);
  Field3D field2(bout::globals::mesh);
  
  field1 = 2.0;
  field2 = 3.0;
  
  FaceField3D field3 = field1 * field2;
  
  // Check field multiplication
  for (int i = 0; i < field3.x().getNx(); ++i) {
    for (int j = 0; j < field3.x().getNy(); ++j) {
      for (int k = 0; k < field3.x().getNz(); ++k) {
        EXPECT_DOUBLE_EQ(field3.x()(i, j, k), 6.0);
        EXPECT_DOUBLE_EQ(field3.y()(i, j, k), 6.0);
        EXPECT_DOUBLE_EQ(field3.z()(i, j, k), 6.0);
      }
    }
  }
}

TEST_F(FaceField3DTest, TimeDeriv) {
  FaceField3D field(bout::globals::mesh);
  
  // Get time derivative
  FieldData* deriv = field.timeDeriv();
  
  // Check it's a FaceField3D
  FaceField3D* face_deriv = dynamic_cast<FaceField3D*>(deriv);
  EXPECT_NE(face_deriv, nullptr);
  
  // Check the derivative was created
  EXPECT_EQ(face_deriv->getMesh(), bout::globals::mesh);
}