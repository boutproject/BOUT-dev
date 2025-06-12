#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "bout/facefield3d.hxx"
#include "bout/mesh.hxx"

/// Global variables for testing
static FakeMesh* test_mesh = nullptr;

class FaceField3DTest : public ::testing::Test {
public:
  FaceField3DTest() { test_mesh = new FakeMesh(3, 3, 3); }
  ~FaceField3DTest() { 
    delete test_mesh; 
    test_mesh = nullptr;
  }
};

TEST_F(FaceField3DTest, Constructor) {
  FaceField3D field(test_mesh);
  
  // Check that components have correct locations
  EXPECT_EQ(field.x().getLocation(), CELL_XLOW);
  EXPECT_EQ(field.y().getLocation(), CELL_YLOW);
  EXPECT_EQ(field.z().getLocation(), CELL_ZLOW);
  
  // Check mesh is set correctly
  EXPECT_EQ(field.getMesh(), test_mesh);
}

TEST_F(FaceField3DTest, AssignmentScalar) {
  FaceField3D field(test_mesh);
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
  FaceField3D field1(test_mesh);
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
  FaceField3D field1(test_mesh);
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
  FaceField3D field1(test_mesh);
  FaceField3D field2(test_mesh);
  
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
  FaceField3D field1(test_mesh);
  FaceField3D field2(test_mesh);
  
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
  FaceField3D field1(test_mesh);
  FaceField3D field2(test_mesh);
  
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
  FaceField3D field(test_mesh);
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
  FaceField3D field(test_mesh);
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
  FaceField3D field(test_mesh);
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
  FaceField3D field1(test_mesh);
  Field3D field2(test_mesh);
  
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
  FaceField3D field(test_mesh);
  
  // Get time derivative
  FieldData* deriv = field.timeDeriv();
  
  // Check it's a FaceField3D
  FaceField3D* face_deriv = dynamic_cast<FaceField3D*>(deriv);
  EXPECT_NE(face_deriv, nullptr);
  
  // Check component derivatives are set up correctly
  EXPECT_EQ(field.x().timeDeriv(), &(face_deriv->x()));
  EXPECT_EQ(field.y().timeDeriv(), &(face_deriv->y()));
  EXPECT_EQ(field.z().timeDeriv(), &(face_deriv->z()));
}