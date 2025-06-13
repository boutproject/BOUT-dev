#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "bout/facefield3d.hxx"
#include "bout/mesh_facefield_comm.hxx"
#include "../fake_mesh.hxx"

/// Global variables for testing
static FakeMesh* test_mesh = nullptr;

class FaceFieldCommTest : public ::testing::Test {
public:
  FaceFieldCommTest() { 
    // Create a mesh with ghost cells
    test_mesh = new FakeMesh(5, 5, 5);
    test_mesh->xstart = 2;
    test_mesh->xend = 4;
    test_mesh->ystart = 2;
    test_mesh->yend = 4;
  }
  ~FaceFieldCommTest() { 
    delete test_mesh; 
    test_mesh = nullptr;
  }
};

TEST_F(FaceFieldCommTest, BasicCommunication) {
  // Test that we can call communicate on a FaceField3D
  FaceField3D field(test_mesh);
  field = 1.0;
  
  // This should compile and not throw
  EXPECT_NO_THROW(communicate(test_mesh, field));
}

TEST_F(FaceFieldCommTest, CommunicateXZ) {
  // Test XZ-only communication
  FaceField3D field(test_mesh);
  field = 2.0;
  
  // This should compile and not throw
  EXPECT_NO_THROW(communicateXZ(test_mesh, field));
}

TEST_F(FaceFieldCommTest, CommunicateYZ) {
  // Test YZ-only communication
  FaceField3D field(test_mesh);
  field = 3.0;
  
  // This should compile and not throw
  EXPECT_NO_THROW(communicateYZ(test_mesh, field));
}

TEST_F(FaceFieldCommTest, ConvenienceMethods) {
  // Test convenience methods that get mesh from field
  FaceField3D field(test_mesh);
  field = 4.0;
  
  // These should compile and not throw
  EXPECT_NO_THROW(communicate(field));
  EXPECT_NO_THROW(communicateXZ(field));
  EXPECT_NO_THROW(communicateYZ(field));
}

TEST_F(FaceFieldCommTest, ComponentLocations) {
  // Verify that components maintain their staggered locations
  // after communication
  FaceField3D field(test_mesh);
  field = 5.0;
  
  // Store original locations
  CELL_LOC x_loc = field.x().getLocation();
  CELL_LOC y_loc = field.y().getLocation();
  CELL_LOC z_loc = field.z().getLocation();
  
  // Communicate
  communicate(field);
  
  // Verify locations unchanged
  EXPECT_EQ(field.x().getLocation(), x_loc);
  EXPECT_EQ(field.y().getLocation(), y_loc);
  EXPECT_EQ(field.z().getLocation(), z_loc);
  
  // Verify they are the expected staggered locations
  EXPECT_EQ(field.x().getLocation(), CELL_XLOW);
  EXPECT_EQ(field.y().getLocation(), CELL_YLOW);
  EXPECT_EQ(field.z().getLocation(), CELL_ZLOW);
}

// Note: More comprehensive multi-processor tests would require
// actual MPI setup and are better suited for integrated tests
// rather than unit tests.