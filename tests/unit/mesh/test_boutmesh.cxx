#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "bout/mesh.hxx"
#include "output.hxx"
#include "unused.hxx"
#include <unordered_map>

class FakeGridDataSource : public GridDataSource {
public:
  FakeGridDataSource(){};
  ~FakeGridDataSource(){};
  bool hasVar(const std::string &UNUSED(name)) { return false; };
  bool get(Mesh *UNUSED(m), int &ival, const std::string &name) {
    if (intvars.count(name)>0) {
      ival = intvars.at(name);
      return true;
    }
    return false;
  };
  bool get(Mesh *UNUSED(m), BoutReal &UNUSED(rval), const std::string &UNUSED(name)) {
    return false;
  }
  bool get(Mesh *UNUSED(m), Field2D &UNUSED(var), const std::string &UNUSED(name),
           BoutReal UNUSED(def) = 0.0) {
    return false;
  }
  bool get(Mesh *UNUSED(m), Field3D &UNUSED(var), const std::string &UNUSED(name),
           BoutReal UNUSED(def) = 0.0) {
    return false;
  }
  bool get(Mesh *UNUSED(m), std::vector<int> &UNUSED(var), const std::string &UNUSED(name),
           int UNUSED(len), int UNUSED(offset) = 0,
           Direction UNUSED(dir) = GridDataSource::X) {
    return false;
  }
  bool get(Mesh *UNUSED(m), std::vector<BoutReal> &UNUSED(var), const std::string &UNUSED(name),
           int UNUSED(len), int UNUSED(offset) = 0,
           Direction UNUSED(dir) = GridDataSource::X) {
    return false;
  }

  std::unordered_map<std::string, int> intvars { {"nx", 6}, {"ny", 7 }, {"nz", 5}};
};

/// Test fixture with a BoutMesh
class BoutMeshTest : public ::testing::Test {
protected:
  static void SetUpTestCase() {
    output_info.disable();
    output_warn.disable();
    output_progress.disable();
  }

  static void TearDownTestCase() {
    output_info.enable();
    output_warn.enable();
    output_progress.enable();
  }

public:
  BoutMeshTest() : source(), localmesh(Mesh::create(&source)) {
    localmesh->StaggerGrids = true;
    output_info.disable();
    localmesh->load();
  }
  FakeGridDataSource source;
  Mesh* localmesh;
};

TEST_F(BoutMeshTest, NullOptionsCheck) {
  // Temporarily turn off outputs to make test quiet
  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
}

TEST_F(BoutMeshTest, AddCoordinatesToMeshCENTRE) {
  EXPECT_NO_THROW(localmesh->addCoordinates(CELL_CENTRE));
  EXPECT_NO_THROW(localmesh->getCoordinates(CELL_CENTRE)->geometry());
  EXPECT_NO_THROW(localmesh->addCoordinates(CELL_YLOW));
  EXPECT_THROW(localmesh->getCoordinates(CELL_CENTRE)->geometry(), BoutException);
}

TEST_F(BoutMeshTest, AddCoordinatesToMeshXLOW) {
  EXPECT_NO_THROW(localmesh->addCoordinates(CELL_XLOW));
  EXPECT_THROW(localmesh->addCoordinates(CELL_XLOW), BoutException);
  EXPECT_THROW(localmesh->getCoordinates(CELL_XLOW)->geometry(), BoutException);
}

TEST_F(BoutMeshTest, AddCoordinatesToMeshYLOW) {
  EXPECT_NO_THROW(localmesh->addCoordinates(CELL_YLOW));
  EXPECT_THROW(localmesh->addCoordinates(CELL_YLOW), BoutException);
  EXPECT_THROW(localmesh->getCoordinates(CELL_YLOW)->geometry(), BoutException);
}

TEST_F(BoutMeshTest, AddCoordinatesToMeshZLOW) {
  EXPECT_NO_THROW(localmesh->addCoordinates(CELL_ZLOW));
  EXPECT_THROW(localmesh->addCoordinates(CELL_ZLOW), BoutException);
  EXPECT_THROW(localmesh->getCoordinates(CELL_ZLOW)->geometry(), BoutException);
}
