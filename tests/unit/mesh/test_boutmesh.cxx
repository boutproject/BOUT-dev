#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "bout/mesh.hxx"
#include "output.hxx"
#include "unused.hxx"

class FakeGridDataSource : public GridDataSource {
public:
  FakeGridDataSource(){};
  ~FakeGridDataSource(){};
  bool hasVar(const string &UNUSED(name)) { return false; };
  bool get(Mesh *UNUSED(m), int &UNUSED(ival), const string &UNUSED(name)) {
    return true;
  };
  bool get(Mesh *UNUSED(m), BoutReal &UNUSED(rval), const string &UNUSED(name)) {
    return true;
  }
  bool get(Mesh *UNUSED(m), Field2D &UNUSED(var), const string &UNUSED(name),
           BoutReal UNUSED(def) = 0.0) {
    return true;
  }
  bool get(Mesh *UNUSED(m), Field3D &UNUSED(var), const string &UNUSED(name),
           BoutReal UNUSED(def) = 0.0) {
    return true;
  }
  bool get(Mesh *UNUSED(m), vector<int> &UNUSED(var), const string &UNUSED(name),
           int UNUSED(len), int UNUSED(offset) = 0,
           Direction UNUSED(dir) = Direction::X) {
    return true;
  }
  bool get(Mesh *UNUSED(m), vector<BoutReal> &UNUSED(var), const string &UNUSED(name),
           int UNUSED(len), int UNUSED(offset) = 0,
           Direction UNUSED(dir) = Direction::X) {
    return true;
  }
};

TEST(BoutMeshTest, NullOptionsCheck) {
  // Temporarily turn off outputs to make test quiet
  output_info.disable();
  output_warn.disable();
  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
  output_info.enable();
  output_warn.enable();
}
