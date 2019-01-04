#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "bout_types.hxx"
#include "output.hxx"
#include "unused.hxx"
#include "bout/griddata.hxx"
#include "bout/mesh.hxx"

class FakeGridDataSource : public GridDataSource {
public:
  FakeGridDataSource(){};
  ~FakeGridDataSource(){};
  bool hasVar(const std::string &UNUSED(name)) { return false; };
  bool get(Mesh *UNUSED(m), int &UNUSED(ival), const std::string &UNUSED(name)) {
    return true;
  };
  bool get(Mesh *UNUSED(m), BoutReal &UNUSED(rval), const std::string &UNUSED(name)) {
    return true;
  }
  bool get(Mesh *UNUSED(m), Field2D &UNUSED(var), const std::string &UNUSED(name),
           BoutReal UNUSED(def) = 0.0) {
    return true;
  }
  bool get(Mesh *UNUSED(m), Field3D &UNUSED(var), const std::string &UNUSED(name),
           BoutReal UNUSED(def) = 0.0) {
    return true;
  }
  bool get(Mesh *UNUSED(m), std::vector<int> &UNUSED(var), const std::string &UNUSED(name),
           int UNUSED(len), int UNUSED(offset) = 0,
           Direction UNUSED(dir) = GridDataSource::X) {
    return true;
  }
  bool get(Mesh *UNUSED(m), std::vector<BoutReal> &UNUSED(var), const std::string &UNUSED(name),
           int UNUSED(len), int UNUSED(offset) = 0,
           Direction UNUSED(dir) = GridDataSource::X) {
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
