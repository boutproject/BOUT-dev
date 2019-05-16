#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "bout/mesh.hxx"
#include "output.hxx"
#include "unused.hxx"

class FakeGridDataSource : public GridDataSource {
public:
  bool hasVar(const std::string& UNUSED(name)) override { return false; };
  bool get(Mesh* UNUSED(m), std::string& UNUSED(sval),
           const std::string& UNUSED(name)) override {
    return true;
  };
  bool get(Mesh* UNUSED(m), int& UNUSED(ival), const std::string& UNUSED(name)) override {
    return true;
  };
  bool get(Mesh* UNUSED(m), BoutReal& UNUSED(rval),
           const std::string& UNUSED(name)) override {
    return true;
  }
  bool get(Mesh* UNUSED(m), Field2D& UNUSED(var), const std::string& UNUSED(name),
           BoutReal UNUSED(def) = 0.0) override {
    return true;
  }
  bool get(Mesh* UNUSED(m), Field3D& UNUSED(var), const std::string& UNUSED(name),
           BoutReal UNUSED(def) = 0.0) override {
    return true;
  }
  bool get(Mesh* UNUSED(m), FieldPerp& UNUSED(var), const std::string& UNUSED(name),
           BoutReal UNUSED(def) = 0.0) override {
    return true;
  }
  bool get(Mesh* UNUSED(m), std::vector<int>& UNUSED(var),
           const std::string& UNUSED(name), int UNUSED(len), int UNUSED(offset) = 0,
           Direction UNUSED(dir) = GridDataSource::X) override {
    return true;
  }
  bool get(Mesh* UNUSED(m), std::vector<BoutReal>& UNUSED(var),
           const std::string& UNUSED(name), int UNUSED(len), int UNUSED(offset) = 0,
           Direction UNUSED(dir) = GridDataSource::X) override {
    return true;
  }
  bool hasXBoundaryGuards(Mesh* UNUSED(m)) override { return false; }
  bool hasYBoundaryGuards() override { return false; }
};

TEST(BoutMeshTest, NullOptionsCheck) {
  // Temporarily turn off outputs to make test quiet
  output_info.disable();
  output_warn.disable();
  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
  output_info.enable();
  output_warn.enable();
}
