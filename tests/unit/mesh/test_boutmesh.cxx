#include "gtest/gtest.h"

#include "bout/mesh.hxx"
#include "../src/mesh/impls/bout/boutmesh.hxx"

#include "test_extras.hxx"

/// Test fixture to make sure the global mesh is our fake one
class BoutMeshTest : public ::testing::Test {
public:
  BoutMeshTest(){};
  static const int nx = 3;
  static const int ny = 5;
  static const int nz = 7;
};

class FakeGridDataSource: public GridDataSource {
public:
  FakeGridDataSource(){};
  virtual ~FakeGridDataSource(){};
  virtual bool hasVar(const string &name) {return false;} ;
  //virtual bool GridDataSource::get(Mesh*, int&, const string&)
  virtual bool get(Mesh *m, int &ival,      const string &name) {
    ival =1;
    return true;};
  virtual bool get(Mesh *m, BoutReal &rval, const string &name) {
    rval = 1;
    return true;}
  virtual bool get(Mesh *m, Field2D &var,   const string &name, BoutReal def=0.0){
    var=def;
    return true;}
  virtual bool get(Mesh *m, Field3D &var,   const string &name, BoutReal def=0.0){
    var=def;
    return true;}
  virtual bool get(Mesh *m, vector<int> &var,      const string &name, int len, int offset=0, Direction dir = GridDataSource::X){
    return true;
  }
  virtual bool get(Mesh *m, vector<BoutReal> &var,      const string &name, int len, int offset=0, Direction dir = GridDataSource::X){
    return true;
  }
};

TEST_F(BoutMeshTest, CreateBoutMesh) {
  EXPECT_NO_THROW(BoutMesh msh(new FakeGridDataSource,nullptr));
  
}
