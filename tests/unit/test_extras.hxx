#ifndef TEST_EXTRAS_H__
#define TEST_EXTRAS_H__

#include "gtest/gtest.h"

#include <iostream>
#include <mpi.h>
#include <vector>

#include "bout/mesh.hxx"
#include "bout/coordinates.hxx"
#include "field3d.hxx"
#include "unused.hxx"

const BoutReal BoutRealTolerance = 1e-15;

/// Does \p str contain \p substring?
::testing::AssertionResult IsSubString(const std::string &str,
                                       const std::string &substring);

/// Is \p field equal to \p number, with a tolerance of \p tolerance?
::testing::AssertionResult IsField3DEqualBoutReal(const Field3D &field, BoutReal number,
                                                  BoutReal tolerance = BoutRealTolerance);

::testing::AssertionResult IsField3DEqualField3D(const Field3D &lhs, const Field3D &rhs,
                                                 const std::string& region = "RGN_ALL",
                                                 BoutReal tolerance = BoutRealTolerance);

/// Is \p field equal to \p number, with a tolerance of \p tolerance?
::testing::AssertionResult IsField2DEqualBoutReal(const Field2D &field, BoutReal number,
                                                  BoutReal tolerance = BoutRealTolerance);

::testing::AssertionResult IsField2DEqualField2D(const Field2D &lhs, const Field2D &rhs,
                                                 const std::string& region = "RGN_ALL",
                                                 BoutReal tolerance = BoutRealTolerance);

/// Is \p field equal to \p number, with a tolerance of \p tolerance?
::testing::AssertionResult IsFieldPerpEqualBoutReal(const FieldPerp &field, BoutReal number,
                                                  BoutReal tolerance = BoutRealTolerance);

void fillField(Field3D& f, std::vector<std::vector<std::vector<BoutReal>>> values);
void fillField(Field2D& f, std::vector<std::vector<BoutReal>> values);

/// Teach googletest how to print SpecificInds
template<IND_TYPE N>
inline std::ostream& operator<< (std::ostream &out, const SpecificInd<N> &index) {
  return out << index.ind;
}

class Options;

/// FakeMesh has just enough information to create fields
///
/// Notes:
///
/// - This is a mesh for a single process, so the global and local
/// indices are the same.
///
/// - There is a single guard cell at each of the start/end x/y grids.
///
/// - Only the **grid** information is assumed to be used -- anything
///   else will likely **not** work!
class FakeMesh : public Mesh {
public:
  FakeMesh(int nx, int ny, int nz) {
    // Mesh only on one process, so global and local indices are the
    // same
    GlobalNx = nx;
    GlobalNy = ny;
    GlobalNz = nz;
    LocalNx = nx;
    LocalNy = ny;
    LocalNz = nz;
    // Small "inner" region
    xstart = 1;
    xend = nx - 2;
    ystart = 1;
    yend = ny - 2;

    StaggerGrids=true;
    // Unused variables
    periodicX = false;
    NXPE = 1;
    PE_XIND = 0;
    StaggerGrids = false;
    IncIntShear = false;
    maxregionblocksize = MAXREGIONBLOCKSIZE;

    coords_map.emplace(CELL_CENTRE, std::unique_ptr<Coordinates>(nullptr));
    coords_map.emplace(CELL_XLOW, std::unique_ptr<Coordinates>(nullptr));
    coords_map.emplace(CELL_YLOW, std::unique_ptr<Coordinates>(nullptr));
    coords_map.emplace(CELL_ZLOW, std::unique_ptr<Coordinates>(nullptr));
  }

  void setCoordinates(Coordinates* coords, CELL_LOC location = CELL_CENTRE) {
    coords_map[location].reset(coords);
  }

  comm_handle send(FieldGroup &UNUSED(g)) { return nullptr; };
  int wait(comm_handle UNUSED(handle)) { return 0; }
  MPI_Request sendToProc(int UNUSED(xproc), int UNUSED(yproc), BoutReal *UNUSED(buffer),
                         int UNUSED(size), int UNUSED(tag)) {
    return MPI_Request();
  }
  comm_handle receiveFromProc(int UNUSED(xproc), int UNUSED(yproc),
                              BoutReal *UNUSED(buffer), int UNUSED(size),
                              int UNUSED(tag)) {
    return nullptr;
  }
  int getNXPE() { return 1; }
  int getNYPE() { return 1; }
  int getXProcIndex() { return 1; }
  int getYProcIndex() { return 1; }
  bool firstX() { return true; }
  bool lastX() { return true; }
  int sendXOut(BoutReal *UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) { return 0; }
  int sendXIn(BoutReal *UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) { return 0; }
  comm_handle irecvXOut(BoutReal *UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) {
    return nullptr;
  }
  comm_handle irecvXIn(BoutReal *UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) {
    return nullptr;
  }
  MPI_Comm getXcomm(int UNUSED(jy)) const { return MPI_COMM_NULL; }
  MPI_Comm getYcomm(int UNUSED(jx)) const { return MPI_COMM_NULL; }
  bool periodicY(int UNUSED(jx)) const { return true; }
  bool periodicY(int UNUSED(jx), BoutReal &UNUSED(ts)) const { return true; }
  bool firstY() const { return true; }
  bool lastY() const { return true; }
  bool firstY(int UNUSED(xpos)) const { return true; }
  bool lastY(int UNUSED(xpos)) const { return true; }
  int UpXSplitIndex() { return 0; }
  int DownXSplitIndex() { return 0; }
  int sendYOutIndest(BoutReal *UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) {
    return 0;
  }
  int sendYOutOutdest(BoutReal *UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) {
    return 0;
  }
  int sendYInIndest(BoutReal *UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) {
    return 0;
  }
  int sendYInOutdest(BoutReal *UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) {
    return 0;
  }
  comm_handle irecvYOutIndest(BoutReal *UNUSED(buffer), int UNUSED(size),
                              int UNUSED(tag)) {
    return nullptr;
  }
  comm_handle irecvYOutOutdest(BoutReal *UNUSED(buffer), int UNUSED(size),
                               int UNUSED(tag)) {
    return nullptr;
  }
  comm_handle irecvYInIndest(BoutReal *UNUSED(buffer), int UNUSED(size),
                             int UNUSED(tag)) {
    return nullptr;
  }
  comm_handle irecvYInOutdest(BoutReal *UNUSED(buffer), int UNUSED(size),
                              int UNUSED(tag)) {
    return nullptr;
  }
  const RangeIterator iterateBndryLowerY() const { return RangeIterator(); }
  const RangeIterator iterateBndryUpperY() const { return RangeIterator(); }
  const RangeIterator iterateBndryLowerOuterY() const { return RangeIterator(); }
  const RangeIterator iterateBndryLowerInnerY() const { return RangeIterator(); }
  const RangeIterator iterateBndryUpperOuterY() const { return RangeIterator(); }
  const RangeIterator iterateBndryUpperInnerY() const { return RangeIterator(); }
  void addBoundary(BoundaryRegion* region) {boundaries.push_back(region);}
  std::vector<BoundaryRegion *> getBoundaries() { return boundaries; }
  std::vector<BoundaryRegionPar *> getBoundariesPar() { return std::vector<BoundaryRegionPar *>(); }
  BoutReal GlobalX(int UNUSED(jx)) const { return 0; }
  BoutReal GlobalY(int UNUSED(jy)) const { return 0; }
  BoutReal GlobalX(BoutReal UNUSED(jx)) const { return 0; }
  BoutReal GlobalY(BoutReal UNUSED(jy)) const { return 0; }
  int XGLOBAL(int UNUSED(xloc)) const { return 0; }
  int YGLOBAL(int UNUSED(yloc)) const { return 0; }

  void initDerivs(Options * opt){
    StaggerGrids=true;
    derivs_init(opt);
  }
private:
  std::vector<BoundaryRegion *> boundaries;
};

/// Test fixture to make sure the global mesh is our fake
/// one. Multiple tests have exactly the same fixture, so use a type
/// alias to make a new test:
///
///     using MyTest = FakeMeshFixture;
class FakeMeshFixture : public ::testing::Test {
public:
  FakeMeshFixture() {
    // Delete any existing mesh
    if (mesh != nullptr) {
      delete mesh;
      mesh = nullptr;
    }
    mesh = new FakeMesh(nx, ny, nz);
    output_info.disable();
    mesh->createDefaultRegions();
    output_info.enable();
  }

  ~FakeMeshFixture() {
    delete mesh;
    mesh = nullptr;
  }

  static constexpr int nx = 3;
  static constexpr int ny = 5;
  static constexpr int nz = 7;
};

#endif //  TEST_EXTRAS_H__
