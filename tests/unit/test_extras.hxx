#ifndef TEST_EXTRAS_H__
#define TEST_EXTRAS_H__

#include "gtest/gtest.h"

#include <functional>
#include <iostream>
#include <mpi.h>
#include <vector>

#include "bout/mesh.hxx"
#include "bout/coordinates.hxx"
#include "field3d.hxx"
#include "unused.hxx"

static constexpr BoutReal BoutRealTolerance{1e-15};
// FFTs have a slightly looser tolerance than other functions
static constexpr BoutReal FFTTolerance{1.e-12};

/// Does \p str contain \p substring?
::testing::AssertionResult IsSubString(const std::string &str,
                                       const std::string &substring);

void fillField(Field3D& f, std::vector<std::vector<std::vector<BoutReal>>> values);
void fillField(Field2D& f, std::vector<std::vector<BoutReal>> values);

/// Enable a function if T is a subclass of Field
template <class T>
using EnableIfField = typename std::enable_if<std::is_base_of<Field, T>::value>::type;

/// Returns a field filled with the result of \p fill_function at each point
/// Arbitrary arguments can be passed to the field constructor
template <class T, class... Args, typename = EnableIfField<T>>
T makeField(const std::function<BoutReal(typename T::ind_type&)>& fill_function,
            Args&&... args) {
  T result{std::forward<Args>(args)...};
  result.allocate();

  for (auto i : result) {
    result[i] = fill_function(i);
  }

  return result;
}

/// Teach googletest how to print SpecificInds
template<IND_TYPE N>
inline std::ostream& operator<< (std::ostream &out, const SpecificInd<N> &index) {
  return out << index.ind;
}

/// Helpers to get the type of a Field as a string
auto inline getFieldType(MAYBE_UNUSED(const Field2D& field)) -> std::string {
  return "Field2D";
}
auto inline getFieldType(MAYBE_UNUSED(const Field3D& field)) -> std::string {
  return "Field3D";
}
auto inline getFieldType(MAYBE_UNUSED(const FieldPerp& field)) -> std::string {
  return "FieldPerp";
}

/// Helpers to get the (x, y, z) index values, along with the
/// single-index of a Field index
auto inline getIndexXYZ(const Ind2D& index) -> std::string {
  std::stringstream ss;
  ss << index.x() << ", " << index.y() << "; [" << index.ind << "]";
  return ss.str();
}
auto inline getIndexXYZ(const Ind3D& index) -> std::string {
  std::stringstream ss;
  ss << index.x() << ", " << index.y() << ", " << index.z() << "; [" << index.ind << "]";
  return ss.str();
}
auto inline getIndexXYZ(const IndPerp& index) -> std::string {
  std::stringstream ss;
  ss << index.x() << ", " << index.y() << ", " << index.z() << "; [" << index.ind << "]";
  return ss.str();
}

/// Is \p field equal to \p reference, with a tolerance of \p tolerance?
template <class T, class U, typename = EnableIfField<T>, typename = EnableIfField<U>>
auto IsFieldEqual(const T& field, const U& reference,
                  const std::string& region = "RGN_ALL",
                  BoutReal tolerance = BoutRealTolerance) -> ::testing::AssertionResult {
  for (auto i : field.getRegion(region)) {
    if (fabs(field[i] - reference[i]) > tolerance) {
      return ::testing::AssertionFailure()
             << getFieldType(field) << "(" << getIndexXYZ(i) << ") == " << field[i]
             << "; Expected: " << reference[i];
    }
  }
  return ::testing::AssertionSuccess();
}

/// Is \p field equal to \p reference, with a tolerance of \p tolerance?
/// Overload for BoutReals
template <class T, typename = EnableIfField<T>>
auto IsFieldEqual(const T& field, BoutReal reference,
                  const std::string& region = "RGN_ALL",
                  BoutReal tolerance = BoutRealTolerance) -> ::testing::AssertionResult {
  for (auto i : field.getRegion(region)) {
    if (fabs(field[i] - reference) > tolerance) {
      return ::testing::AssertionFailure()
             << getFieldType(field) << "(" << getIndexXYZ(i) << ") == " << field[i]
             << "; Expected: " << reference;
    }
  }
  return ::testing::AssertionSuccess();
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
  }

  void setCoordinates(std::shared_ptr<Coordinates> coords, CELL_LOC location = CELL_CENTRE) {
    coords_map[location] = coords;
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
  BoutReal GlobalX(int jx) const { return jx; }
  BoutReal GlobalY(int jy) const { return jy; }
  BoutReal GlobalX(BoutReal jx) const { return jx; }
  BoutReal GlobalY(BoutReal jy) const { return jy; }
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
    if (bout::globals::mesh != nullptr) {
      delete bout::globals::mesh;
      bout::globals::mesh = nullptr;
    }
    bout::globals::mesh = new FakeMesh(nx, ny, nz);
    output_info.disable();
    bout::globals::mesh->createDefaultRegions();
    output_info.enable();
  }

  virtual ~FakeMeshFixture() {
    delete bout::globals::mesh;
    bout::globals::mesh = nullptr;
  }

  static constexpr int nx = 3;
  static constexpr int ny = 5;
  static constexpr int nz = 7;
};

#endif //  TEST_EXTRAS_H__
