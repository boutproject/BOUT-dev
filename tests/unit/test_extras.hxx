#ifndef TEST_EXTRAS_H__
#define TEST_EXTRAS_H__

#include "gtest/gtest.h"

#include <numeric>
#include <functional>
#include <iostream>
#include <vector>

#include "boutcomm.hxx"
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

using bout::utils::EnableIfField;

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
template <class T, class U, typename = EnableIfField<T, U>>
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

/// Disable a ConditionalOutput during a scope; reenable it on
/// exit. You must give the variable a name!
///
///     {
///       WithQuietoutput quiet{output};
///       // output disabled during this scope
///     }
///     // output now enabled
class WithQuietOutput {
public:
  explicit WithQuietOutput(ConditionalOutput& output_in) : output(output_in) {
    output.disable();
  }

  ~WithQuietOutput() { output.enable(); }
  ConditionalOutput& output;
};

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
    OffsetX = 0;
    OffsetY = 0;
    OffsetZ = 0;

    // Small "inner" region
    xstart = 1;
    xend = nx - 2;
    ystart = 1;
    yend = ny - 2;
    zstart = 0;
    zend = nz - 1;

    StaggerGrids=false;
    
    // Unused variables
    periodicX = false;
    NXPE = 1;
    PE_XIND = 0;
    IncIntShear = false;
    maxregionblocksize = MAXREGIONBLOCKSIZE;

    // Need some options for parallelTransform
    options = Options::getRoot();
  }

  void setCoordinates(std::shared_ptr<Coordinates> coords, CELL_LOC location = CELL_CENTRE) {
    coords_map[location] = coords;
  }

  void setGridDataSource(GridDataSource* source_in) {
    source = source_in;
  }

  // Use this if the FakeMesh needs x- and y-boundaries
  void createBoundaries() {
    addBoundary(new BoundaryRegionXIn("core", ystart, yend, this));
    addBoundary(new BoundaryRegionXOut("sol", ystart, yend, this));
    addBoundary(new BoundaryRegionYUp("upper_target", xstart, xend, this));
    addBoundary(new BoundaryRegionYDown("lower_target", xstart, xend, this));
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
  std::pair<bool, BoutReal> hasBranchCutLower(int UNUSED(jx)) const {
    return std::make_pair(false, 0.);
  }
  std::pair<bool, BoutReal> hasBranchCutUpper(int UNUSED(jx)) const {
    return std::make_pair(false, 0.);
  }
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
  int XLOCAL(int UNUSED(xglo)) const { return 0; }
  int YLOCAL(int UNUSED(yglo)) const { return 0; }

  void initDerivs(Options * opt){
    StaggerGrids=true;
    derivs_init(opt);
  }

  void createBoundaryRegions() {
    addRegion2D("RGN_LOWER_Y",
                Region<Ind2D>(0, LocalNx - 1, 0, ystart - 1, 0, 0, LocalNy, 1));
    addRegion3D("RGN_LOWER_Y", Region<Ind3D>(0, LocalNx - 1, 0, ystart - 1, 0,
                                             LocalNz - 1, LocalNy, LocalNz));
    addRegion2D("RGN_UPPER_Y",
                Region<Ind2D>(0, LocalNx - 1, yend + 1, LocalNy - 1, 0, 0, LocalNy, 1));
    addRegion3D("RGN_UPPER_Y", Region<Ind3D>(0, LocalNx - 1, yend + 1, LocalNy - 1, 0,
                                             LocalNz - 1, LocalNy, LocalNz));
    addRegion2D("RGN_INNER_X",
                Region<Ind2D>(0, xstart - 1, 0, LocalNy - 1, 0, 0, LocalNy, 1));
    addRegion3D("RGN_INNER_X", Region<Ind3D>(0, xstart - 1, 0, LocalNy - 1, 0,
                                             LocalNz - 1, LocalNy, LocalNz));
    addRegion2D("RGN_OUTER_X",
                Region<Ind2D>(xend + 1, LocalNx - 1, 0, LocalNy - 1, 0, 0, LocalNy, 1));
    addRegion3D("RGN_OUTER_X", Region<Ind3D>(xend + 1, LocalNx - 1, 0, LocalNy - 1, 0,
                                             LocalNz - 1, LocalNy, LocalNz));

    const auto boundary_names = {"RGN_LOWER_Y", "RGN_UPPER_Y", "RGN_INNER_X",
                                 "RGN_OUTER_X"};

    // Sum up and get unique points in the boundaries defined above
    addRegion2D("RGN_BNDRY",
                std::accumulate(begin(boundary_names), end(boundary_names),
                                Region<Ind2D>{},
                                [this](Region<Ind2D>& a, const std::string& b) {
                                  return a + getRegion2D(b);
                                })
                    .unique());

    addRegion3D("RGN_BNDRY",
                std::accumulate(begin(boundary_names), end(boundary_names),
                                Region<Ind3D>{},
                                [this](Region<Ind3D>& a, const std::string& b) {
                                  return a + getRegion3D(b);
                                })
                    .unique());
  }

private:
  std::vector<BoundaryRegion *> boundaries;
};

/// FakeGridDataSource provides a non-null GridDataSource* source to use with FakeMesh, to
/// allow testing of methods that use 'source' - in particular allowing
/// source->hasXBoundaryGuards and source->hasXBoundaryGuards to be called.
class FakeGridDataSource : public GridDataSource {
  bool hasVar(const std::string& UNUSED(name)) { return false; }

  bool get(Mesh* UNUSED(m), std::string& UNUSED(sval), const std::string& UNUSED(name)) {
    return false;
  }
  bool get(Mesh* UNUSED(m), int& UNUSED(ival), const std::string& UNUSED(name)) {
    return false;
  }
  bool get(Mesh* UNUSED(m), BoutReal& UNUSED(rval), const std::string& UNUSED(name)) {
    return false;
  }
  bool get(Mesh* UNUSED(m), Field2D& UNUSED(var), const std::string& UNUSED(name),
      BoutReal UNUSED(def) = 0.0) {
    return false;
  }
  bool get(Mesh* UNUSED(m), Field3D& UNUSED(var), const std::string& UNUSED(name),
      BoutReal UNUSED(def) = 0.0) {
    return false;
  }
  bool get(Mesh* UNUSED(m), FieldPerp& UNUSED(var), const std::string& UNUSED(name),
      BoutReal UNUSED(def) = 0.0) {
    return false;
  }

  bool get(Mesh* UNUSED(m), std::vector<int>& UNUSED(var),
      const std::string& UNUSED(name), int UNUSED(len), int UNUSED(offset) = 0,
      Direction UNUSED(dir) = GridDataSource::X) {
    return false;
  }
  bool get(Mesh* UNUSED(m), std::vector<BoutReal>& UNUSED(var),
      const std::string& UNUSED(name), int UNUSED(len), int UNUSED(offset) = 0,
      Direction UNUSED(dir) = GridDataSource::X) {
    return false;
  }

  bool hasXBoundaryGuards(Mesh* UNUSED(m)) { return true; }

  bool hasYBoundaryGuards() { return true; }
};

/// Test fixture to make sure the global mesh is our fake
/// one. Also initialize the global mesh_staggered for use in tests with
/// staggering. Multiple tests have exactly the same fixture, so use a type
/// alias to make a new test:
///
///     using MyTest = FakeMeshFixture;
class FakeMeshFixture : public ::testing::Test {
public:
  FakeMeshFixture() {
    WithQuietOutput quiet{output_info};

    delete bout::globals::mesh;
    bout::globals::mesh = new FakeMesh(nx, ny, nz);
    bout::globals::mesh->createDefaultRegions();
    static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(nullptr);
    test_coords = std::make_shared<Coordinates>(
        bout::globals::mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
        false);
    static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(test_coords);
    static_cast<FakeMesh*>(bout::globals::mesh)->setGridDataSource(
        new FakeGridDataSource());
    // May need a ParallelTransform to create fields, because create3D calls
    // fromFieldAligned
    test_coords->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));

    delete mesh_staggered;
    mesh_staggered = new FakeMesh(nx, ny, nz);
    mesh_staggered->StaggerGrids = true;
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr);
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_XLOW);
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_YLOW);
    static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_ZLOW);
    mesh_staggered->createDefaultRegions();

    test_coords_staggered = std::make_shared<Coordinates>(
        mesh_staggered, Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        BoutReal{1.0}, Field2D{1.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered}, false);
    test_coords_staggered->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*mesh_staggered));
  }

  virtual ~FakeMeshFixture() {
    delete bout::globals::mesh;
    bout::globals::mesh = nullptr;
    delete mesh_staggered;
    mesh_staggered = nullptr;
  }

  static constexpr int nx = 3;
  static constexpr int ny = 5;
  static constexpr int nz = 7;

  Mesh* mesh_staggered = nullptr;

  std::shared_ptr<Coordinates> test_coords{nullptr};
  std::shared_ptr<Coordinates> test_coords_staggered{nullptr};
};

#endif //  TEST_EXTRAS_H__
