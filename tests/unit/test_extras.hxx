#ifndef TEST_EXTRAS_H__
#define TEST_EXTRAS_H__

#include "gtest/gtest.h"

#include <numeric>
#include <functional>
#include <iostream>
#include <vector>

#include "boutcomm.hxx"
#include "field3d.hxx"
#include "unused.hxx"
#include "bout/coordinates.hxx"
#include "bout/mesh.hxx"
#include "bout/mpi_wrapper.hxx"
#include "bout/operatorstencil.hxx"

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
    mpi = bout::globals::mpi;
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

  comm_handle send(FieldGroup& UNUSED(g)) override { return nullptr; }
  comm_handle sendX(FieldGroup& UNUSED(g), comm_handle UNUSED(handle) = nullptr,
                    bool UNUSED(disable_corners) = false) override {
    return nullptr;
  }
  comm_handle sendY(FieldGroup& UNUSED(g), comm_handle UNUSED(handle) = nullptr) override
  {
    return nullptr;
  }
  int wait(comm_handle UNUSED(handle)) override { return 0; }
  MPI_Request sendToProc(int UNUSED(xproc), int UNUSED(yproc), BoutReal* UNUSED(buffer),
                         int UNUSED(size), int UNUSED(tag)) override {
    return MPI_Request();
  }
  comm_handle receiveFromProc(int UNUSED(xproc), int UNUSED(yproc),
                              BoutReal* UNUSED(buffer), int UNUSED(size),
                              int UNUSED(tag)) override {
    return nullptr;
  }
  int getNXPE() override { return 1; }
  int getNYPE() override { return 1; }
  int getXProcIndex() override { return 1; }
  int getYProcIndex() override { return 1; }
  bool firstX() const override { return true; }
  bool lastX() const override { return true; }
  int sendXOut(BoutReal* UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) override {
    return 0;
  }
  int sendXIn(BoutReal* UNUSED(buffer), int UNUSED(size), int UNUSED(tag)) override {
    return 0;
  }
  comm_handle irecvXOut(BoutReal* UNUSED(buffer), int UNUSED(size),
                        int UNUSED(tag)) override {
    return nullptr;
  }
  comm_handle irecvXIn(BoutReal* UNUSED(buffer), int UNUSED(size),
                       int UNUSED(tag)) override {
    return nullptr;
  }
  MPI_Comm getXcomm(int UNUSED(jy)) const override { return BoutComm::get(); }
  MPI_Comm getYcomm(int UNUSED(jx)) const override { return BoutComm::get(); }
  bool periodicY(int UNUSED(jx)) const override { return true; }
  bool periodicY(int UNUSED(jx), BoutReal& UNUSED(ts)) const override { return true; }
  int numberOfYBoundaries() const override { return 1; }
  std::pair<bool, BoutReal> hasBranchCutLower(int UNUSED(jx)) const override {
    return std::make_pair(false, 0.);
  }
  std::pair<bool, BoutReal> hasBranchCutUpper(int UNUSED(jx)) const override {
    return std::make_pair(false, 0.);
  }
  bool firstY() const override { return true; }
  bool lastY() const override { return true; }
  bool firstY(int UNUSED(xpos)) const override { return true; }
  bool lastY(int UNUSED(xpos)) const override { return true; }
  int UpXSplitIndex() override { return 0; }
  int DownXSplitIndex() override { return 0; }
  int sendYOutIndest(BoutReal* UNUSED(buffer), int UNUSED(size),
                     int UNUSED(tag)) override {
    return 0;
  }
  int sendYOutOutdest(BoutReal* UNUSED(buffer), int UNUSED(size),
                      int UNUSED(tag)) override {
    return 0;
  }
  int sendYInIndest(BoutReal* UNUSED(buffer), int UNUSED(size),
                    int UNUSED(tag)) override {
    return 0;
  }
  int sendYInOutdest(BoutReal* UNUSED(buffer), int UNUSED(size),
                     int UNUSED(tag)) override {
    return 0;
  }
  comm_handle irecvYOutIndest(BoutReal* UNUSED(buffer), int UNUSED(size),
                              int UNUSED(tag)) override {
    return nullptr;
  }
  comm_handle irecvYOutOutdest(BoutReal* UNUSED(buffer), int UNUSED(size),
                               int UNUSED(tag)) override {
    return nullptr;
  }
  comm_handle irecvYInIndest(BoutReal* UNUSED(buffer), int UNUSED(size),
                             int UNUSED(tag)) override {
    return nullptr;
  }
  comm_handle irecvYInOutdest(BoutReal* UNUSED(buffer), int UNUSED(size),
                              int UNUSED(tag)) override {
    return nullptr;
  }
  RangeIterator iterateBndryLowerY() const override {
    return RangeIterator(xstart, xend);
  }
  RangeIterator iterateBndryUpperY() const override {
    return RangeIterator(xstart, xend);
  }
  RangeIterator iterateBndryLowerOuterY() const override { return RangeIterator(); }
  RangeIterator iterateBndryLowerInnerY() const override { return RangeIterator(); }
  RangeIterator iterateBndryUpperOuterY() const override { return RangeIterator(); }
  RangeIterator iterateBndryUpperInnerY() const override { return RangeIterator(); }
  void addBoundary(BoundaryRegion* region) override { boundaries.push_back(region); }
  std::vector<BoundaryRegion*> getBoundaries() override { return boundaries; }
  std::vector<BoundaryRegionPar*> getBoundariesPar() override {
    return std::vector<BoundaryRegionPar*>();
  }
  BoutReal GlobalX(int jx) const override { return jx; }
  BoutReal GlobalY(int jy) const override { return jy; }
  BoutReal GlobalX(BoutReal jx) const override { return jx; }
  BoutReal GlobalY(BoutReal jy) const override { return jy; }
  int getGlobalXIndex(int) const override { return 0; }
  int getGlobalXIndexNoBoundaries(int) const override { return 0; }
  int getGlobalYIndex(int y) const override { return y; }
  int getGlobalYIndexNoBoundaries(int y) const override { return y; }
  int getGlobalZIndex(int) const override { return 0; }
  int getGlobalZIndexNoBoundaries(int) const override { return 0; }
  int getLocalXIndex(int) const override { return 0; }
  int getLocalXIndexNoBoundaries(int) const override { return 0; }
  int getLocalYIndex(int y) const override { return y; }
  int getLocalYIndexNoBoundaries(int y) const override { return y; }
  int getLocalZIndex(int) const override { return 0; }
  int getLocalZIndexNoBoundaries(int) const override { return 0; }

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
                Region<Ind2D>(0, xstart - 1, ystart, yend, 0, 0, LocalNy, 1));
    addRegion3D("RGN_INNER_X", Region<Ind3D>(0, xstart - 1, ystart, yend, 0, LocalNz - 1,
                                             LocalNy, LocalNz));
    addRegionPerp("RGN_INNER_X",
                  Region<IndPerp>(0, xstart - 1, 0, 0, 0, LocalNz - 1, 1, LocalNz));
    addRegion2D("RGN_OUTER_X",
                Region<Ind2D>(xend + 1, LocalNx - 1, ystart, yend, 0, 0, LocalNy, 1));
    addRegion3D("RGN_OUTER_X", Region<Ind3D>(xend + 1, LocalNx - 1, ystart, yend, 0,
                                             LocalNz - 1, LocalNy, LocalNz));
    addRegionPerp("RGN_OUTER_X", Region<IndPerp>(xend + 1, LocalNx - 1, 0, 0, 0,
                                                 LocalNz - 1, 1, LocalNz));

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
    addRegionPerp("RGN_BNDRY",
                  getRegionPerp("RGN_INNER_X") + getRegionPerp("RGN_OUTER_X"));
  }

  // Make this public so we can test it
  using Mesh::msg_len;

private:
  std::vector<BoundaryRegion *> boundaries;
};

/// FakeGridDataSource provides a non-null GridDataSource* source to use with FakeMesh, to
/// allow testing of methods that use 'source' - in particular allowing
/// source->hasXBoundaryGuards and source->hasXBoundaryGuards to be called.
/// By passing in a set of values, the mesh->get routines can be used.
class FakeGridDataSource : public GridDataSource {
public:
  FakeGridDataSource() {}
  /// Constructor setting values which can be fetched from this source
  FakeGridDataSource(Options& values) : values(values) {}

  /// Take an rvalue (e.g. initializer list), convert to lvalue and delegate constructor
  FakeGridDataSource(Options&& values) : FakeGridDataSource(values) {}

  bool hasVar(const std::string& UNUSED(name)) override { return false; }

  bool get(MAYBE_UNUSED(Mesh* m), std::string& sval, const std::string& name,
           const std::string& def = "") override {
    if (values[name].isSet()) {
      sval = values[name].as<std::string>();
      return true;
    }
    sval = def;
    return false;
  }
  bool get(MAYBE_UNUSED(Mesh* m), int& ival, const std::string& name,
           int def = 0) override {
    if (values[name].isSet()) {
      ival = values[name].as<int>();
      return true;
    }
    ival = def;
    return false;
  }
  bool get(MAYBE_UNUSED(Mesh* m), BoutReal& rval, const std::string& name,
           BoutReal def = 0.0) override {
    if (values[name].isSet()) {
      rval = values[name].as<BoutReal>();
      return true;
    }
    rval = def;
    return false;
  }
  bool get(Mesh* mesh, Field2D& fval, const std::string& name, BoutReal def = 0.0,
           CELL_LOC location = CELL_DEFAULT) override {
    if (values[name].isSet()) {
      fval = values[name].as(Field2D(0.0, mesh));
      return true;
    }
    fval = def;
    return false;
  }
  bool get(Mesh* mesh, Field3D& fval, const std::string& name, BoutReal def = 0.0,
           CELL_LOC location = CELL_DEFAULT) override {
    if (values[name].isSet()) {
      fval = values[name].as(Field3D(0.0, mesh));
      return true;
    }
    fval = def;
    return false;
  }
  bool get(Mesh* mesh, FieldPerp& fval, const std::string& name, BoutReal def = 0.0,
           CELL_LOC location = CELL_DEFAULT) override {
    if (values[name].isSet()) {
      fval = values[name].as(FieldPerp(0.0, mesh));
      return true;
    }
    fval = def;
    return false;
  }

  bool get(MAYBE_UNUSED(Mesh* m), MAYBE_UNUSED(std::vector<int>& var),
           MAYBE_UNUSED(const std::string& name), MAYBE_UNUSED(int len),
           MAYBE_UNUSED(int def) = 0, Direction = GridDataSource::X) override {
    throw BoutException("Not Implemented");
    return false;
  }
  bool get(MAYBE_UNUSED(Mesh* m), MAYBE_UNUSED(std::vector<BoutReal>& var),
           MAYBE_UNUSED(const std::string& name), MAYBE_UNUSED(int len),
           MAYBE_UNUSED(int def) = 0,
           Direction UNUSED(dir) = GridDataSource::X) override {
    throw BoutException("Not Implemented");
    return false;
  }

  bool hasXBoundaryGuards(Mesh* UNUSED(m)) override { return true; }

  bool hasYBoundaryGuards() override { return true; }
private:
  Options values; ///< Store values to be returned by get()
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
    WithQuietOutput quiet_info{output_info};
    WithQuietOutput quiet_warn{output_warn};

    delete bout::globals::mesh;
    bout::globals::mpi = new MpiWrapper();
    bout::globals::mesh = new FakeMesh(nx, ny, nz);
    bout::globals::mesh->createDefaultRegions();
    static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(nullptr);
    test_coords = std::make_shared<Coordinates>(
        bout::globals::mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0});
    // No call to Coordinates::geometry() needed here
    static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(test_coords);
    static_cast<FakeMesh*>(bout::globals::mesh)->setGridDataSource(
        new FakeGridDataSource());
    // May need a ParallelTransform to create fields, because create3D calls
    // fromFieldAligned
    test_coords->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));
    static_cast<FakeMesh*>(bout::globals::mesh)->createBoundaryRegions();

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
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered});
    // No call to Coordinates::geometry() needed here
    test_coords_staggered->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*mesh_staggered));
  }

  ~FakeMeshFixture() override {
    delete bout::globals::mesh;
    bout::globals::mesh = nullptr;
    delete mesh_staggered;
    mesh_staggered = nullptr;
    delete bout::globals::mpi;
    bout::globals::mpi = nullptr;

    Options::cleanup();
  }

  static constexpr int nx = 3;
  static constexpr int ny = 5;
  static constexpr int nz = 7;

  Mesh* mesh_staggered = nullptr;

  std::shared_ptr<Coordinates> test_coords{nullptr};
  std::shared_ptr<Coordinates> test_coords_staggered{nullptr};
};

/// Returns a stencil object which indicates that non-boundary cells
/// depend on all of their neighbours to a depth of one, including
/// corners.
template <class T>
OperatorStencil<T> squareStencil(Mesh* localmesh) {
  OperatorStencil<T> stencil;
  IndexOffset<T> zero;
  std::set<IndexOffset<T>> offsets = {
      zero,
      zero.xp(),
      zero.xm(),
  };
  if (!std::is_same<T, IndPerp>::value) {
    offsets.insert(zero.yp());
    offsets.insert(zero.ym());
    offsets.insert(zero.xp().yp());
    offsets.insert(zero.xp().ym());
    offsets.insert(zero.xm().yp());
    offsets.insert(zero.xm().ym());
  }
  if (!std::is_same<T, Ind2D>::value) {
    offsets.insert(zero.zp());
    offsets.insert(zero.zm());
    offsets.insert(zero.xp().zp());
    offsets.insert(zero.xp().zm());
    offsets.insert(zero.xm().zp());
    offsets.insert(zero.xm().zm());
  }
  if (std::is_same<T, Ind3D>::value) {
    offsets.insert(zero.yp().zp());
    offsets.insert(zero.yp().zm());
    offsets.insert(zero.ym().zp());
    offsets.insert(zero.ym().zm());
  }
  std::vector<IndexOffset<T>> offsetsVec(offsets.begin(), offsets.end());
  stencil.add(
      [localmesh](T ind) -> bool {
        return (localmesh->xstart <= ind.x() && ind.x() <= localmesh->xend
                && (std::is_same<T, IndPerp>::value
                    || (localmesh->ystart <= ind.y() && ind.y() <= localmesh->yend))
                && (std::is_same<T, Ind2D>::value
                    || (localmesh->zstart <= ind.z() && ind.z() <= localmesh->zend)));
      },
      offsetsVec);
  stencil.add([](T UNUSED(ind)) -> bool { return true; }, {zero});
  return stencil;
}

/// Returns a stencil object which indicates that non-boundary cells
/// depend on all of their neighbours to a depth of one, excluding
/// corners.
template <class T>
OperatorStencil<T> starStencil(Mesh* localmesh) {
  OperatorStencil<T> stencil;
  IndexOffset<T> zero;
  std::set<IndexOffset<T>> offsets = {
      zero,
      zero.xp(),
      zero.xm(),
  };
  if (!std::is_same<T, IndPerp>::value) {
    offsets.insert(zero.yp());
    offsets.insert(zero.ym());
  }
  if (!std::is_same<T, Ind2D>::value) {
    offsets.insert(zero.zp());
    offsets.insert(zero.zm());
  }
  std::vector<IndexOffset<T>> offsetsVec(offsets.begin(), offsets.end());
  stencil.add(
      [localmesh](T ind) -> bool {
        return (localmesh->xstart <= ind.x() && ind.x() <= localmesh->xend
                && (std::is_same<T, IndPerp>::value
                    || (localmesh->ystart <= ind.y() && ind.y() <= localmesh->yend))
                && (std::is_same<T, Ind2D>::value
                    || (localmesh->zstart <= ind.z() && ind.z() <= localmesh->zend)));
      },
      offsetsVec);
  stencil.add([](T UNUSED(ind)) -> bool { return true; }, {zero});
  return stencil;
}

#endif //  TEST_EXTRAS_H__
