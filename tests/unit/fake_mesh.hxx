#pragma once

#include <bout/boundary_region.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/coordinates.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/fieldgroup.hxx>
#include <bout/globals.hxx>
#include <bout/griddata.hxx>
#include <bout/mesh.hxx>
#include <bout/mpi_wrapper.hxx>
#include <bout/region.hxx>
#include <bout/sys/range.hxx>
#include <bout/unused.hxx>

#include <iterator>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

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
  FakeMesh(int nx, int ny, int nz, MpiWrapper& mpi_in) {
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

    StaggerGrids = false;

    // Unused variables
    periodicX = false;
    NXPE = 1;
    PE_XIND = 0;
    IncIntShear = false;
    maxregionblocksize = MAXREGIONBLOCKSIZE;

    // Need some options for parallelTransform
    options = Options::getRoot();
    mpi = &mpi_in;
  }
  FakeMesh(int nx, int ny, int nz) : FakeMesh(nx, ny, nz, *bout::globals::mpi) {}

  void setCoordinates(std::shared_ptr<Coordinates> coords,
                      CELL_LOC location = CELL_CENTRE) {
    coords_map[location] = coords;
  }

  void setGridDataSource(GridDataSource* source_in) { source = source_in; }

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
  comm_handle sendY(FieldGroup& UNUSED(g),
                    comm_handle UNUSED(handle) = nullptr) override {
    return nullptr;
  }
  int wait(comm_handle UNUSED(handle)) override { return 0; }
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
  std::vector<std::shared_ptr<BoundaryRegionPar>>
  getBoundariesPar(BoundaryParType UNUSED(type)) override {
    return std::vector<std::shared_ptr<BoundaryRegionPar>>();
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

  void initDerivs(Options* opt) {
    StaggerGrids = true;
    derivs_init(opt);
  }

  void createBoundaryRegions() {
    addRegion2D("RGN_LOWER_Y",
                Region<Ind2D>(0, LocalNx - 1, 0, ystart - 1, 0, 0, LocalNy, 1));
    addRegion3D("RGN_LOWER_Y", Region<Ind3D>(0, LocalNx - 1, 0, ystart - 1, 0,
                                             LocalNz - 1, LocalNy, LocalNz));
    addRegion2D("RGN_LOWER_Y_THIN", getRegion2D("RGN_LOWER_Y"));
    addRegion3D("RGN_LOWER_Y_THIN", getRegion3D("RGN_LOWER_Y"));

    addRegion2D("RGN_UPPER_Y",
                Region<Ind2D>(0, LocalNx - 1, yend + 1, LocalNy - 1, 0, 0, LocalNy, 1));
    addRegion3D("RGN_UPPER_Y", Region<Ind3D>(0, LocalNx - 1, yend + 1, LocalNy - 1, 0,
                                             LocalNz - 1, LocalNy, LocalNz));
    addRegion2D("RGN_UPPER_Y_THIN", getRegion2D("RGN_UPPER_Y"));
    addRegion3D("RGN_UPPER_Y_THIN", getRegion3D("RGN_UPPER_Y"));

    addRegion2D("RGN_INNER_X",
                Region<Ind2D>(0, xstart - 1, ystart, yend, 0, 0, LocalNy, 1));
    addRegion3D("RGN_INNER_X", Region<Ind3D>(0, xstart - 1, ystart, yend, 0, LocalNz - 1,
                                             LocalNy, LocalNz));
    addRegionPerp("RGN_INNER_X",
                  Region<IndPerp>(0, xstart - 1, 0, 0, 0, LocalNz - 1, 1, LocalNz));
    addRegionPerp("RGN_INNER_X_THIN",
                  Region<IndPerp>(0, xstart - 1, 0, 0, 0, LocalNz - 1, 1, LocalNz));
    addRegion2D("RGN_INNER_X_THIN", getRegion2D("RGN_INNER_X"));
    addRegion3D("RGN_INNER_X_THIN", getRegion3D("RGN_INNER_X"));

    addRegion2D("RGN_OUTER_X",
                Region<Ind2D>(xend + 1, LocalNx - 1, ystart, yend, 0, 0, LocalNy, 1));
    addRegion3D("RGN_OUTER_X", Region<Ind3D>(xend + 1, LocalNx - 1, ystart, yend, 0,
                                             LocalNz - 1, LocalNy, LocalNz));
    addRegionPerp("RGN_OUTER_X", Region<IndPerp>(xend + 1, LocalNx - 1, 0, 0, 0,
                                                 LocalNz - 1, 1, LocalNz));
    addRegionPerp("RGN_OUTER_X_THIN", Region<IndPerp>(xend + 1, LocalNx - 1, 0, 0, 0,
                                                      LocalNz - 1, 1, LocalNz));
    addRegion2D("RGN_OUTER_X_THIN", getRegion2D("RGN_OUTER_X"));
    addRegion3D("RGN_OUTER_X_THIN", getRegion3D("RGN_OUTER_X"));

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
  std::vector<BoundaryRegion*> boundaries;
};

/// FakeGridDataSource provides a non-null GridDataSource* source to use with FakeMesh, to
/// allow testing of methods that use 'source' - in particular allowing
/// source->hasXBoundaryGuards and source->hasXBoundaryGuards to be called.
/// By passing in a set of values, the mesh->get routines can be used.
class FakeGridDataSource : public GridDataSource {
public:
  FakeGridDataSource() {}
  /// Constructor setting values which can be fetched from this source
  FakeGridDataSource(const Options& values) : values(values.copy()) {}

  /// Take an rvalue (e.g. initializer list), convert to lvalue and delegate constructor
  FakeGridDataSource(Options&& values) : FakeGridDataSource(values) {}

  bool hasVar(const std::string& UNUSED(name)) override { return false; }

  bool get([[maybe_unused]] Mesh* m, std::string& sval, const std::string& name,
           const std::string& def = "") override {
    if (values[name].isSet()) {
      sval = values[name].as<std::string>();
      return true;
    }
    sval = def;
    return false;
  }
  bool get([[maybe_unused]] Mesh* m, int& ival, const std::string& name,
           int def = 0) override {
    if (values[name].isSet()) {
      ival = values[name].as<int>();
      return true;
    }
    ival = def;
    return false;
  }
  bool get([[maybe_unused]] Mesh* m, BoutReal& rval, const std::string& name,
           BoutReal def = 0.0) override {
    if (values[name].isSet()) {
      rval = values[name].as<BoutReal>();
      return true;
    }
    rval = def;
    return false;
  }
  bool get(Mesh* mesh, Field2D& fval, const std::string& name, BoutReal def = 0.0,
           CELL_LOC UNUSED(location) = CELL_DEFAULT) override {
    if (values[name].isSet()) {
      fval = values[name].as(Field2D(0.0, mesh));
      return true;
    }
    fval = def;
    return false;
  }
  bool get(Mesh* mesh, Field3D& fval, const std::string& name, BoutReal def = 0.0,
           CELL_LOC UNUSED(location) = CELL_DEFAULT) override {
    if (values[name].isSet()) {
      fval = values[name].as(Field3D(0.0, mesh));
      return true;
    }
    fval = def;
    return false;
  }
  bool get(Mesh* mesh, FieldPerp& fval, const std::string& name, BoutReal def = 0.0,
           CELL_LOC UNUSED(location) = CELL_DEFAULT) override {
    if (values[name].isSet()) {
      fval = values[name].as(FieldPerp(0.0, mesh));
      return true;
    }
    fval = def;
    return false;
  }

  bool get([[maybe_unused]] Mesh* m, [[maybe_unused]] std::vector<int>& var,
           [[maybe_unused]] const std::string& name, [[maybe_unused]] int len,
           [[maybe_unused]] int def = 0, Direction = GridDataSource::X) override {
    throw BoutException("Not Implemented");
    return false;
  }
  bool get([[maybe_unused]] Mesh* m, [[maybe_unused]] std::vector<BoutReal>& var,
           [[maybe_unused]] const std::string& name, [[maybe_unused]] int len,
           [[maybe_unused]] int def = 0,
           Direction UNUSED(dir) = GridDataSource::X) override {
    throw BoutException("Not Implemented");
    return false;
  }

  bool hasXBoundaryGuards(Mesh* UNUSED(m)) override { return true; }

  bool hasYBoundaryGuards() override { return true; }

private:
  Options values; ///< Store values to be returned by get()
};
