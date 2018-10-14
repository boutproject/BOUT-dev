#ifndef TEST_EXTRAS_H__
#define TEST_EXTRAS_H__

#include "gtest/gtest.h"

#include <iostream>
#include <mpi.h>

#include "bout/mesh.hxx"
#include "field3d.hxx"
#include "unused.hxx"

const BoutReal BoutRealTolerance = 1e-15;

/// Does \p str contain \p substring?
::testing::AssertionResult IsSubString(const std::string &str,
                                       const std::string &substring);

/// Is \p field equal to \p number, with a tolerance of \p tolerance?
::testing::AssertionResult IsField3DEqualBoutReal(const Field3D &field, BoutReal number,
                                                  BoutReal tolerance = BoutRealTolerance);

/// Is \p field equal to \p number, with a tolerance of \p tolerance?
::testing::AssertionResult IsField2DEqualBoutReal(const Field2D &field, BoutReal number,
                                                  BoutReal tolerance = BoutRealTolerance);

/// Is \p field equal to \p number, with a tolerance of \p tolerance?
::testing::AssertionResult IsFieldPerpEqualBoutReal(const FieldPerp &field, BoutReal number,
                                                  BoutReal tolerance = BoutRealTolerance);


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
  bool periodicY(int UNUSED(jx)) { return true; }
  bool periodicY(int UNUSED(jx), BoutReal &UNUSED(ts), CELL_LOC UNUSED(location)) { return true; }
  bool hasBranchCut() const { return false; }
  bool hasBranchCutDown(int UNUSED(jx)) const { return false; }
  bool hasBranchCutUp(int UNUSED(jx)) const { return false; }
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
  vector<BoundaryRegion *> getBoundaries() { return boundaries; }
  vector<BoundaryRegionPar *> getBoundariesPar() { return vector<BoundaryRegionPar *>(); }
  BoutReal GlobalX(int UNUSED(jx)) const { return 0; }
  BoutReal GlobalY(int UNUSED(jy)) const { return 0; }
  BoutReal GlobalZ(int UNUSED(jz)) const { return 0; }
  BoutReal GlobalX(BoutReal UNUSED(jx)) const { return 0; }
  BoutReal GlobalY(BoutReal UNUSED(jy)) const { return 0; }
  BoutReal GlobalZ(BoutReal UNUSED(jz)) const { return 0; }
  int XGLOBAL(int UNUSED(xloc)) const { return 0; }
  int YGLOBAL(int UNUSED(yloc)) const { return 0; }

  void initDerivs(Options * opt){
    StaggerGrids=true;
    derivs_init(opt);
  }
private:
  vector<BoundaryRegion *> boundaries;
};



#endif //  TEST_EXTRAS_H__
