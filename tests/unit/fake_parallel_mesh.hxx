#ifndef FAKE_PARALLEL_MESH_H__
#define FAKE_PARALLEL_MESH_H__
#include "gtest/gtest.h"
#include <iostream>

#include <functional>
#include <map>
#include <numeric>
#include <vector>

#include "../../src/mesh/impls/bout/boutmesh.hxx"
#include "boutcomm.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "fieldperp.hxx"
#include "unused.hxx"
#include "bout/coordinates.hxx"
#include "bout/fieldgroup.hxx"
#include "bout/mesh.hxx"

class Options;

/// FakeMesh has just enough information to create fields and test how
/// processes communicate with each other. Rather than actually
/// creating multiple processes, it creates multiple mesh objects
/// which represent the subdomain on each processor. These can then be
/// used to create "local" fields. If registered with their corresponding
/// mesh, the local fields will be able to communicate guard regions
/// with each other.
///
/// Notes:
///
/// - All processes must be of uniform size.
///
/// - There is a single guard cell at each of the start/end x/y grids.
///
/// - Only the **grid** and **some of the communication** information
///   is assumed to be used -- anything else will likely **not** work!
///
/// - Only a simple, rectangular topology is available
///
class FakeParallelMesh : public BoutMesh {
public:
  FakeParallelMesh(int nx, int ny, int nz, int nxpe, int nype, int pe_xind, int pe_yind)
      : BoutMesh((nxpe * (nx - 2)) + 2, nype * ny, nz, 1, 1, nxpe, nype, pe_xind,
                 pe_yind),
        yUpMesh(nullptr), yDownMesh(nullptr), xInMesh(nullptr), xOutMesh(nullptr) {
    StaggerGrids = false;
    periodicX = false;
    IncIntShear = false;
    calcParallelSlices_on_communicate = true;
    options = Options::getRoot();
    wait_any_count = -1;
  }

  void initDerivs(Options* opt) {
    StaggerGrids = true;
    derivs_init(opt);
  }

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

  comm_handle send(FieldGroup& g) override {
    overlapHandleMemory(yUpMesh, yDownMesh, xInMesh, xOutMesh);
    std::vector<int> ids;
    int i = 0;
    for (auto f : g) {
      ids.push_back(registeredFields[f]);
      ASSERT1(registeredFieldIds[ids[i]] == f);
      i++;
    }
    if (yUpMesh != nullptr && yUpMesh != this) {
      FieldGroup yUpGroup = makeGroup(yUpMesh, ids);
      yUpMesh->parentSend(yUpGroup);
    }
    if (yDownMesh != nullptr && yDownMesh != this) {
      FieldGroup yDownGroup = makeGroup(yDownMesh, ids);
      yDownMesh->parentSend(yDownGroup);
    }
    if (xInMesh != nullptr && xInMesh != this) {
      FieldGroup xInGroup = makeGroup(xInMesh, ids);
      xInMesh->parentSend(xInGroup);
    }
    if (xOutMesh != nullptr && xOutMesh != this) {
      FieldGroup xOutGroup = makeGroup(xOutMesh, ids);
      xOutMesh->parentSend(xOutGroup);
    }
    return parentSend(g);
  }

  // Need to override this functions to trick mesh into communicating for
  // FieldPerp type
  void communicate(FieldPerp& f) override {
    int nin = xstart;              // Number of x points in inner guard cell
    int nout = LocalNx - xend - 1; // Number of x points in outer guard cell

    if (registeredFieldPerps.count(&f) != 0) {
      int id = registeredFieldPerps[&f];

      if (xInMesh != nullptr && xInMesh->registeredFieldPerpIds.count(id) != 0) {
        FieldPerp* xInField = xInMesh->registeredFieldPerpIds[id];
        for (int i = 0; i < nin * LocalNz; i++) {
          IndPerp ind(i, 1, LocalNz);
          f[ind] = (*xInField)[ind.xp(xend - xstart + 1)];
        }
      }

      if (xOutMesh != nullptr && xOutMesh->registeredFieldPerpIds.count(id) != 0) {
        FieldPerp* xOutField = xOutMesh->registeredFieldPerpIds[id];
        for (int i = 0; i < nout * LocalNz; i++) {
          IndPerp ind((xend + 1) * LocalNz + i, 1, LocalNz);
          f[ind] = (*xOutField)[ind.xm(xend - xstart + 1)];
        }
      }
    }
  }

  virtual int localSize3D() override { return local3D; }

  virtual int localSize2D() override { return local2D; }

  virtual int localSizePerp() override { return localPerp; }

  virtual int globalStartIndex3D() override { return start3D; }

  virtual int globalStartIndex2D() override { return start2D; }

  virtual int globalStartIndexPerp() override { return startPerp; }

  /// Use these methods to let the mesh know that this field has been
  /// created with it. It can then check in with its sibling meshes
  /// (representing other processors) to see if a corresponding field
  /// has been created for them which can be used to communicate guard
  /// cells with.
  void registerField(FieldData& f, int id) {
    registeredFields.emplace(&f, id);
    registeredFieldIds.emplace(id, &f);
  }
  void registerField(FieldPerp& f, int id) {
    registeredFieldPerps.emplace(&f, id);
    registeredFieldPerpIds.emplace(id, &f);
  }

  friend std::vector<FakeParallelMesh> createFakeProcessors(int nx, int ny, int nz,
                                                            int nxpe, int nype);

private:
  FakeParallelMesh *yUpMesh, *yDownMesh, *xInMesh, *xOutMesh;
  std::map<FieldData*, int> registeredFields;
  std::map<int, FieldData*> registeredFieldIds;
  std::map<FieldPerp*, int> registeredFieldPerps;
  std::map<int, FieldPerp*> registeredFieldPerpIds;

  int local3D, local2D, localPerp;
  int start3D, start2D, startPerp;
  int wait_any_count;

  comm_handle parentSend(FieldGroup& g) { return BoutMesh::send(g); }
  virtual int MPI_Irecv(void* UNUSED(buf), int UNUSED(count),
                        MPI_Datatype UNUSED(datatype), int UNUSED(source),
                        int UNUSED(tag), MPI_Comm UNUSED(comm),
                        MPI_Request* UNUSED(request)) override {
    return 0;
  }
  virtual int MPI_Isend(const void* UNUSED(buf), int UNUSED(count),
                        MPI_Datatype UNUSED(datatype), int UNUSED(dest), int UNUSED(tag),
                        MPI_Comm UNUSED(comm), MPI_Request* UNUSED(request)) override {
    return 0;
  }
  virtual int MPI_Send(const void* UNUSED(buf), int UNUSED(count),
                       MPI_Datatype UNUSED(datatype), int UNUSED(dest), int UNUSED(tag),
                       MPI_Comm UNUSED(comm)) override {
    return 0;
  }
  virtual int MPI_Wait(MPI_Request* UNUSED(request),
                       MPI_Status* UNUSED(status)) override {
    return 0;
  }
  virtual int MPI_Waitany(int UNUSED(count), MPI_Request UNUSED(array_of_requests[]),
                          int* indx, MPI_Status* UNUSED(status)) override {
    if (yUpMesh && wait_any_count < 0 && UpXSplitIndex() > 0) {
      *indx = wait_any_count = 0;
    } else if (yDownMesh && wait_any_count < 1 && UpXSplitIndex() == 0) {
      *indx = wait_any_count = 1;
    } else if (yDownMesh && wait_any_count < 2 && DownXSplitIndex() > 0) {
      *indx = wait_any_count = 2;
    } else if (yDownMesh && wait_any_count < 3 && DownXSplitIndex() == 0) {
      *indx = wait_any_count = 3;
    } else if (xInMesh && wait_any_count < 4) {
      *indx = wait_any_count = 4;
    } else if (xOutMesh && wait_any_count < 5) {
      *indx = wait_any_count = 5;
    } else {
      *indx = MPI_UNDEFINED;
      wait_any_count = -1;
    }
    return 0;
  }

  FieldGroup makeGroup(FakeParallelMesh* m, const std::vector<int> ids) {
    FieldGroup g;
    for (int i : ids) {
      ASSERT1(m->registeredFieldIds.count(i) != 0);
      g.add(*m->registeredFieldIds[i]);
    }
    return g;
  }
};

std::vector<FakeParallelMesh> createFakeProcessors(int nx, int ny, int nz, int nxpe,
                                                   int nype) {
  std::shared_ptr<Coordinates> test_coords{nullptr};
  std::vector<FakeParallelMesh> meshes;
  for (int i = 0; i < nxpe; i++) {
    for (int j = 0; j < nype; j++) {
      meshes.push_back(FakeParallelMesh(nx, ny, nz, nxpe, nype, i, j));
      meshes[j + i * nype].createDefaultRegions();
      bout::globals::mesh = &meshes[j + i * nype];
      static_cast<FakeParallelMesh*>(bout::globals::mesh)->setCoordinates(nullptr);
      test_coords = std::make_shared<Coordinates>(
          bout::globals::mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0},
          Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
          Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
          Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, false);
      static_cast<FakeParallelMesh*>(&meshes[j + i * nype])->setCoordinates(test_coords);
      test_coords->setParallelTransform(
          bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));
    }
  }
  int start3 = 0, start2 = 0, startP = 0;
  for (int i = 0; i < nxpe; i++) {
    for (int j = 0; j < nype; j++) {
      meshes.at(j + i * nype).local3D = (nx - 2) * ny * nz;
      meshes.at(j + i * nype).local2D = (nx - 2) * ny;
      meshes.at(j + i * nype).localPerp = (nx - 2) * nz;
      if (i == 0) {
        meshes.at(j + i * nype).local3D += ny * nz;
        meshes.at(j + i * nype).local2D += ny;
        meshes.at(j + i * nype).localPerp += nz;
      } else {
        meshes.at(j + i * nype).xInMesh = &meshes.at(j + (i - 1) * nype);
      }
      if (i == nxpe - 1) {
        meshes.at(j + i * nype).local3D += ny * nz;
        meshes.at(j + i * nype).local2D += ny;
        meshes.at(j + i * nype).localPerp += nz;
      } else {
        meshes.at(j + i * nype).xOutMesh = &meshes.at(j + (i + 1) * nype);
      }
      if (j == 0) {
        meshes.at(j + i * nype).local3D += (nx - 2) * nz;
      } else {
        meshes.at(j + i * nype).yDownMesh = &meshes.at(j - 1 + i * nype);
      }
      if (j == nype - 1) {
        meshes.at(j + i * nype).local3D += (nx - 2) * nz;
      } else {
        meshes.at(j + i * nype).yUpMesh = &meshes.at(j + 1 + i * nype);
      }
      meshes.at(j + i * nype).start3D = start3;
      meshes.at(j + i * nype).start2D = start2;
      meshes.at(j + i * nype).startPerp = startP;
      start3 += meshes.at(j + i * nype).local3D;
      start2 += meshes.at(j + i * nype).local2D;
      startP += meshes.at(j + i * nype).localPerp;
    }
    meshes.at(nype - 1 + i * nype).yUpMesh = &meshes.at(i * nype);
    meshes.at(i * nype).yDownMesh = &meshes.at(nype - 1 + i * nype);
  }
  return meshes;
}

#endif // FAKE_PARALLEL_MESH_H__
