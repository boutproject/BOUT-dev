#ifndef FAKE_PARALLEL_MESH_H__
#define FAKE_PARALLEL_MESH_H__
#include <iostream>
#include "gtest/gtest.h"

#include <vector>
#include <numeric>
#include <functional>
#include <map>

#include "boutcomm.hxx"
#include "bout/mesh.hxx"
#include "bout/coordinates.hxx"
#include "bout/fieldgroup.hxx"
#include "../../src/mesh/impls/bout/boutmesh.hxx"
#include "fieldperp.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "unused.hxx"

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
  FakeParallelMesh(int nx, int ny, int nz, int nxpe, int nype, int pe_yind, int pe_xind) :
    BoutMesh(nx, ny, nz, 1, 1, nxpe, nype, pe_yind, pe_xind), yUpMesh(nullptr),
    yDownMesh(nullptr), xInMesh(nullptr), xOutMesh(nullptr)
  {
    StaggerGrids=false;
    periodicX = false;
    IncIntShear = false;
    maxregionblocksize = MAXREGIONBLOCKSIZE;
    calcParallelSlices_on_communicate = true;
    options = Options::getRoot();
  }

  void initDerivs(Options * opt){
    StaggerGrids=true;
    derivs_init(opt);
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

  comm_handle send(FieldGroup &g) override {
    for (FieldData *f : g.get()) {
      if (registeredFields.count(f) != 0) {
        int id = registeredFields[f];
        ASSERT1(registeredFieldIds[id] == f);
        if (yUpMesh != nullptr && yUpMesh->registeredFieldIds.count(id) != 0) {
	  FieldGroup yUpGroup(*yUpMesh->registeredFieldIds[id]);
          yUpMesh->parentSend(yUpGroup);
        }
        if (yDownMesh != nullptr && yDownMesh->registeredFieldIds.count(id) != 0) {
	  FieldGroup yDownGroup(*yDownMesh->registeredFieldIds[id]);
          yDownMesh->parentSend(yDownGroup);
        }
        if (xInMesh != nullptr && xInMesh->registeredFieldIds.count(id) != 0) {
	  FieldGroup xInGroup(*xInMesh->registeredFieldIds[id]);
          xInMesh->parentSend(xInGroup);
        }
        if (xOutMesh != nullptr && xOutMesh->registeredFieldIds.count(id) != 0) {
	  FieldGroup xOutGroup(*xOutMesh->registeredFieldIds[id]);
          xOutMesh->parentSend(xOutGroup);
        }
      }
    }
    return parentSend(g);
  }

  // Need to override this functions to trick mesh into communicating for
  // FieldPerp type
  void communicate(FieldPerp &f) override {
    int nin = xstart; // Number of x points in inner guard cell
    int nout = LocalNx-xend-1; // Number of x points in outer guard cell

    if (registeredFieldPerps.count(&f) != 0) {
      int id = registeredFieldPerps[&f];

      if (xInMesh != nullptr && xInMesh->registeredFieldPerpIds.count(id) != 0) {
        FieldPerp *xInField = xInMesh->registeredFieldPerpIds[id];
        for (int i = 0; i < nin*LocalNz; i++) {
          f[IndPerp(i)] = (*xInField)[IndPerp(xend-nout+1 + i)];
        }
      }

      if (xOutMesh != nullptr && xOutMesh->registeredFieldPerpIds.count(id) != 0) {
        FieldPerp *xOutField = xOutMesh->registeredFieldPerpIds[id];
        for (int i = 0; i < nin * LocalNz; i++) {
          f[IndPerp(xend + 1 + i)] = (*xOutField)[IndPerp(xstart + i - 1)];
        }
      }
    }
    
  }

  virtual int localSize3D() override {
    return local3D;
  }

  virtual int localSize2D() override {
    return local2D;
  }

  virtual int localSizePerp() override {
    return localPerp;
  }

  virtual int globalStartIndex3D() override {
    std::cout << "globalStartIndex3D";
    return start3D;
  }

  virtual int globalStartIndex2D() override {
    std::cout << "globalStartIndex2D";
    return start2D;
  }

  virtual int globalStartIndexPerp() override {
    std::cout << "globalStartIndexPerp";
    return startPerp;
  }

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
  std::map<FieldData*, int, std::less<void> > registeredFields;
  std::map<int, FieldData*> registeredFieldIds;
  std::map<FieldPerp*, int, std::less<void> > registeredFieldPerps;
  std::map<int, FieldPerp*> registeredFieldPerpIds;

  int local3D, local2D, localPerp;
  int start3D, start2D, startPerp;

  comm_handle parentSend(FieldGroup &g) {
    return BoutMesh::send(g);
  }
  virtual int MPI_Irecv(void *UNUSED(buf), int UNUSED(count), MPI_Datatype UNUSED(datatype),
                        int UNUSED(source), int UNUSED(tag), MPI_Comm UNUSED(comm),
                        MPI_Request *UNUSED(request)) override {
    return 0;
  }
  virtual int MPI_Isend(const void *UNUSED(buf), int UNUSED(count),
                        MPI_Datatype UNUSED(datatype), int UNUSED(dest), int UNUSED(tag),
                        MPI_Comm UNUSED(comm), MPI_Request *UNUSED(request)) override {
    return 0;
  }
  virtual int MPI_Send(const void *UNUSED(buf), int UNUSED(count),
                       MPI_Datatype UNUSED(datatype), int UNUSED(dest), int UNUSED(tag),
                       MPI_Comm UNUSED(comm)) override {
    return 0;
  }
  virtual int MPI_Wait(MPI_Request *UNUSED(request), MPI_Status *UNUSED(status)) override {
    return 0;
  }
  virtual int MPI_Waitany(int UNUSED(count), MPI_Request UNUSED(array_of_requests[]),
                          int *UNUSED(indx), MPI_Status *UNUSED(status)) override {
    return 0;
  }
};

std::vector<FakeParallelMesh> createFakeProcessors(int nx, int ny, int nz,
						   int nxpe, int nype) {
  std::vector<FakeParallelMesh> meshes;
  for (int i = 0; i < nxpe; i++) {
    for (int j = 0; j < nxpe; j++) {
      meshes.push_back(FakeParallelMesh(nx, ny, nz, nxpe, nype, i, j));
    }
  }
  int start3 = 0, start2 = 0, startP = 0;
  for (int j = 0; j < nxpe; j++) {
    for (int i = 0; i < nxpe; i++) {
      meshes[j*nype + i].local3D = (nx - 2) * (ny - 2) * nz;
      meshes[j*nype + i].local2D = (nx - 2) * (ny - 2);
      meshes[j*nype + i].localPerp = (nx - 2) * nz;
      if (i > 0) {
        meshes[j*nype + i].local3D += (ny - 2) * nz;
        meshes[j*nype + i].local2D += (ny - 2);
        meshes[j*nype + i].localPerp += (nx - 2) * nz;
        meshes[j*nype + i].xInMesh = &meshes[j*nype + i - 1];
      }
      if (i < nxpe - 1) {
        meshes[j*nype + i].local3D += (ny - 2) * nz;
        meshes[j*nype + i].local2D += (ny - 2);
        meshes[j*nype + i].localPerp += (nx - 2) * nz;
	meshes[j*nype + i].xOutMesh = &meshes[j*nype + i + 1];
      }
      if (j > 0) {
        meshes[j*nype + i].local3D += (nx - 2) * nz;
        meshes[j*nype + i].local2D += (nx - 2);
	meshes[j*nype + i].yUpMesh = &meshes[(j - 1)*nype + i];
      }
      if (j < nype - 1) {
        meshes[j*nype + i].local3D += (nx - 2) * nz;
        meshes[j*nype + i].local2D += (nx - 2);
	meshes[j*nype + i].yDownMesh = &meshes[(j + 1)*nype + i];
      }
      meshes[j*nype + i].start3D = start3;
      meshes[j*nype + i].start2D = start2;
      meshes[j*nype + i].startPerp = startP;
      start3 += meshes[j*nype + i].start3D;
      start2 += meshes[j*nype + i].start2D;
      startP += meshes[j*nype + i].startPerp;
    }
  }
  return meshes;
}


#endif // FAKE_PARALLEL_MESH_H__
