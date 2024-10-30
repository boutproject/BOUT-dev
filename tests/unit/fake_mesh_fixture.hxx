#pragma once

#include <gtest/gtest.h>

#include <bout/boundary_region.hxx>
#include <bout/boutcomm.hxx>
#include <bout/coordinates.hxx>
#include <bout/griddata.hxx>
#include <bout/mesh.hxx>
#include <bout/mpi_wrapper.hxx>
#include <bout/operatorstencil.hxx>
#include <bout/unused.hxx>

#include "fake_mesh.hxx"

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
        bout::globals::mesh, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
        Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0});

    // Set some auxilliary variables
    // Usually set in geometry()
    // Note: For testing these are set to non-zero values
    test_coords->G1 = test_coords->G2 = test_coords->G3 = 0.1;

    // Set nonuniform corrections
    test_coords->non_uniform = true;
    test_coords->d1_dx = test_coords->d1_dy = 0.2;
    test_coords->d1_dz = 0.0;
#if BOUT_USE_METRIC_3D
    test_coords->Bxy.splitParallelSlices();
    test_coords->Bxy.yup() = test_coords->Bxy.ydown() = test_coords->Bxy;
#endif

    // No call to Coordinates::geometry() needed here
    static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(test_coords);
    static_cast<FakeMesh*>(bout::globals::mesh)
        ->setGridDataSource(new FakeGridDataSource());
    // May need a ParallelTransform to create fields, because create3D calls
    // fromFieldAligned
    test_coords->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));
    dynamic_cast<FakeMesh*>(bout::globals::mesh)->createBoundaryRegions();

    delete mesh_staggered;
    mesh_staggered = new FakeMesh(nx, ny, nz);
    mesh_staggered->StaggerGrids = true;
    dynamic_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr);
    dynamic_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_XLOW);
    dynamic_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_YLOW);
    dynamic_cast<FakeMesh*>(mesh_staggered)->setCoordinates(nullptr, CELL_ZLOW);
    mesh_staggered->createDefaultRegions();

    test_coords_staggered = std::make_shared<Coordinates>(
        mesh_staggered, Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{1.0, mesh_staggered}, Field2D{1.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{0.0, mesh_staggered}, Field2D{0.0, mesh_staggered},
        Field2D{0.0, mesh_staggered});

    // Set some auxilliary variables
    test_coords_staggered->G1 = test_coords_staggered->G2 = test_coords_staggered->G3 =
        0.1;

    // Set nonuniform corrections
    test_coords_staggered->non_uniform = true;
    test_coords_staggered->d1_dx = test_coords_staggered->d1_dy = 0.2;
    test_coords_staggered->d1_dz = 0.0;
#if BOUT_USE_METRIC_3D
    test_coords_staggered->Bxy.splitParallelSlices();
    test_coords_staggered->Bxy.yup() = test_coords_staggered->Bxy.ydown() =
        test_coords_staggered->Bxy;
#endif

    // No call to Coordinates::geometry() needed here
    test_coords_staggered->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*mesh_staggered));

    // Set all coordinates to the same Coordinates object for now
    dynamic_cast<FakeMesh*>(mesh_staggered)->setCoordinates(test_coords_staggered);
    dynamic_cast<FakeMesh*>(mesh_staggered)
        ->setCoordinates(test_coords_staggered, CELL_XLOW);
    dynamic_cast<FakeMesh*>(mesh_staggered)
        ->setCoordinates(test_coords_staggered, CELL_YLOW);
    dynamic_cast<FakeMesh*>(mesh_staggered)
        ->setCoordinates(test_coords_staggered, CELL_ZLOW);
  }

  FakeMeshFixture(const FakeMeshFixture&) = delete;
  FakeMeshFixture& operator=(const FakeMeshFixture&) = delete;
  FakeMeshFixture(FakeMeshFixture&&) = delete;
  FakeMeshFixture& operator=(FakeMeshFixture&&) = delete;

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
