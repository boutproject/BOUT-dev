#pragma once

#include <gtest/gtest.h>
#include <memory>

#include <bout/boundary_region.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutcomm.hxx>
#include <bout/build_defines.hxx>
#include <bout/coordinates.hxx>
#include <bout/globals.hxx>
#include <bout/griddata.hxx>
#include <bout/mesh.hxx>
#include <bout/mpi_wrapper.hxx>
#include <bout/operatorstencil.hxx>
#include <bout/output.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/utils.hxx>

#include "fake_mesh.hxx" // IWYU pragma: export

/// Test fixture to make sure the global mesh is our fake
/// one. Also initialize the global mesh_staggered for use in tests with
/// staggering. Multiple tests have exactly the same fixture, so use a type
/// alias to make a new test:
///
///     using MyTest = FakeMeshFixture;
///
/// Type alias `FakeMeshFixture = FakeMeshFixture_tmpl<3, 5, 7>`
/// is used as a shim to allow FakeMeshFixture to be used with default values for nx, ny, nz.
/// Use this template class directly to use different sized grid:
///
///     using MyTest = FakeMeshFixture_tmpl<7, 9, 11>;
template <int NX, int NY, int NZ>
class FakeMeshFixture_tmpl : public ::testing::Test {
public:
  FakeMeshFixture_tmpl()
      : mesh_m(NX, NY, NZ, mpi), mesh_staggered_m(NX, NY, NZ, mpi),
        mesh_staggered(&mesh_staggered_m) {

    bout::globals::mpi = &mpi;
    bout::globals::mesh = &mesh_m;
    bout::globals::mesh->createDefaultRegions();
    mesh_m.setCoordinates(nullptr);
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
    mesh_m.setCoordinates(test_coords);
    mesh_m.setGridDataSource(new FakeGridDataSource());
    // May need a ParallelTransform to create fields, because create3D calls
    // fromFieldAligned
    test_coords->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));
    mesh_m.createBoundaryRegions();

    mesh_staggered_m.StaggerGrids = true;
    mesh_staggered_m.setCoordinates(nullptr);
    mesh_staggered_m.setCoordinates(nullptr, CELL_XLOW);
    mesh_staggered_m.setCoordinates(nullptr, CELL_YLOW);
    mesh_staggered_m.setCoordinates(nullptr, CELL_ZLOW);
    mesh_staggered_m.createDefaultRegions();

    test_coords_staggered = std::make_shared<Coordinates>(
        &mesh_staggered_m, Field2D{1.0, &mesh_staggered_m},
        Field2D{1.0, &mesh_staggered_m}, Field2D{1.0, &mesh_staggered_m},
        Field2D{1.0, &mesh_staggered_m}, Field2D{1.0, &mesh_staggered_m},
        Field2D{1.0, &mesh_staggered_m}, Field2D{1.0, &mesh_staggered_m},
        Field2D{1.0, &mesh_staggered_m}, Field2D{0.0, &mesh_staggered_m},
        Field2D{0.0, &mesh_staggered_m}, Field2D{0.0, &mesh_staggered_m},
        Field2D{1.0, &mesh_staggered_m}, Field2D{1.0, &mesh_staggered_m},
        Field2D{1.0, &mesh_staggered_m}, Field2D{0.0, &mesh_staggered_m},
        Field2D{0.0, &mesh_staggered_m}, Field2D{0.0, &mesh_staggered_m},
        Field2D{0.0, &mesh_staggered_m}, Field2D{0.0, &mesh_staggered_m});

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
        bout::utils::make_unique<ParallelTransformIdentity>(mesh_staggered_m));

    // Set all coordinates to the same Coordinates object for now
    mesh_staggered_m.setCoordinates(test_coords_staggered);
    mesh_staggered_m.setCoordinates(test_coords_staggered, CELL_XLOW);
    mesh_staggered_m.setCoordinates(test_coords_staggered, CELL_YLOW);
    mesh_staggered_m.setCoordinates(test_coords_staggered, CELL_ZLOW);
  }

  FakeMeshFixture_tmpl(const FakeMeshFixture_tmpl&) = delete;
  FakeMeshFixture_tmpl& operator=(const FakeMeshFixture_tmpl&) = delete;
  FakeMeshFixture_tmpl(FakeMeshFixture_tmpl&&) = delete;
  FakeMeshFixture_tmpl& operator=(FakeMeshFixture_tmpl&&) = delete;

  ~FakeMeshFixture_tmpl() override {
    bout::globals::mesh = nullptr;
    bout::globals::mpi = nullptr;

    Options::cleanup();
  }

  static constexpr int nx = NX;
  static constexpr int ny = NY;
  static constexpr int nz = NZ;

private:
  std::shared_ptr<Coordinates> test_coords{nullptr};
  std::shared_ptr<Coordinates> test_coords_staggered{nullptr};

  WithQuietOutput quiet_info{output_info};
  WithQuietOutput quiet_warn{output_warn};
  MpiWrapper mpi;

  FakeMesh mesh_m;
  FakeMesh mesh_staggered_m;

public:
  // Public pointer to our staggered mesh
  Mesh* mesh_staggered; // NOLINT
};

using FakeMeshFixture = FakeMeshFixture_tmpl<3, 5, 7>;
