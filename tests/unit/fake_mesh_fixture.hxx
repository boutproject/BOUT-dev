#pragma once

#include <gtest/gtest.h>
#include <memory>

#include "bout/build_config.hxx"
#include "bout/metric_tensor.hxx"
#include <bout/boundary_region.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutcomm.hxx>
#include <bout/build_defines.hxx>
#include <bout/coordinates.hxx>
#include <bout/field2d.hxx>
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
template <int NX, int NY, int NZ, bool FCI = false>
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
    mesh_m.setCoordinates(test_coords);

    // Set nonuniform corrections
    test_coords->setNon_uniform(true);
    test_coords->setD1_dx(0.2);
    test_coords->setD1_dy(0.2);
    test_coords->setD1_dz(0.0);

    if (bout::build::use_metric_3d) {
      bout::FieldMetric mutable_Bxy = test_coords->Bxy();
      mutable_Bxy.splitParallelSlices();
      mutable_Bxy.yup() = test_coords->Bxy();
      mutable_Bxy.ydown() = test_coords->Bxy();
      test_coords->setBxy(mutable_Bxy);
    }

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
    mesh_staggered_m.setCoordinates(test_coords_staggered);

    // Set nonuniform corrections
    test_coords_staggered->setNon_uniform(true);
    test_coords_staggered->setD1_dx(0.2);
    test_coords_staggered->setD1_dy(0.2);
    test_coords_staggered->setD1_dz(0.0);

    if (bout::build::use_metric_3d) {
      bout::FieldMetric mutable_Bxy = test_coords_staggered->Bxy();
      mutable_Bxy.splitParallelSlices();
      test_coords_staggered->setBxy(mutable_Bxy);

      mutable_Bxy = test_coords_staggered->Bxy();
      mutable_Bxy.yup() = test_coords_staggered->Bxy();
      mutable_Bxy.ydown() = test_coords_staggered->Bxy();
      test_coords_staggered->setBxy(mutable_Bxy);
    }

    test_coords_staggered->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(mesh_staggered_m));

    // Set all coordinates to the same Coordinates object for now
    mesh_staggered_m.setCoordinates(test_coords_staggered);
    mesh_staggered_m.setCoordinates(test_coords_staggered, CELL_XLOW);
    mesh_staggered_m.setCoordinates(test_coords_staggered, CELL_YLOW);
    mesh_staggered_m.setCoordinates(test_coords_staggered, CELL_ZLOW);

    if constexpr (FCI) {
      mesh_m.getCoordinates()->setParallelTransform(
          bout::utils::make_unique<MockParallelTransform>(mesh_m, false));
    }
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
using FakeMeshFixtureFCI = FakeMeshFixture_tmpl<3, 5, 7, true>;
