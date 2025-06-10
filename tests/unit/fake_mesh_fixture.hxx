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
class FakeMeshFixture : public ::testing::Test {
public:
    FakeMeshFixture(int nx_ = nx, int ny_ = ny, int nz_ = nz) {
        WithQuietOutput quiet_info{output_info};
        WithQuietOutput quiet_warn{output_warn};

        delete bout::globals::mesh;
        bout::globals::mpi = new MpiWrapper();
        bout::globals::mesh = new FakeMesh(nx_, ny_, nz_);
        bout::globals::mesh->createDefaultRegions();
        static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(nullptr);
        test_coords = std::make_shared<Coordinates>(
                bout::globals::mesh, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
                Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
                Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
                Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0});
        static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(test_coords);

        // Set nonuniform corrections
        test_coords->setNon_uniform(true);
        test_coords->setD1_dx(0.2);
        test_coords->setD1_dy(0.2);
        test_coords->setD1_dz(0.0);
#if BOUT_USE_METRIC_3D

        FieldMetric mutable_Bxy = test_coords->Bxy();
        mutable_Bxy.splitParallelSlices();
        test_coords->setBxy(mutable_Bxy);

        mutable_Bxy = test_coords->Bxy();
        mutable_Bxy.yup() = test_coords->Bxy();
        mutable_Bxy.ydown() = test_coords->Bxy();
        test_coords->setBxy(mutable_Bxy);

#endif

        static_cast<FakeMesh*>(bout::globals::mesh)->setCoordinates(test_coords);
        static_cast<FakeMesh*>(bout::globals::mesh)
                ->setGridDataSource(new FakeGridDataSource());
        // May need a ParallelTransform to create fields, because create3D calls
        // fromFieldAligned
        test_coords->setParallelTransform(
                bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));
        dynamic_cast<FakeMesh*>(bout::globals::mesh)->createBoundaryRegions();

        delete mesh_staggered;
        mesh_staggered = new FakeMesh(nx_, ny_, nz_);
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
        static_cast<FakeMesh*>(mesh_staggered)->setCoordinates(test_coords_staggered);

        // Set nonuniform corrections
        test_coords_staggered->setNon_uniform(true);
        test_coords_staggered->setD1_dx(0.2);
        test_coords_staggered->setD1_dy(0.2);
        test_coords_staggered->setD1_dz(0.0);
#if BOUT_USE_METRIC_3D

        mutable_Bxy = test_coords_staggered->Bxy();
        mutable_Bxy.splitParallelSlices();
        test_coords_staggered->setBxy(mutable_Bxy);

        mutable_Bxy = test_coords_staggered->Bxy();
        mutable_Bxy.yup() = test_coords_staggered->Bxy();
        mutable_Bxy.ydown() = test_coords_staggered->Bxy();
        test_coords_staggered->setBxy(mutable_Bxy);

#endif

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