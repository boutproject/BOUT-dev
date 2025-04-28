#include "gtest/gtest.h"
#include "bout/coordinates.hxx"
#include "bout/mesh.hxx"
#include "fake_mesh_fixture.hxx"
#include <bout/tokamak_coordinates.hxx>


using bout::globals::mesh;

class CoordinateTransformTest : public FakeMeshFixture {
public:
    using FieldMetric = Coordinates::FieldMetric;

    CoordinateTransformTest() : FakeMeshFixture() {}
};

TEST_F(CoordinateTransformTest, CylindricalToCartesian) {

    auto tokamak_options = bout::TokamakOptions(*mesh);

    for (int i = 0; i < tokamak_options.Rxy.getNx(); i++) {
        for (int j = 0; j < tokamak_options.Rxy.getNy(); j++) {
            tokamak_options.Rxy(i, j) = ((float) i + 1) / 1000 * ((float) j + 1) / 1000;
        }
    }

    bout::Coordinates3D cartesian_coords = tokamak_options.CylindricalCoordinatesToCartesian();

    for (int jx = 0; jx < mesh->xstart; jx++) {
        for (int jy = 0; jy < mesh->ystart; jy++) {
            for (int jz = 0; jz < mesh->LocalNz; jz++) {

                auto actual_x = cartesian_coords.x(jx, jy, jz);
                auto actual_y = cartesian_coords.y(jx, jy, jz);
                auto actual_z = cartesian_coords.z(jx, jy, jz);

                auto expected_x = tokamak_options.Rxy(jx, jy) * cos(tokamak_options.toroidal_angle(jx, jy, jz));
                auto expected_y = tokamak_options.Rxy(jx, jy) * sin(tokamak_options.toroidal_angle(jx, jy, jz));
                auto expected_z = tokamak_options.Zxy(jx, jy);

                EXPECT_EQ(actual_x, expected_x);
                EXPECT_EQ(actual_y, expected_y);
                EXPECT_EQ(actual_z, expected_z);
            }
        }
    }
}
