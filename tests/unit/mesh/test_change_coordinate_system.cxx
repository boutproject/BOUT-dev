#include "gtest/gtest.h"
#include "bout/coordinates.hxx"
#include "bout/mesh.hxx"
#include "fake_mesh_fixture.hxx"
#include "bout/constants.hxx"
#include <bout/tokamak_coordinates.hxx>


using bout::globals::mesh;

class CoordinateTransformTest : public FakeMeshFixture {
public:
    using FieldMetric = Coordinates::FieldMetric;

    CoordinateTransformTest() : FakeMeshFixture() {}
};

TEST_F(CoordinateTransformTest, CylindricalToCartesian) {

    double R0 = 2.0;  // major radius
    double r[9] = {0.10, 0.15, 0.20, 0.25, 0.30};  // minor radius
    double theta[3] = {1.07712, 3.17151, 5.26591};  // poloidal angle

    auto tokamak_options = bout::TokamakOptions(*mesh);

    for (int i = 0; i < tokamak_options.Rxy.getNx(); i++) {
        for (int j = 0; j < tokamak_options.Zxy.getNy(); j++) {
            tokamak_options.Rxy(i, j) = R0 + r[i] * std::cos(theta[j]);
            tokamak_options.Zxy(i, j) = r[i] * std::sin(theta[j]);
        }
    }

    bout::Coordinates3D cartesian_coords = tokamak_options.CylindricalCoordinatesToCartesian();

    for (int jx = 0; jx < mesh->xend; jx++) {
        for (int jy = 0; jy < mesh->yend; jy++) {
            for (int jz = 0; jz < mesh->LocalNz; jz++) {

                auto actual_x = cartesian_coords.x(jx, jy, jz);
                auto actual_y = cartesian_coords.y(jx, jy, jz);
                auto actual_z = cartesian_coords.z(jx, jy, jz);

                auto expected_x = tokamak_options.Rxy(jx, jy) * std::cos(tokamak_options.toroidal_angle(jx, jy, jz));
                auto expected_y = tokamak_options.Rxy(jx, jy) * std::sin(tokamak_options.toroidal_angle(jx, jy, jz));
                auto expected_z = tokamak_options.Zxy(jx, jy);

                EXPECT_EQ(actual_x, expected_x);
                EXPECT_EQ(actual_y, expected_y);
                EXPECT_EQ(actual_z, expected_z);
            }
        }
    }
}
