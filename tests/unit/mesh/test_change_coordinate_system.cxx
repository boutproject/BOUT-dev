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

    const double R0 = 2.0;  // major radius
    const std::array<double, 5> r_values = {0.10, 0.15, 0.20, 0.25, 0.30};  // minor radius
    const std::array<double, 4> theta_values = {0.0, 1.07712, 3.17151, 5.26591};  // poloidal angle

    auto tokamak_options = bout::TokamakOptions(*mesh);

    int i = 0;
    for (auto r: r_values) {
        int j = 0;
        for (auto theta: theta_values) {
            tokamak_options.Rxy(i, j) = R0 + r * std::cos(theta);
            tokamak_options.Zxy(i, j) = r * std::sin(theta);
            j++;
        }
        i++;
    }

    bout::Coordinates3D cartesian_coords = tokamak_options.CylindricalCoordinatesToCartesian();

    for (int jx = 0; jx < mesh->xend; jx++) {
        for (int jy = 0; jy < mesh->yend; jy++) {
            for (int jz = 0; jz < mesh->LocalNz; jz++) {

                auto actual_x = cartesian_coords.x(jx, jy, jz);
                auto actual_y = cartesian_coords.y(jx, jy, jz);
                auto actual_z = cartesian_coords.z(jx, jy, jz);

                auto expected_x = tokamak_options.Rxy(jx, jy) * std::cos(tokamak_options.toroidal_angles[jz]);
                auto expected_y = tokamak_options.Rxy(jx, jy) * std::sin(tokamak_options.toroidal_angles[jz]);
                auto expected_z = tokamak_options.Zxy(jx, jy);

                EXPECT_EQ(actual_x, expected_x);
                EXPECT_EQ(actual_y, expected_y);
                EXPECT_EQ(actual_z, expected_z);
            }
        }
    }
}
