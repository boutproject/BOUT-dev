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

    // arrange

    // Set up test values
    // Calculate cylindrical coordinates (Rxy, Zxy)
    // from (2D) orthogonal poloidal coordinates (r, theta)

    const double R0 = 2.0;  // major radius
    const std::array<double, nx> r_values = {0.1, 0.2, 0.3};  // minor radius
    const std::array<double, ny> theta_values = {  // poloidal angle
            0.0,
            PI / 2,
            PI,
            3 * PI / 2,
            2 * PI};

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

    // act
    bout::Coordinates3D cartesian_coords = tokamak_options.CylindricalCoordinatesToCartesian();

    // assert
    const auto max_r = *std::max_element(begin(r_values), end(r_values));
    const auto expected_max_x = R0 + max_r;
    // With nz=7, there is no toroidal coordinate point at exactly pi/2; the nearest point is at 2/7 * 2pi
    const auto expected_max_y = (R0 + max_r) * std::sin(TWOPI * 2 / 7);
    const auto expected_max_z = max_r;

    // With nz=7, there is no toroidal coordinate point at exactly pi; the nearest point is at 3/7 * 2pi
    const auto expected_min_x = -1 * (R0 + max_r) * std::cos(TWOPI / 7 / 2);
    const auto expected_min_y = -1 * (R0 + max_r) * std::sin(TWOPI * 2 / 7);
    const auto expected_min_z = -1 * expected_max_z;

    const auto actual_max_x = max(cartesian_coords.x, false, "RGN_ALL");
    const auto actual_max_y = max(cartesian_coords.y, false, "RGN_ALL");
    const auto actual_max_z = max(cartesian_coords.z, false, "RGN_ALL");
    const auto actual_min_x = min(cartesian_coords.x, false, "RGN_ALL");
    const auto actual_min_y = min(cartesian_coords.y, false, "RGN_ALL");
    const auto actual_min_z = min(cartesian_coords.z, false, "RGN_ALL");

    EXPECT_EQ(expected_max_x, actual_max_x);
    EXPECT_EQ(expected_max_y, actual_max_y);
    EXPECT_EQ(expected_max_z, actual_max_z);
    EXPECT_EQ(expected_min_x, actual_min_x);
    EXPECT_EQ(expected_min_y, actual_min_y);
    EXPECT_EQ(expected_min_z, actual_min_z);
}
