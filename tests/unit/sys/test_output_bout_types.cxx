#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"
#include "gtest/gtest.h"

#include "bout/field2d.hxx"
#include "bout/globals.hxx"
#include "bout/output_bout_types.hxx"

#include <string>

#include <fmt/format.h>

using FormatFieldTest = FakeMeshFixture_tmpl<3, 5, 2>;

TEST_F(FormatFieldTest, Field2D) {
  Field2D f{bout::globals::mesh};

  fillField(f, {{0., 1., 2., 3., 4.}, {5., 6., 7., 8., 9.}, {10., 11., 12., 13., 14.}});

  const auto out = fmt::format("{}", f);

  const std::string expected =
      R"((0, 0): 0; (0, 1): 1; (0, 2): 2; (0, 3): 3; (0, 4): 4;
(1, 0): 5; (1, 1): 6; (1, 2): 7; (1, 3): 8; (1, 4): 9;
(2, 0): 10; (2, 1): 11; (2, 2): 12; (2, 3): 13; (2, 4): 14;

)";
  EXPECT_EQ(out, expected);
}

TEST_F(FormatFieldTest, Field2DSpec) {
  Field2D f{bout::globals::mesh};

  fillField(f, {{0., 1., 2., 3., 4.}, {5., 6., 7., 8., 9.}, {10., 11., 12., 13., 14.}});

  const auto out = fmt::format("{:3.1e}", f);

  const std::string expected =
      R"((0, 0): 0.0e+00; (0, 1): 1.0e+00; (0, 2): 2.0e+00; (0, 3): 3.0e+00; (0, 4): 4.0e+00;
(1, 0): 5.0e+00; (1, 1): 6.0e+00; (1, 2): 7.0e+00; (1, 3): 8.0e+00; (1, 4): 9.0e+00;
(2, 0): 1.0e+01; (2, 1): 1.1e+01; (2, 2): 1.2e+01; (2, 3): 1.3e+01; (2, 4): 1.4e+01;

)";
  EXPECT_EQ(out, expected);
}

TEST_F(FormatFieldTest, Field3D) {
  Field3D f{bout::globals::mesh};

  fillField(f, {{{0., 1}, {2., 3}, {4., 5}, {6., 7}, {8., 9}},
                {{10., 11}, {12., 13}, {14., 15}, {16., 17}, {18., 19}},
                {{20., 21}, {22., 23}, {24., 25}, {26., 27}, {28., 29}}});

  const auto out = fmt::format("{}", f);

  const std::string expected =
      R"((0, 0, 0): 0; (0, 0, 1): 1;
(0, 1, 0): 2; (0, 1, 1): 3;
(0, 2, 0): 4; (0, 2, 1): 5;
(0, 3, 0): 6; (0, 3, 1): 7;
(0, 4, 0): 8; (0, 4, 1): 9;

(1, 0, 0): 10; (1, 0, 1): 11;
(1, 1, 0): 12; (1, 1, 1): 13;
(1, 2, 0): 14; (1, 2, 1): 15;
(1, 3, 0): 16; (1, 3, 1): 17;
(1, 4, 0): 18; (1, 4, 1): 19;

(2, 0, 0): 20; (2, 0, 1): 21;
(2, 1, 0): 22; (2, 1, 1): 23;
(2, 2, 0): 24; (2, 2, 1): 25;
(2, 3, 0): 26; (2, 3, 1): 27;
(2, 4, 0): 28; (2, 4, 1): 29;


)";
  EXPECT_EQ(out, expected);
}

TEST_F(FormatFieldTest, Field3DSpec) {
  Field3D f{bout::globals::mesh};

  fillField(f, {{{0., 1}, {2., 3}, {4., 5}, {6., 7}, {8., 9}},
                {{10., 11}, {12., 13}, {14., 15}, {16., 17}, {18., 19}},
                {{20., 21}, {22., 23}, {24., 25}, {26., 27}, {28., 29}}});

  const auto out = fmt::format("{:3.1e}", f);

  const std::string expected =
      R"((0, 0, 0): 0.0e+00; (0, 0, 1): 1.0e+00;
(0, 1, 0): 2.0e+00; (0, 1, 1): 3.0e+00;
(0, 2, 0): 4.0e+00; (0, 2, 1): 5.0e+00;
(0, 3, 0): 6.0e+00; (0, 3, 1): 7.0e+00;
(0, 4, 0): 8.0e+00; (0, 4, 1): 9.0e+00;

(1, 0, 0): 1.0e+01; (1, 0, 1): 1.1e+01;
(1, 1, 0): 1.2e+01; (1, 1, 1): 1.3e+01;
(1, 2, 0): 1.4e+01; (1, 2, 1): 1.5e+01;
(1, 3, 0): 1.6e+01; (1, 3, 1): 1.7e+01;
(1, 4, 0): 1.8e+01; (1, 4, 1): 1.9e+01;

(2, 0, 0): 2.0e+01; (2, 0, 1): 2.1e+01;
(2, 1, 0): 2.2e+01; (2, 1, 1): 2.3e+01;
(2, 2, 0): 2.4e+01; (2, 2, 1): 2.5e+01;
(2, 3, 0): 2.6e+01; (2, 3, 1): 2.7e+01;
(2, 4, 0): 2.8e+01; (2, 4, 1): 2.9e+01;


)";
  EXPECT_EQ(out, expected);
}

TEST_F(FormatFieldTest, FieldPerp) {
  Field3D f{bout::globals::mesh};

  fillField(f, {{{0., 1}, {2., 3}, {4., 5}, {6., 7}, {8., 9}},
                {{10., 11}, {12., 13}, {14., 15}, {16., 17}, {18., 19}},
                {{20., 21}, {22., 23}, {24., 25}, {26., 27}, {28., 29}}});

  FieldPerp g = sliceXZ(f, 0);

  const auto out = fmt::format("{}", g);

  const std::string expected =
      R"((0, 0): 0; (0, 1): 1;
(1, 0): 10; (1, 1): 11;
(2, 0): 20; (2, 1): 21;

)";
  EXPECT_EQ(out, expected);
}

TEST_F(FormatFieldTest, FieldPerpSpec) {
  Field3D f{bout::globals::mesh};

  fillField(f, {{{0., 1}, {2., 3}, {4., 5}, {6., 7}, {8., 9}},
                {{10., 11}, {12., 13}, {14., 15}, {16., 17}, {18., 19}},
                {{20., 21}, {22., 23}, {24., 25}, {26., 27}, {28., 29}}});

  FieldPerp g = sliceXZ(f, 0);

  const auto out = fmt::format("{:3.1e}", g);

  const std::string expected =
      R"((0, 0): 0.0e+00; (0, 1): 1.0e+00;
(1, 0): 1.0e+01; (1, 1): 1.1e+01;
(2, 0): 2.0e+01; (2, 1): 2.1e+01;

)";
  EXPECT_EQ(out, expected);
}
