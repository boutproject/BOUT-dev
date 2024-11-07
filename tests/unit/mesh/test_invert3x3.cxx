#include "../../src/mesh/invert3x3.hxx"

#include "gtest/gtest.h"

TEST(Invert3x3Test, Identity) {
  Matrix<BoutReal> input(3, 3);
  input = 0;
  for (int i = 0; i < 3; i++) {
    input(i, i) = 1.0;
  }
  auto expected = input;
  bout::invert3x3(input);

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      EXPECT_EQ(input(i, j), expected(i, j));
    }
  }
}

TEST(Invert3x3Test, InvertTwice) {
  std::vector<BoutReal> rawDataMat = {0.05567105, 0.92458227, 0.19954631,
                                      0.28581972, 0.54009039, 0.13234403,
                                      0.8841194,  0.161224,   0.74853209};
  std::vector<BoutReal> rawDataInv = {-2.48021781, 4.27410022,  -0.09449605,
                                      0.6278449,   0.87275842,  -0.32168092,
                                      2.79424897,  -5.23628123, 1.51684677};

  Matrix<BoutReal> input(3, 3);
  Matrix<BoutReal> expected(3, 3);

  int counter = 0;
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      input(i, j) = rawDataMat[counter];
      expected(i, j) = rawDataInv[counter];
      counter++;
    }
  }

  // Invert twice to check if we get back to where we started
  bout::invert3x3(input);

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      // Note we only check to single tolerance here
      EXPECT_FLOAT_EQ(input(i, j), expected(i, j));
    }
  }
}

TEST(Invert3x3Test, Singular) {
  Matrix<BoutReal> input(3, 3);
  input = 0;
  auto result = bout::invert3x3(input);
  EXPECT_TRUE(result.has_value());
}

TEST(Invert3x3Test, BadCondition) {
  Matrix<BoutReal> input(3, 3);

  input = 0.;
  input(0, 0) = 1.0e-16;
  input(1, 1) = 1.0;
  input(2, 2) = 1.0;
  EXPECT_TRUE(bout::invert3x3(input).has_value());

  // not quite bad enough condition
  input = 0.;
  input(0, 0) = 1.0e-12;
  input(1, 1) = 1.0;
  input(2, 2) = 1.0;
  EXPECT_FALSE(bout::invert3x3(input).has_value());
}
