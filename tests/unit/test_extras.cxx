#include "test_extras.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "fieldperp.hxx"
#include "gtest/gtest.h"

#include <cmath>

// Need to provide a redundant declaration because C++
constexpr int FakeMeshFixture::nx;
constexpr int FakeMeshFixture::ny;
constexpr int FakeMeshFixture::nz;

::testing::AssertionResult IsSubString(const std::string &str,
                                       const std::string &substring) {
  if (str.find(substring) != std::string::npos) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure() << '"' << substring << "\" not found in " << str;
  }
}

::testing::AssertionResult IsField3DEqualBoutReal(const Field3D &field, BoutReal number,
                                                  BoutReal tolerance) {
  const auto &region = field.getMesh()->getRegion3D("RGN_ALL");
  BOUT_FOR_SERIAL(i, region) {
    if (fabs(field[i] - number) > tolerance) {
      return ::testing::AssertionFailure()
             << "Field3D(" << i.x() << ", " << i.y() << ", " << i.z()
             << ") == " << field[i] << "; Expected: " << number;
    }
  }

  return ::testing::AssertionSuccess();
}

::testing::AssertionResult IsField2DEqualBoutReal(const Field2D &field, BoutReal number,
                                                  BoutReal tolerance) {
  const auto &region = field.getMesh()->getRegion2D("RGN_ALL");
  BOUT_FOR_SERIAL(i, region) {
    if (fabs(field[i] - number) > tolerance) {
      return ::testing::AssertionFailure()
             << "Field2D(" << i.x() << ", " << i.y() << ") == " << field[i]
             << "; Expected: " << number;
    }
  }

  return ::testing::AssertionSuccess();
}

::testing::AssertionResult IsFieldPerpEqualBoutReal(const FieldPerp &field,
                                                    BoutReal number, BoutReal tolerance) {
  const auto &region = field.getMesh()->getRegionPerp("RGN_ALL");
  BOUT_FOR_SERIAL(i, region) {
    if (fabs(field[i] - number) > tolerance) {
      return ::testing::AssertionFailure()
             << "FieldPerp(" << i.x() << ", " << i.z() << ") == " << field[i]
             << "; Expected: " << number;
    }
  }

  return ::testing::AssertionSuccess();
}
