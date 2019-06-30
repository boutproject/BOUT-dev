/*
 * Implements the identity transform for parallel derivatives
 * 
 * By default fields are stored so that Y is field aligned, so X-Z are not
 * orthogonal.
 *
 */

#include "bout/mesh.hxx"
#include "bout/paralleltransform.hxx"

void ParallelTransformIdentity::calcParallelSlices(Field3D& f) {
  if (f.getDirectionY() == YDirectionType::Aligned) {
    // Cannot calculate parallel slices for field-aligned fields, so just return without
    // setting yup or ydown
    return;
  }

  f.splitParallelSlices();

  for (int i = 0; i < f.getMesh()->ystart; ++i) {
    f.yup(i) = f;
    f.ydown(i) = f;
  }
}

void ParallelTransformIdentity::checkInputGrid() {
  std::string parallel_transform;
  if (!mesh.get(parallel_transform, "parallel_transform")) {
    if (parallel_transform != "identity") {
      throw BoutException("Incorrect parallel transform type '"+parallel_transform+"' used "
          "to generate metric components for ParallelTransformIdentity. Should be "
          "'identity'.");
    }
  } // else: parallel_transform variable not found in grid input, indicates older input
    //       file so must rely on the user having ensured the type is correct
}

void ParallelTransformIdentity::applyTwistShift(Field3D& f, bool twist_shift_enabled) {
  // All Field3Ds require twist-shift, because all are effectively field-aligned, but
  // allow twist-shift to be turned off by twist_shift_enabled
  if (twist_shift_enabled) {
    // Lower boundary
    // Note "TwistShiftDown" Region is empty if no twist-shift is required at the lower
    // boundary of this processor
    BOUT_FOR(i, f.getRegion2D("TwistShiftDown")) {
      int x = i.x();
      ParallelTransform::shiftZ(f, x, i.y(), ShiftAngle[x]);
    }

    // Upper boundary
    // Note "TwistShiftUp" Region is empty if no twist-shift is required at the upper
    // boundary of this processor
    BOUT_FOR(i, f.getRegion2D("TwistShiftUp")) {
      int x = i.x();
      ParallelTransform::shiftZ(f, x, i.y(), -ShiftAngle[x]);
    }
  }
}
