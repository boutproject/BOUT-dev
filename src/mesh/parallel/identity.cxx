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
  std::string coordinates_type = "";
  if (!mesh.get(coordinates_type, "coordinates_type")) {
    if (coordinates_type != "field_aligned") {
      throw BoutException("Incorrect coordinate system type '"+coordinates_type+"' used "
          "to generate metric components for ParallelTransformIdentity. Should be "
          "'field_aligned'.");
    }
  } // else: coordinate_system variable not found in grid input, indicates older input
    //       file so must rely on the user having ensured the type is correct
}
