/*
 * Implements the identity transform for parallel derivatives
 * 
 * By default fields are stored so that Y is field aligned, so X-Z are not
 * orthogonal.
 *
 */

#include <bout/mesh.hxx>
#include <bout/paralleltransform.hxx>

void ParallelTransformIdentity::checkInputGrid() {
  std::string coordinates_type = "";
  if (mesh.get(coordinates_type, "coordinates_type")) {
    // coordinate_system variable not found in grid input
    return;
  } else {
    if (coordinates_type != "field_aligned") {
      throw BoutException("Incorrect coordinate system type "+coordinates_type+" used "
          "to generate metric components for ParallelTransformIdentity. Should be "
          "'field_aligned.");
    }
  }
}
