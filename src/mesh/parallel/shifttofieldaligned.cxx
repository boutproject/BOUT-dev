/*
 * Implements the shifted metric method for parallel derivatives by
 * transforming to globally field-aligned coordinates (in contrast to
 * ShiftedMetric where the coordinates are only locally field-aligned).
 *
 * By default fields are stored so that X-Z are orthogonal,
 * and so not aligned in Y.
 *
 */

#include "shifttofieldaligned.hxx"

#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <fft.hxx>
#include <interpolation.hxx>

namespace {
void fixzShiftBoundaries(Field2D& zShift) {
  Mesh* localmesh = zShift.getMesh();

  localmesh->communicate(zShift);

  for (int x = 0; x < localmesh->LocalNx; x++) {
    if (localmesh->hasBranchCutDown(x)) {
      // ShiftAngle is the total toroidal shift of a field line after one poloidal turn
      BoutReal ShiftAngle = 0.;
      localmesh->periodicY(x, ShiftAngle);
      for (int y = 0; y < localmesh->ystart; y++) {
        zShift(x, y) -= ShiftAngle;
      }
    }
    if (localmesh->hasBranchCutUp(x)) {
      // ShiftAngle is the total toroidal shift of a field line after one poloidal turn
      BoutReal ShiftAngle = 0.;
      localmesh->periodicY(x, ShiftAngle);
      for (int y = localmesh->yend + 1; y < localmesh->LocalNy; y++) {
        zShift(x, y) += ShiftAngle;
      }
    }
  }

  zShift.applyBoundary("free_o3");
}
} // namespace

ShiftToFieldAligned::ShiftToFieldAligned(Mesh& mesh_in) {
  // Must *not* twist-shift when communicating
  bool twistshift = Options::root()["TwistShift"].withDefault(false);
  if (twistshift) {
    throw BoutException("ShiftToFieldAligned requires TwistShift=false");
  }

  Field2D zShift(&mesh_in);
  // Read the zShift angle from the mesh
  if (mesh_in.get(zShift, "zShift")) {
    // No zShift variable. Try qinty in BOUT grid files
    mesh_in.get(zShift, "qinty");
  }

  init(mesh_in, zShift);
}

void ShiftToFieldAligned::init(Mesh& mesh_in, const Field2D& zShift) {
  // always include CELL_CENTRE implementation
  implementations.emplace(CELL_CENTRE, Implementation(mesh_in, CELL_CENTRE, zShift));

  if (mesh_in.StaggerGrids) {
    // Don't know which locations we will need, so populate all those that can
    // be interpolated (for now at least)
    implementations.emplace(CELL_XLOW, Implementation(mesh_in, CELL_XLOW, zShift));
    implementations.emplace(CELL_XLOW, Implementation(mesh_in, CELL_YLOW, zShift));
    implementations.emplace(CELL_XLOW, Implementation(mesh_in, CELL_ZLOW, zShift));
  }
}

ShiftToFieldAligned::Implementation::Implementation(Mesh& mesh_in,
                                                    const CELL_LOC location_in,
                                                    const Field2D& zShift_in)
    : mesh(mesh_in), location(location_in), zShift(&mesh_in) {

  zShift = interp_to(zShift_in, location);
  fixzShiftBoundaries(zShift);
  cachePhases();
}

void ShiftToFieldAligned::Implementation::cachePhases() {
  // If we wanted to be efficient we could move the following cached phase setup
  // into the relevant shifting routines (with a bool to record if the phases
  // have already been cached) so that we only calculate the phase if we
  // actually call a relevant shift routine -- however that would make calls to
  // toFieldAligned and fromFieldAligned non-const, so it's better to keep this
  // in the constructor

  // As we're attached to a mesh we can expect the z direction to
  // not change once we've been created so precalculate the complex
  // phases used in transformations
  nmodes = mesh.LocalNz / 2 + 1;
  BoutReal zlength = mesh.getCoordinates()->zlength();

  // Allocate storage for our 3d vector structures.
  fromAlignedPhs = Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, mesh.LocalNz);
  toAlignedPhs = Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, mesh.LocalNz);

  // To/From field aligned phases
  BOUT_FOR(i, mesh.getRegion2D("RGN_ALL")) {
    for (int jz = 0; jz < nmodes; jz++) {
      BoutReal kwave = jz * 2.0 * PI / zlength; // wave number is 1/[rad]
      fromAlignedPhs(i.x(), i.y(), jz) =
        dcomplex(cos(kwave * zShift[i]), -sin(kwave * zShift[i]));
      toAlignedPhs(i.x(), i.y(), jz) =
        dcomplex(cos(kwave * zShift[i]), sin(kwave * zShift[i]));
    }
  }
}

/*!
 * Get the shifted field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftToFieldAligned::Implementation::toFieldAligned(const Field3D& f,
                                                                  const REGION region) {
  return shiftZ(f, toAlignedPhs, region);
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftToFieldAligned::Implementation::fromFieldAligned(const Field3D& f,
                                                                    const REGION region) {
  return shiftZ(f, fromAlignedPhs, region);
}

const Field3D ShiftToFieldAligned::Implementation::shiftZ(
    const Field3D& f, const Tensor<dcomplex>& phs, const REGION region) {
  ASSERT1(&mesh == f.getMesh());
  ASSERT1(f.getLocation() == location);

  if (mesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&mesh);
  result.allocate();
  result.setLocation(location);

  BOUT_FOR(i, mesh.getRegion2D(REGION_STRING(region))) {
    shiftZ(&f(i, 0), &phs(i.x(), i.y(), 0), &result(i, 0));
  }

  return result;
}

void ShiftToFieldAligned::Implementation::shiftZ(const BoutReal* in,
                                                 const dcomplex* phs,
                                                 BoutReal* out) {
  Array<dcomplex> cmplx(nmodes);

  // Take forward FFT
  rfft(in, mesh.LocalNz, &cmplx[0]);

  // Following is an algorithm approach to write a = a*b where a and b are
  // vectors of dcomplex.
  //  std::transform(cmplxOneOff.begin(),cmplxOneOff.end(), ptr.begin(),
  //                 cmplxOneOff.begin(), std::multiplies<dcomplex>());

  for (int jz = 1; jz < nmodes; jz++) {
    cmplx[jz] *= phs[jz];
  }

  irfft(&cmplx[0], mesh.LocalNz, out); // Reverse FFT
}
