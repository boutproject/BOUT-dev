/*
 * Implements the shifted metric method for parallel derivatives
 *
 * By default fields are stored so that X-Z are orthogonal,
 * and so not aligned in Y.
 *
 */

#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include "bout/paralleltransform.hxx"
#include <fft.hxx>

#include <cmath>

#include <output.hxx>

ShiftedMetric::ShiftedMetric(Mesh& m, CELL_LOC location_in, Field2D zShift_,
    BoutReal zlength_in)
    : ParallelTransform(m), location(location_in), zShift(std::move(zShift_)),
      zlength(zlength_in) {
  ASSERT1(zShift.getLocation() == location);
  // check the coordinate system used for the grid data source
  checkInputGrid();

  // TwistShift should not be set for derivatives to be correct at the jump where
  // poloidal angle theta goes 2pi->0. zShift has been corrected for the jump
  // already in Coordinates::Coordinates
  bool twistshift = Options::root()["TwistShift"]
                        .doc("Enable twist-shift boundary condition in core region?")
                        .withDefault(false);
  if (twistshift) {
    throw BoutException("ShiftedMetric requires the option TwistShift=false");
  }

  cachePhases();
}

void ShiftedMetric::checkInputGrid() {
  std::string coordinates_type = "";
  if (!mesh.get(coordinates_type, "coordinates_type")) {
    if (coordinates_type != "orthogonal") {
      throw BoutException("Incorrect coordinate system type '" + coordinates_type
                          + "' used to generate metric components for ShiftedMetric. "
                            "Should be 'orthogonal'.");
    }
  } // else: coordinate_system variable not found in grid input, indicates older input
    //       file so must rely on the user having ensured the type is correct
}

void ShiftedMetric::cachePhases() {
  // If we wanted to be efficient we could move the following cached phase setup
  // into the relevant shifting routines (with static bool first protection)
  // so that we only calculate the phase if we actually call a relevant shift
  // routine -- however as we're only going to do this initialisation once I
  // think it's cleaner to put it in the constructor here.

  // As we're attached to a mesh we can expect the z direction to
  // not change once we've been created so precalculate the complex
  // phases used in transformations
  nmodes = (mesh.zend + 1 - mesh.zstart) / 2 + 1;

  // Allocate storage for our 3d phase information.
  fromAlignedPhs = Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, nmodes);
  toAlignedPhs = Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, nmodes);

  // To/From field aligned phases
  BOUT_FOR(i, mesh.getRegion2D("RGN_ALL")) {
    int ix = i.x();
    int iy = i.y();
    for (int jz = 0; jz < nmodes; jz++) {
      BoutReal kwave = jz * 2.0 * PI / zlength; // wave number is 1/[rad]
      fromAlignedPhs(ix, iy, jz) =
          dcomplex(cos(kwave * zShift[i]), -sin(kwave * zShift[i]));
      toAlignedPhs(ix, iy, jz) = dcomplex(cos(kwave * zShift[i]), sin(kwave * zShift[i]));
    }
  }

  // Allocate space for parallel slice caches: y-guard cells in each
  // direction
  parallel_slice_phases.resize(mesh.ystart * 2);

  // Careful with the indices/offsets! Offsets are 1-indexed (as 0
  // would be the original slice), and Mesh::ystart is the number of
  // guard cells. The parallel slice vector stores the offsets as
  //    {+1, ..., +n, -1, ..., -n}
  // Once parallel_slice_phases is initialised though, each element
  // stores its phase and offset, so we don't need to faff about after
  // this
  for (int i = 0; i < mesh.ystart; ++i) {
    // NOTE: std::vector constructor here takes a **copy** of the
    // Array! We *must* call `Array::ensureUnique` on each element
    // before using it!
    parallel_slice_phases[i].phase_shift =
        Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, nmodes);
    parallel_slice_phases[i].y_offset = i + 1;

    // Backwards parallel slices
    parallel_slice_phases[mesh.ystart + i].phase_shift =
        Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, nmodes);
    parallel_slice_phases[mesh.ystart + i].y_offset = -(i + 1);
  }

  // Parallel slice phases -- note we don't shift in the boundaries/guards
  for (auto& slice : parallel_slice_phases) {
    BOUT_FOR(i, mesh.getRegion2D("RGN_NOY")) {

      int ix = i.x();
      int iy = i.y();
      BoutReal slice_shift = zShift[i] - zShift[i.yp(slice.y_offset)];

      for (int jz = 0; jz < nmodes; jz++) {
        // wave number is 1/[rad]
        BoutReal kwave = jz * 2.0 * PI / zlength;

        slice.phase_shift(ix, iy, jz) =
            dcomplex(cos(kwave * slice_shift), -sin(kwave * slice_shift));
      }
    }
  }
}

/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftedMetric::toFieldAligned(const Field3D& f, const REGION region) {
  switch (f.getDirectionY()) {
  case (YDirectionType::Standard):
    return shiftZ(f, toAlignedPhs, YDirectionType::Aligned, region);
  case (YDirectionType::Aligned):
    // f is already in field-aligned coordinates
    return f;
  default:
    throw BoutException("Unrecognized y-direction type for Field3D passed to "
                        "ShiftedMetric::toFieldAligned");
    // This should never happen, but use 'return f' to avoid compiler warnings
    return f;
  }
}
const FieldPerp ShiftedMetric::toFieldAligned(const FieldPerp& f, const REGION region) {
  switch (f.getDirectionY()) {
  case (YDirectionType::Standard):
    return shiftZ(f, toAlignedPhs, YDirectionType::Aligned, region);
  case (YDirectionType::Aligned):
    // f is already in field-aligned coordinates
    return f;
  default:
    throw BoutException("Unrecognized y-direction type for FieldPerp passed to "
                        "ShiftedMetric::toFieldAligned");
    // This should never happen, but use 'return f' to avoid compiler warnings
    return f;
  }
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftedMetric::fromFieldAligned(const Field3D& f, const REGION region) {
  switch (f.getDirectionY()) {
  case (YDirectionType::Aligned):
    return shiftZ(f, fromAlignedPhs, YDirectionType::Standard, region);
  case (YDirectionType::Standard):
    // f is already in orthogonal coordinates
    return f;
  default:
    throw BoutException("Unrecognized y-direction type for Field3D passed to "
                        "ShiftedMetric::toFieldAligned");
    // This should never happen, but use 'return f' to avoid compiler warnings
    return f;
  }
}
const FieldPerp ShiftedMetric::fromFieldAligned(const FieldPerp& f, const REGION region) {
  switch (f.getDirectionY()) {
  case (YDirectionType::Aligned):
    return shiftZ(f, fromAlignedPhs, YDirectionType::Standard, region);
  case (YDirectionType::Standard):
    // f is already in orthogonal coordinates
    return f;
  default:
    throw BoutException("Unrecognized y-direction type for FieldPerp passed to "
                        "ShiftedMetric::toFieldAligned");
    // This should never happen, but use 'return f' to avoid compiler warnings
    return f;
  }
}

const Field3D ShiftedMetric::shiftZ(const Field3D& f, const Tensor<dcomplex>& phs,
                                    const YDirectionType y_direction_out,
                                    const REGION region) const {
  ASSERT1(f.getMesh() == &mesh);
  ASSERT1(f.getLocation() == location);

  if (mesh.zend == mesh.zstart)
    return f; // Shifting makes no difference

  Field3D result{emptyFrom(f).setDirectionY(y_direction_out)};

  BOUT_FOR(i, mesh.getRegion2D(toString(region))) {
    shiftZ(&f(i, mesh.zstart), &phs(i.x(), i.y(), mesh.zstart), &result(i, mesh.zstart));
  }

  return result;
}

const FieldPerp ShiftedMetric::shiftZ(const FieldPerp& f, const Tensor<dcomplex>& phs,
                                      const YDirectionType y_direction_out,
                                      const REGION UNUSED(region)) const {
  ASSERT1(f.getMesh() == &mesh);
  ASSERT1(f.getLocation() == location);

  if (mesh.zend == mesh.zstart)
    return f; // Shifting makes no difference

  FieldPerp result{emptyFrom(f).setDirectionY(y_direction_out)};

  const int y = f.getIndex();
  for (int i = mesh.xstart; i <= mesh.xend; ++i) {
    shiftZ(&f(i, mesh.zstart), &phs(i, y, mesh.zstart), &result(i, mesh.zstart));
  }

  return result;
}

void ShiftedMetric::shiftZ(const BoutReal* in, const dcomplex* phs, BoutReal* out) const {
  Array<dcomplex> cmplx(nmodes);

  // Take forward FFT
  rfft(in, mesh.zend + 1 - mesh.zstart, &cmplx[0]);

  // Following is an algorithm approach to write a = a*b where a and b are
  // vectors of dcomplex.
  //  std::transform(cmplxOneOff.begin(),cmplxOneOff.end(), ptr.begin(),
  //		 cmplxOneOff.begin(), std::multiplies<dcomplex>());

  for (int jz = 1; jz < nmodes; jz++) {
    cmplx[jz] *= phs[jz];
  }

  irfft(&cmplx[0], mesh.zend + 1 - mesh.zstart, out); // Reverse FFT
}

void ShiftedMetric::calcParallelSlices(Field3D& f) {
  if (f.getDirectionY() == YDirectionType::Aligned) {
    // Cannot calculate parallel slices for field-aligned fields, so return without
    // setting yup or ydown
    return;
  }

  auto results = shiftZ(f, parallel_slice_phases);

  ASSERT3(results.size() == parallel_slice_phases.size());

  f.splitParallelSlices();

  for (std::size_t i = 0; i < results.size(); ++i) {
    f.ynext(parallel_slice_phases[i].y_offset) = std::move(results[i]);
  }
}

std::vector<Field3D>
ShiftedMetric::shiftZ(const Field3D& f,
                      const std::vector<ParallelSlicePhase>& phases) const {
  ASSERT1(f.getMesh() == &mesh);
  ASSERT1(f.getLocation() == location);
  ASSERT1(f.getDirectionY() == YDirectionType::Standard);

  const int number_z_points = (mesh.zend + 1 - mesh.zstart);
  const int nmodes = number_z_points / 2 + 1;

  // FFT in Z of input field at each (x, y) point
  Matrix<Array<dcomplex>> f_fft(mesh.LocalNx, mesh.LocalNy);
  f_fft = Array<dcomplex>(nmodes);

  BOUT_FOR(i, mesh.getRegion2D("RGN_ALL")) {
    int ix = i.x();
    int iy = i.y();
    f_fft(ix, iy).ensureUnique();
    rfft(&f(i, mesh.zstart), number_z_points, &(f_fft(ix, iy)[mesh.zstart]));
  }

  std::vector<Field3D> results{};

  for (auto& phase : phases) {
    // In C++17 std::vector::emplace_back returns a reference, which
    // would be very useful here!
    results.emplace_back(&mesh);
    auto& current_result = results.back();
    current_result.allocate();
    current_result.setLocation(f.getLocation());

    BOUT_FOR(i, mesh.getRegion2D("RGN_NOY")) {
      // Deep copy the FFT'd field
      int ix = i.x();
      int iy = i.y();

      Array<dcomplex> shifted_temp(f_fft(ix, iy + phase.y_offset));
      shifted_temp.ensureUnique();

      for (int jz = 1; jz < nmodes; ++jz) {
        shifted_temp[jz] *= phase.phase_shift(ix, iy, jz);
      }

      irfft(&(shifted_temp[mesh.zstart]), number_z_points,
            &current_result(i.yp(phase.y_offset), mesh.zstart));
    }
  }

  return results;
}

// Old approach retained so we can still specify a general zShift
const Field3D ShiftedMetric::shiftZ(const Field3D& f, const Field2D& zangle,
                                    const REGION region) const {
  ASSERT1(&mesh == f.getMesh());
  ASSERT1(f.getLocation() == zangle.getLocation());
  if (mesh.zend == mesh.zstart)
    return f; // Shifting makes no difference

  Field3D result{emptyFrom(f)};

  // We only use methods in ShiftedMetric to get fields for parallel operations
  // like interp_to or DDY.
  // Therefore we don't need x-guard cells, so do not set them.
  // (Note valgrind complains about corner guard cells if we try to loop over
  // the whole grid, because zShift is not initialized in the corner guard
  // cells.)
  BOUT_FOR(i, mesh.getRegion2D(toString(region))) {
    shiftZ(&f(i, mesh.zstart), mesh.zend + 1 - mesh.zstart, zangle[i],
           &result(i, mesh.zstart));
  }

  return result;
}

void ShiftedMetric::shiftZ(const BoutReal* in, int len, BoutReal zangle,
                           BoutReal* out) const {
  int nmodes = len / 2 + 1;

  // Complex array used for FFTs
  Array<dcomplex> cmplxLoc(nmodes);

  // Take forward FFT
  rfft(in, len, &cmplxLoc[0]);

  // Apply phase shift
  for (int jz = 1; jz < nmodes; jz++) {
    BoutReal kwave = jz * 2.0 * PI / zlength; // wave number is 1/[rad]
    cmplxLoc[jz] *= dcomplex(cos(kwave * zangle), -sin(kwave * zangle));
  }

  irfft(&cmplxLoc[0], len, out); // Reverse FFT
}
