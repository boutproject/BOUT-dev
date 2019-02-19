/*
 * Implements the shifted metric method for parallel derivatives
 * 
 * By default fields are stored so that X-Z are orthogonal,
 * and so not aligned in Y.
 *
 */

#include <bout/paralleltransform.hxx>
#include <bout/mesh.hxx>
#include <fft.hxx>
#include <bout/constants.hxx>

#include <cmath>

#include <output.hxx>

ShiftedMetric::ShiftedMetric(Mesh &m) : ParallelTransform(m), zShift(&m) {
  // check the coordinate system used for the grid data source
  checkInputGrid();

  // Read the zShift angle from the mesh
  if (mesh.get(zShift, "zShift")) {
    // No zShift variable. Try qinty in BOUT grid files
    mesh.get(zShift, "qinty");
  }

  // TwistShift needs to be set for derivatives to be correct at the jump where
  // poloidal angle theta goes 2pi->0
  bool twistshift = Options::root()["TwistShift"].withDefault(false);
  bool shift_without_twist = Options::root()["ShiftWithoutTwist"].withDefault(false);
  if (!twistshift and !shift_without_twist) {
    throw BoutException(
        "ShiftedMetric usually requires the option TwistShift=true\n"
        "    Set ShiftWithoutTwist=true to use ShiftedMetric without TwistShift");
  }

  cachePhases();
}

ShiftedMetric::ShiftedMetric(Mesh &m, Field2D zShift_) : ParallelTransform(m), zShift(std::move(zShift_)) {
  // check the coordinate system used for the grid data source
  checkInputGrid();

  cachePhases();
}

void ShiftedMetric::checkInputGrid() {
  std::string coordinates_type = "";
  if (!mesh.get(coordinates_type, "coordinates_type")) {
    if (coordinates_type != "orthogonal") {
      throw BoutException("Incorrect coordinate system type "+coordinates_type+" used "
          "to generate metric components for ShiftedMetric. Should be 'orthogonal.");
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
  int nmodes = mesh.LocalNz / 2 + 1;
  BoutReal zlength = mesh.getCoordinates()->zlength();

  // Allocate storage for our 3d vector structures.
  // This could be made more succinct but this approach is fairly
  // verbose --> transparent
  fromAlignedPhs.resize(mesh.LocalNx);
  toAlignedPhs.resize(mesh.LocalNx);

  for (int jx = 0; jx < mesh.LocalNx; jx++) {
    fromAlignedPhs[jx].resize(mesh.LocalNy);
    toAlignedPhs[jx].resize(mesh.LocalNy);

    for (int jy = 0; jy < mesh.LocalNy; jy++) {
      fromAlignedPhs[jx][jy].resize(nmodes);
      toAlignedPhs[jx][jy].resize(nmodes);
    }
  }

  // To/From field aligned phases
  for (int jx = 0; jx < mesh.LocalNx; jx++) {
    for (int jy = 0; jy < mesh.LocalNy; jy++) {
      for (int jz = 0; jz < nmodes; jz++) {
        BoutReal kwave = jz * 2.0 * PI / zlength; // wave number is 1/[rad]
        fromAlignedPhs[jx][jy][jz] =
            dcomplex(cos(kwave * zShift(jx, jy)), -sin(kwave * zShift(jx, jy)));
        toAlignedPhs[jx][jy][jz] =
            dcomplex(cos(kwave * zShift(jx, jy)), sin(kwave * zShift(jx, jy)));
      }
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
        arr3Dvec(mesh.LocalNx,
                 std::vector<Array<dcomplex>>(mesh.LocalNy, Array<dcomplex>(nmodes)));
    parallel_slice_phases[i].y_offset = i + 1;

    // Backwards parallel slices
    parallel_slice_phases[mesh.ystart + i].phase_shift =
        arr3Dvec(mesh.LocalNx,
                 std::vector<Array<dcomplex>>(mesh.LocalNy, Array<dcomplex>(nmodes)));
    parallel_slice_phases[mesh.ystart + i].y_offset = -(i + 1);
  }

  // Parallel slice phases -- note we don't shift in the boundaries/guards
  for (auto& slice : parallel_slice_phases) {
    for (int jx = 0; jx < mesh.LocalNx; jx++) {
      for (int jy = mesh.ystart; jy <= mesh.yend; jy++) {

        slice.phase_shift[jx][jy].ensureUnique();
        BoutReal slice_shift = zShift(jx, jy) - zShift(jx, jy + slice.y_offset);

        for (int jz = 0; jz < nmodes; jz++) {
          // wave number is 1/[rad]
          BoutReal kwave = jz * 2.0 * PI / zlength;

          slice.phase_shift[jx][jy][jz] =
              dcomplex(cos(kwave * slice_shift), -sin(kwave * slice_shift));
        }
      }
    }
  }
}

/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftedMetric::toFieldAligned(const Field3D &f) {
  return shiftZ(f, toAlignedPhs);
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftedMetric::fromFieldAligned(const Field3D &f) {
  return shiftZ(f, fromAlignedPhs);
}

const Field3D ShiftedMetric::shiftZ(const Field3D &f, const arr3Dvec &phs) const {
  ASSERT1(&mesh == f.getMesh());
  if(mesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&mesh);
  result.allocate();
  result.setLocation(f.getLocation());

  for(int jx=0;jx<mesh.LocalNx;jx++) {
    for(int jy=0;jy<mesh.LocalNy;jy++) {
      shiftZ(f(jx,jy), phs[jx][jy], result(jx,jy));
    }
  }
  
  return result;

}

void ShiftedMetric::shiftZ(const BoutReal* in, const Array<dcomplex>& phs,
                           BoutReal* out) const {

  int nmodes = mesh.LocalNz / 2 + 1;
  Array<dcomplex> cmplx(nmodes);

  // Take forward FFT
  rfft(in, mesh.LocalNz, &cmplx[0]);

  // Following is an algorithm approach to write a = a*b where a and b are
  // vectors of dcomplex.
  //  std::transform(cmplxOneOff.begin(),cmplxOneOff.end(), ptr.begin(),
  //		 cmplxOneOff.begin(), std::multiplies<dcomplex>());

  for (int jz = 1; jz < nmodes; jz++) {
    cmplx[jz] *= phs[jz];
  }

  irfft(&cmplx[0], mesh.LocalNz, out); // Reverse FFT
}


void ShiftedMetric::calcYUpDown(Field3D& f) {

  auto results = shiftZ(f, parallel_slice_phases);

  ASSERT3(results.size() == parallel_slice_phases.size());

  f.splitYupYdown();

  for (std::size_t i = 0; i < results.size(); ++i) {
    f.ynext(parallel_slice_phases[i].y_offset) = std::move(results[i]);
  }
}

std::vector<Field3D>
ShiftedMetric::shiftZ(const Field3D& f,
                      const std::vector<ParallelSlicePhase>& phases) const {

  const int nmodes = mesh.LocalNz / 2 + 1;

  // FFT in Z of input field at each (x, y) point
  arr3Dvec f_fft(mesh.LocalNx,
                 std::vector<Array<dcomplex>>(mesh.LocalNy, Array<dcomplex>(nmodes)));

  for (int jx = 0; jx < mesh.LocalNx; jx++) {
    for (int jy = 0; jy < mesh.LocalNy; jy++) {
      f_fft[jx][jy].ensureUnique();
      rfft(f(jx, jy), mesh.LocalNz, f_fft[jx][jy].begin());
    }
  }

  std::vector<Field3D> results{};

  for (auto& phase : phases) {
    // In C++17 std::vector::emplace_back returns a reference, which
    // would be very useful here!
    results.emplace_back(&mesh);
    auto& current_result = results.back();
    current_result.allocate();
    current_result.setLocation(f.getLocation());

    for (int jx = 0; jx < mesh.LocalNx; jx++) {
      for (int jy = mesh.ystart; jy <= mesh.yend; jy++) {

        // Deep copy the FFT'd field
        Array<dcomplex> shifted_temp(f_fft[jx][jy + phase.y_offset]);
        shifted_temp.ensureUnique();

        for (int jz = 1; jz < nmodes; ++jz) {
          shifted_temp[jz] *= phase.phase_shift[jx][jy][jz];
        }

        irfft(shifted_temp.begin(), mesh.LocalNz,
              current_result(jx, jy + phase.y_offset));
      }
    }
  }

  return results;
}

//Old approach retained so we can still specify a general zShift
const Field3D ShiftedMetric::shiftZ(const Field3D &f, const Field2D &zangle) const {
  ASSERT1(&mesh == f.getMesh());
  if(mesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&mesh);
  result.allocate();
  result.setLocation(f.getLocation());

  for(int jx=0;jx<mesh.LocalNx;jx++) {
    for(int jy=0;jy<mesh.LocalNy;jy++) {
      shiftZ(f(jx,jy), mesh.LocalNz, zangle(jx,jy), result(jx,jy));
    }
  }
  
  return result;
}

void ShiftedMetric::shiftZ(const BoutReal* in, int len, BoutReal zangle, BoutReal* out) const {
  int nmodes = len / 2 + 1;

  // Complex array used for FFTs
  Array<dcomplex> cmplxLoc(nmodes);

  // Take forward FFT
  rfft(in, len, &cmplxLoc[0]);

  // Apply phase shift
  BoutReal zlength = mesh.getCoordinates()->zlength();
  for (int jz = 1; jz < nmodes; jz++) {
    BoutReal kwave = jz * 2.0 * PI / zlength; // wave number is 1/[rad]
    cmplxLoc[jz] *= dcomplex(cos(kwave * zangle), -sin(kwave * zangle));
  }

  irfft(&cmplxLoc[0], len, out); // Reverse FFT
}
