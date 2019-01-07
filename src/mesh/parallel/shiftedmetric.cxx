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

ShiftedMetric::ShiftedMetric(Mesh &mesh_in)
    : ParallelTransform(mesh_in), zShift(&thismesh) {
  // Read the zShift angle from the mesh
  
  if(thismesh.get(zShift, "zShift")) {
    // No zShift variable. Try qinty in BOUT grid files
    thismesh.get(zShift, "qinty");
  }

  // TwistShift needs to be set for derivatives to be correct at the jump where
  // poloidal angle theta goes 2pi->0
  bool twistshift = Options::root()["TwistShift"].withDefault(false);
  bool shift_without_twist = Options::root()["ShiftWithoutTwist"].withDefault(false);
  if (!twistshift and !shift_without_twist) {
    throw BoutException("ShiftedMetric usually requires the option TwistShift=true\n"
        "    Set ShiftWithoutTwist=true to use ShiftedMetric without TwistShift");
  }

  //If we wanted to be efficient we could move the following cached phase setup
  //into the relevant shifting routines (with static bool first protection)
  //so that we only calculate the phase if we actually call a relevant shift 
  //routine -- however as we're only going to do this initialisation once I 
  //think it's cleaner to put it in the constructor here.

  //As we're attached to a mesh we can expect the z direction to
  //not change once we've been created so precalculate the complex
  //phases used in transformations
  nmodes = thismesh.LocalNz / 2 + 1;
  BoutReal zlength = thismesh.getCoordinates()->zlength();

  // Allocate storage for our 3d phase information.
  fromAlignedPhs = Tensor<dcomplex>(thismesh.LocalNx, thismesh.LocalNy, nmodes);
  toAlignedPhs = Tensor<dcomplex>(thismesh.LocalNx, thismesh.LocalNy, nmodes);

  yupPhs = Tensor<dcomplex>(thismesh.LocalNx, thismesh.LocalNy, nmodes);
  ydownPhs = Tensor<dcomplex>(thismesh.LocalNx, thismesh.LocalNy, nmodes);

  //To/From field aligned phases
  BOUT_FOR(i, thismesh.getRegion2D("RGN_ALL")) {
    for(int jz=0;jz<nmodes;jz++) {
      BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
      fromAlignedPhs(i.x(), i.y(), jz) =
          dcomplex(cos(kwave * zShift[i]), -sin(kwave * zShift[i]));
      toAlignedPhs(i.x(), i.y(), jz) =
          dcomplex(cos(kwave * zShift[i]), sin(kwave * zShift[i]));
    }
  }

  //Yup/Ydown phases -- note we don't shift in the boundaries/guards
  BOUT_FOR(i, thismesh.getRegion2D("RGN_NOY")) {
    BoutReal yupShift = zShift[i] - zShift[i.yp()];
    BoutReal ydownShift = zShift[i] - zShift[i.ym()];

    for(int jz=0;jz<nmodes;jz++) {
      BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]

      yupPhs(i.x(), i.y(), jz) = dcomplex(cos(kwave * yupShift), -sin(kwave * yupShift));
      ydownPhs(i.x(), i.y(), jz) =
          dcomplex(cos(kwave * ydownShift), -sin(kwave * ydownShift));
    }
  }

}

/*!
 * Calculate the Y up and down fields
 */
void ShiftedMetric::calcYUpDown(Field3D &f) {
  ASSERT1(&thismesh == f.getMesh());

  f.splitYupYdown();
  
  Field3D& yup = f.yup();
  yup.allocate();

  BOUT_FOR(i, thismesh.getRegion2D("RGN_NOY")) {
    shiftZ(&f(i.yp(), 0), &yupPhs(i.x(), i.y(), 0), &yup(i.yp(), 0));
  }

  Field3D& ydown = f.ydown();
  ydown.allocate();

  BOUT_FOR(i, thismesh.getRegion2D("RGN_NOY")) {
    shiftZ(&f(i.ym(), 0), &ydownPhs(i.x(), i.y(), 0), &ydown(i.ym(), 0));
  }
}
  
/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftedMetric::toFieldAligned(const Field3D &f, const REGION region) {
  ASSERT2(f.getCoordinateSystem() == CoordinateSystem::Orthogonal);
  Field3D result = shiftZ(f, toAlignedPhs, region);
  result.setCoordinateSystem(CoordinateSystem::FieldAligned);
  return result;
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftedMetric::fromFieldAligned(const Field3D &f, const REGION region) {
  ASSERT2(f.getCoordinateSystem() == CoordinateSystem::FieldAligned);
  Field3D result = shiftZ(f, fromAlignedPhs, region);
  result.setCoordinateSystem(f.getMesh()->getCoordinateSystem());
  return result;
}

const Field3D ShiftedMetric::shiftZ(const Field3D& f, const Tensor<dcomplex>& phs,
                                    const REGION region) {
  ASSERT1(&thismesh == f.getMesh());
  if(thismesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&thismesh);
  result.allocate();

  BOUT_FOR(i, thismesh.getRegion2D(REGION_STRING(region))) {
    shiftZ(&f(i, 0), &phs(i.x(), i.y(), 0), &result(i, 0));
  }
  
  return result;

}

void ShiftedMetric::shiftZ(const BoutReal* in, const dcomplex* phs, BoutReal* out) {
  Array<dcomplex> cmplx(nmodes);

  // Take forward FFT
  rfft(in, thismesh.LocalNz, &cmplx[0]);

  //Following is an algorithm approach to write a = a*b where a and b are
  //vectors of dcomplex.
  //  std::transform(cmplxOneOff.begin(),cmplxOneOff.end(), ptr.begin(), 
  //		 cmplxOneOff.begin(), std::multiplies<dcomplex>());

  for(int jz=1;jz<nmodes;jz++) {
    cmplx[jz] *= phs[jz];
  }

  irfft(&cmplx[0], thismesh.LocalNz, out); // Reverse FFT
}

//Old approach retained so we can still specify a general zShift
const Field3D ShiftedMetric::shiftZ(const Field3D &f, const Field2D &zangle, const REGION region) {
  ASSERT1(&thismesh == f.getMesh());
  ASSERT1(f.getLocation() == zangle.getLocation());
  if(thismesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&thismesh);
  result.allocate();

  // We only use methods in ShiftedMetric to get fields for parallel operations
  // like interp_to or DDY.
  // Therefore we don't need x-guard cells, so do not set them.
  // (Note valgrind complains about corner guard cells if we try to loop over
  // the whole grid, because zShift is not initialized in the corner guard
  // cells.)
  BOUT_FOR(i, thismesh.getRegion2D(REGION_STRING(region))) {
    shiftZ(&f(i, 0), thismesh.LocalNz, zangle[i], &result(i, 0));
  }
  
  return result;
}

void ShiftedMetric::shiftZ(const BoutReal *in, int len, BoutReal zangle,  BoutReal *out) {
  int nmodes = len/2 + 1;

  // Complex array used for FFTs
  Array<dcomplex> cmplxLoc(nmodes);

  // Take forward FFT
  rfft(in, len, &cmplxLoc[0]);
  
  // Apply phase shift
  BoutReal zlength = thismesh.getCoordinates()->zlength();
  for(int jz=1;jz<nmodes;jz++) {
    BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
    cmplxLoc[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(&cmplxLoc[0], len, out); // Reverse FFT
}
