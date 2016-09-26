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

ShiftedMetric::ShiftedMetric(Mesh &m) : mesh(m) {
  // Read the zShift angle from the mesh
  
  if(mesh.get(zShift, "zShift")) {
    // No zShift variable. Try qinty in BOUT grid files
    mesh.get(zShift, "qinty");
  }
}

/*!
 * Calculate the Y up and down fields
 */
void ShiftedMetric::calcYUpDown(Field3D &f) {
  f.splitYupYdown();
  
  Field3D& yup = f.yup();
  yup.allocate();

  for(int jx=0;jx<mesh.LocalNx;jx++) {
    for(int jy=mesh.ystart;jy<=mesh.yend;jy++) {
      shiftZ(&(f(jx,jy+1,0)), mesh.LocalNz, zShift(jx,jy) - zShift(jx,jy+1), &(yup(jx,jy+1,0)));
    }
  }

  Field3D& ydown = f.ydown();
  ydown.allocate();

  for(int jx=0;jx<mesh.LocalNx;jx++) {
    for(int jy=mesh.ystart;jy<=mesh.yend;jy++) {
      shiftZ(&(f(jx,jy-1,0)), mesh.LocalNz, zShift(jx,jy) - zShift(jx,jy-1), &(ydown(jx,jy-1,0)));
    }
  }
}
  
/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftedMetric::toFieldAligned(const Field3D &f) {
  return shiftZ(f, -zShift);
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftedMetric::fromFieldAligned(const Field3D &f) {
  return shiftZ(f, zShift);
}

const Field3D ShiftedMetric::shiftZ(const Field3D f, const Field2D zangle) {
  if(mesh.LocalNz == 1)
    return f; // Shifting makes no difference
  
  Field3D result;
  result.allocate();

  for(int jx=0;jx<mesh.LocalNx;jx++) {
    for(int jy=0;jy<mesh.LocalNy;jy++) {
      shiftZ(f(jx,jy), mesh.LocalNz, zangle(jx,jy), result(jx,jy));
    }
  }
  
  return result;
}

void ShiftedMetric::shiftZ(const BoutReal *in, int len, BoutReal zangle,  BoutReal *out) {
  int nmodes = len/2 + 1;

  // Complex array used for FFTs
  cmplx.resize(nmodes);
  
  // Take forward FFT
  rfft(in, len, &cmplx[0]);
  
  // Apply phase shift
  BoutReal zlength = mesh.coordinates()->zlength();
  for(int jz=1;jz<nmodes;jz++) {
    BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
    cmplx[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(&cmplx[0], len, out); // Reverse FFT
}
