/*
 * Shared methods for all ParallelTransform classes
 *
 */

#include <bout/constants.hxx>
#include "bout/mesh.hxx"
#include "bout/paralleltransform.hxx"
#include <fft.hxx>

/*
 * Use FFT to shift by an angle in the Z direction
 */
void ParallelTransform::shiftZ(Field3D &f, int jx, int jy, double zangle) {
  TRACE("shiftZ");
  checkData(f);
  f.allocate(); // Ensure that f is unique
  Mesh *localmesh = f.getMesh();

  int ncz = localmesh->LocalNz;
  if(ncz == 1)
    return; // Shifting doesn't do anything

  Array<dcomplex> v(ncz/2 + 1);

  rfft(&(f(jx,jy,0)), ncz, v.begin()); // Forward FFT

  BoutReal zlength = f.getCoordinates()->zlength();

  // Apply phase shift
  for(int jz=1;jz<=ncz/2;jz++) {
    BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
    v[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(v.begin(), ncz, &(f(jx,jy,0))); // Reverse FFT
}

void ParallelTransform::shiftZ(Field3D &var, double zangle, const std::string& rgn) {
  const auto region_str = toString(rgn);

  // Only allow a whitelist of regions for now
  ASSERT2(region_str == "RGN_ALL" || region_str == "RGN_NOBNDRY" ||
          region_str == "RGN_NOX" || region_str == "RGN_NOY");

  const Region<Ind2D> &region = var.getRegion2D(region_str);

  // Could be OpenMP if shiftZ(Field3D, int, int, double) didn't throw
  BOUT_FOR_SERIAL(i, region) {
    shiftZ(var, i.x(), i.y(), zangle);
  }
}
