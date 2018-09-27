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

class Datafile;

ShiftedMetric::ShiftedMetric(Mesh &m) : mesh(m), zShift(&m) {
  // Read the zShift angle from the mesh
  
  if(mesh.get(zShift, "zShift")) {
    // No zShift variable. Try qinty in BOUT grid files
    mesh.get(zShift, "qinty");
  }

  //If we wanted to be efficient we could move the following cached phase setup
  //into the relevant shifting routines (with static bool first protection)
  //so that we only calculate the phase if we actually call a relevant shift 
  //routine -- however as we're only going to do this initialisation once I 
  //think it's cleaner to put it in the constructor here.

  //As we're attached to a mesh we can expect the z direction to
  //not change once we've been created so precalculate the complex
  //phases used in transformations
  int nmodes = mesh.LocalNz/2 + 1;
  BoutReal zlength = mesh.coordinates()->zlength();

  //Allocate storage for complex intermediate
  cmplx.resize(nmodes);
  std::fill(cmplx.begin(), cmplx.end(), 0.0);

  //Allocate storage for our 3d vector structures.
  //This could be made more succinct but this approach is fairly
  //verbose --> transparent
  fromAlignedPhs.resize(mesh.LocalNx);
  toAlignedPhs.resize(mesh.LocalNx);
  
  yupPhs.resize(mesh.LocalNx);
  ydownPhs.resize(mesh.LocalNx);

  for(int jx=0;jx<mesh.LocalNx;jx++){
    fromAlignedPhs[jx].resize(mesh.LocalNy);
    toAlignedPhs[jx].resize(mesh.LocalNy);

    yupPhs[jx].resize(mesh.LocalNy);
    ydownPhs[jx].resize(mesh.LocalNy);
    for(int jy=0;jy<mesh.LocalNy;jy++){
      fromAlignedPhs[jx][jy].resize(nmodes);
      toAlignedPhs[jx][jy].resize(nmodes);
      
      yupPhs[jx][jy].resize(nmodes);
      ydownPhs[jx][jy].resize(nmodes);
    }
  }
	
  //To/From field aligned phases
  for(int jx=0;jx<mesh.LocalNx;jx++){
    for(int jy=0;jy<mesh.LocalNy;jy++){
      for(int jz=0;jz<nmodes;jz++) {
  	BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
  	fromAlignedPhs[jx][jy][jz] = dcomplex(cos(kwave*zShift(jx,jy)) , -sin(kwave*zShift(jx,jy)));
  	toAlignedPhs[jx][jy][jz] =   dcomplex(cos(kwave*zShift(jx,jy)) ,  sin(kwave*zShift(jx,jy)));
      }
    }
  }

  //Yup/Ydown phases -- note we don't shift in the boundaries/guards
  for(int jx=0;jx<mesh.LocalNx;jx++){
    for(int jy=mesh.ystart;jy<=mesh.yend;jy++){
      BoutReal yupShift = zShift(jx,jy) - zShift(jx,jy+1);
      BoutReal ydownShift = zShift(jx,jy) - zShift(jx,jy-1);
      
      for(int jz=0;jz<nmodes;jz++) {
  	BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]

  	yupPhs[jx][jy][jz] = dcomplex(cos(kwave*yupShift) , -sin(kwave*yupShift));
  	ydownPhs[jx][jy][jz] = dcomplex(cos(kwave*ydownShift) , -sin(kwave*ydownShift));
      }
    }
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
      shiftZ(&(f(jx,jy+1,0)), yupPhs[jx][jy], &(yup(jx,jy+1,0)));
    }
  }

  Field3D& ydown = f.ydown();
  ydown.allocate();

  for(int jx=0;jx<mesh.LocalNx;jx++) {
    for(int jy=mesh.ystart;jy<=mesh.yend;jy++) {
      shiftZ(&(f(jx,jy-1,0)), ydownPhs[jx][jy], &(ydown(jx,jy-1,0)));
    }
  }
}
  
/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftedMetric::toFieldAligned(const Field3D &f, const REGION region) {
  return shiftZ(f, toAlignedPhs, region);
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftedMetric::fromFieldAligned(const Field3D &f, const REGION region) {
  return shiftZ(f, fromAlignedPhs, region);
}

const Field3D ShiftedMetric::shiftZ(const Field3D &f, const arr3Dvec &phs, const REGION region) {
  ASSERT1(&mesh == f.getMesh());
  ASSERT1(region == RGN_NOX || region == RGN_NOBNDRY); // Never calculate x-guard cells here
  if(mesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&mesh);
  result.allocate();

  for (const auto &i : f.region2D(region)) {
    shiftZ(f(i.x,i.y), phs[i.x][i.y], result(i.x,i.y));
  }
  
  return result;

}

void ShiftedMetric::shiftZ(const BoutReal *in, const std::vector<dcomplex> &phs, BoutReal *out) {
  // Take forward FFT
  rfft(in, mesh.LocalNz, &cmplx[0]);

  //Following is an algorithm approach to write a = a*b where a and b are
  //vectors of dcomplex.
  //  std::transform(cmplxOneOff.begin(),cmplxOneOff.end(), ptr.begin(), 
  //		 cmplxOneOff.begin(), std::multiplies<dcomplex>());

  const int nmodes = cmplx.size();
  for(int jz=1;jz<nmodes;jz++) {
    cmplx[jz] *= phs[jz];
  }

  irfft(&cmplx[0], mesh.LocalNz, out); // Reverse FFT
}

//Old approach retained so we can still specify a general zShift
const Field3D ShiftedMetric::shiftZ(const Field3D &f, const Field2D &zangle, const REGION region) {
  ASSERT1(&mesh == f.getMesh());
  ASSERT1(region == RGN_NOX || region == RGN_NOBNDRY); // Never calculate x-guard cells here
  ASSERT1(f.getLocation() == zangle.getLocation());
  if(mesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&mesh);
  result.allocate();

  // We only use methods in ShiftedMetric to get fields for parallel operations
  // like interp_to or DDY.
  // Therefore we don't need x-guard cells, so do not set them.
  // (Note valgrind complains about corner guard cells if we try to loop over
  // the whole grid, because zShift is not initialized in the corner guard
  // cells.)
  for(const auto &i : f.region2D(region)) {
    shiftZ(f(i.x, i.y), mesh.LocalNz, zangle(i.x,i.y), result(i.x, i.y));
  }
  
  return result;
}

void ShiftedMetric::shiftZ(const BoutReal *in, int len, BoutReal zangle,  BoutReal *out) {
  int nmodes = len/2 + 1;

  // Complex array used for FFTs
  cmplxLoc.resize(nmodes);
  
  // Take forward FFT
  rfft(in, len, &cmplxLoc[0]);
  
  // Apply phase shift
  BoutReal zlength = mesh.coordinates()->zlength();
  for(int jz=1;jz<nmodes;jz++) {
    BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
    cmplxLoc[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(&cmplxLoc[0], len, out); // Reverse FFT
}

void ShiftedMetric::outputVars(Datafile &file) {
  file.add(zShift, "zShift", 0);
}
