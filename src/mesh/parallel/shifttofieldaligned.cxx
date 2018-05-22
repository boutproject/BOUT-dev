/*
 * Implements the shifted metric method for parallel derivatives
 *
 * By default fields are stored so that X-Z are orthogonal,
 * and so not aligned in Y.
 *
 */

#include <bout/paralleltransform.hxx>
#include <bout/mesh.hxx>
#include <interpolation.hxx>
#include <fft.hxx>
#include <bout/constants.hxx>

#include <cmath>
#include <msg_stack.hxx>

#include <output.hxx>

ShiftToFieldAligned::ShiftToFieldAligned(Mesh &m) :
  mesh(m), zShift_CENTRE(&m), zShift_XLOW(&m), zShift_YLOW(&m),
  has_toAligned_CENTRE(false), has_toAligned_XLOW(false), has_toAligned_YLOW(false),
  has_fromAligned_CENTRE(false), has_fromAligned_XLOW(false), has_fromAligned_YLOW(false) {

  // Read the zShift angle from the mesh
  if(mesh.get(zShift_CENTRE, "zShift")) {
    // No zShift variable. Try qinty in BOUT grid files
    mesh.get(zShift_CENTRE, "qinty");
  }

  if (mesh.xstart >=2) {
    // Can interpolate in x-direction
    // Calculate staggered field for zShift and apply boundary conditions
    zShift_XLOW = interp_to(zShift_CENTRE, CELL_XLOW, RGN_ALL);
    zShift_XLOW.applyBoundary("neumann"); // Set boundary guard cells to closest grid cell value
  }
  if (mesh.ystart >=2) {
    // Can interpolate in y-direction
    // Calculate staggered field for zShift and apply boundary conditions
    zShift_YLOW = interp_to(zShift_CENTRE, CELL_YLOW, RGN_ALL);
    zShift_YLOW.applyBoundary("neumann"); // Set boundary guard cells to closest grid cell value
  }

  int nmodes = mesh.LocalNz/2 + 1;
  //Allocate storage for complex intermediate
  cmplx = Array<dcomplex>(nmodes);
  std::fill(cmplx.begin(), cmplx.end(), 0.0);
}

//As we're attached to a mesh we can expect the z direction to not change
//once we've been created so cache the complex phases used in transformations
//the first time they are needed
Matrix< Array<dcomplex> > ShiftToFieldAligned::getFromAlignedPhs(CELL_LOC location) {
  switch (location) {
  case CELL_CENTRE: {
    if (!has_toAligned_CENTRE) {
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.getCoordinates()->zlength();

      has_toAligned_CENTRE = true;
      fromAlignedPhs_CENTRE = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : fromAlignedPhs_CENTRE) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=mesh.xstart; jx<=mesh.xend; jx++) {
        for(int jy=0;jy<mesh.LocalNy;jy++) {
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            fromAlignedPhs_CENTRE(jx, jy)[jz] = dcomplex(cos(kwave*zShift_CENTRE(jx,jy)) , -sin(kwave*zShift_CENTRE(jx,jy)));
          }
        }
      }
    }
    return fromAlignedPhs_CENTRE;
    break;
  }
  case CELL_XLOW: {
    if (!has_toAligned_XLOW) {
      ASSERT1(mesh.xstart>=2); //otherwise we cannot interpolate in the x-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.getCoordinates()->zlength();

      has_toAligned_XLOW = true;
      fromAlignedPhs_XLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : fromAlignedPhs_XLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=mesh.xstart; jx<=mesh.xend; jx++) {
        for(int jy=0;jy<mesh.LocalNy;jy++) {
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            fromAlignedPhs_XLOW(jx, jy)[jz] = dcomplex(cos(kwave*zShift_XLOW(jx,jy)), -sin(kwave*zShift_XLOW(jx,jy)));
          }
        }
      }
    }
    return fromAlignedPhs_XLOW;
    break;
  }
  case CELL_YLOW: {
    if (!has_toAligned_YLOW) {
      ASSERT1(mesh.ystart>=2); //otherwise we cannot interpolate in the y-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.getCoordinates()->zlength();

      has_toAligned_YLOW = true;
      fromAlignedPhs_YLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : fromAlignedPhs_YLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=mesh.xstart; jx<=mesh.xend; jx++) {
        for(int jy=0;jy<mesh.LocalNy;jy++) {
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            fromAlignedPhs_YLOW(jx, jy)[jz] = dcomplex(cos(kwave*zShift_YLOW(jx,jy)), -sin(kwave*zShift_YLOW(jx,jy)));
          }
        }
      }
    }
    return fromAlignedPhs_YLOW;
    break;
  }
  case CELL_ZLOW: {
    // shifts don't depend on z, so are the same for CELL_ZLOW as for CELL_CENTRE
    return getFromAlignedPhs(CELL_CENTRE);
    break;
  }
  default: {
    // This should never happen
    throw BoutException("Unsupported stagger of phase shifts\n"
                        " - don't know how to interpolate to %s",strLocation(location));
    break;
  }
  };
}

Matrix< Array<dcomplex> > ShiftToFieldAligned::getToAlignedPhs(CELL_LOC location) {
  switch (location) {
  case CELL_CENTRE: {
    if (!has_fromAligned_CENTRE) {
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.getCoordinates()->zlength();

      has_fromAligned_CENTRE = true;
      toAlignedPhs_CENTRE = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : toAlignedPhs_CENTRE) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=mesh.xstart; jx<=mesh.xend; jx++) {
        for(int jy=0;jy<mesh.LocalNy;jy++) {
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            toAlignedPhs_CENTRE(jx, jy)[jz] = dcomplex(cos(kwave*zShift_CENTRE(jx,jy)), sin(kwave*zShift_CENTRE(jx,jy)));
          }
        }
      }
    }
    return toAlignedPhs_CENTRE;
    break;
  }
  case CELL_XLOW: {
    if (!has_fromAligned_XLOW) {
      ASSERT1(mesh.xstart>=2); //otherwise we cannot interpolate in the x-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.getCoordinates()->zlength();

      has_fromAligned_XLOW = true;
      toAlignedPhs_XLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : toAlignedPhs_XLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=mesh.xstart; jx<=mesh.xend; jx++) {
        for(int jy=0;jy<mesh.LocalNy;jy++) {
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            toAlignedPhs_XLOW(jx, jy)[jz] = dcomplex(cos(kwave*zShift_XLOW(jx,jy)), sin(kwave*zShift_XLOW(jx,jy)));
          }
        }
      }
    }
    return toAlignedPhs_XLOW;
    break;
  }
  case CELL_YLOW: {
    if (!has_fromAligned_YLOW) {
      ASSERT1(mesh.ystart>=2); //otherwise we cannot interpolate in the y-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.getCoordinates()->zlength();

      has_fromAligned_YLOW = true;
      toAlignedPhs_YLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : toAlignedPhs_YLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=mesh.xstart; jx<=mesh.xend; jx++) {
        for(int jy=0;jy<mesh.LocalNy;jy++) {
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            toAlignedPhs_YLOW(jx, jy)[jz] = dcomplex(cos(kwave*zShift_YLOW(jx,jy)), sin(kwave*zShift_YLOW(jx,jy)));
          }
        }
      }
    }
    return toAlignedPhs_YLOW;
    break;
  }
  case CELL_ZLOW: {
    // shifts don't depend on z, so are the same for CELL_ZLOW as for CELL_CENTRE
    return getToAlignedPhs(CELL_CENTRE);
    break;
  }
  default: {
    // This should never happen
    throw BoutException("Unsupported stagger of phase shifts\n"
                        " - don't know how to interpolate to %s",strLocation(location));
    break;
  }
  };
}


/*!
 * Calculate the field-aligned field
 */
void ShiftToFieldAligned::calcYUpDown(Field3D &f, REGION region) {
  ASSERT1(region == RGN_NOX || region == RGN_NOBNDRY);
  ASSERT1(&mesh == f.getMesh());

  // We only use methods in ShiftToFieldAligned to get fields for parallel operations
  // like interp_to or DDY.
  // Therefore we don't need x-guard cells, so do not set them.
  // (Note valgrind complains about corner guard cells if we try to loop over
  // the whole grid, because zShift is not initialized in the corner guard
  // cells.)

  Field3D& f_fa = f.fieldAligned();
  f_fa = shiftZ(f, getToAlignedPhs(f.getLocation()), region);
  f.setHasFieldAligned(true);
}

/*!
 * Get the shifted field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftToFieldAligned::toFieldAligned(const Field3D &f, const REGION region) {
  ASSERT1(region==RGN_NOX || region==RGN_NOBNDRY);
  if (f.hasFieldAligned()) {
    return f.fieldAligned();
  } else {
#if CHECK > 1
    static int count = 0;

    if (count<100) {
      // Get a trace of where we were called from
      std::string message = msg_stack.getDump();

      output<<"Warning:"<<endl;
      output<<"Called toFieldAligned on a field without a field_fa member (i.e. f.hasFieldAligned()==false)."<<endl;
      output<<"This may be suboptimal."<<endl;
      output<<"Was called from:"<<endl;
      output<<message<<endl;

      count++;
    } else if (count == 100) {
      output<<"More than 100 warnings from ShiftToFieldAligned::toFieldAligned."<<endl
        <<"Suppressing further output from here."<<endl;
      count++;
    }
#endif
    return shiftZ(f, getToAlignedPhs(f.getLocation()), region);
  }
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftToFieldAligned::fromFieldAligned(const Field3D &f, const REGION region) {
  return shiftZ(f, getFromAlignedPhs(f.getLocation()), region);
}

const Field3D ShiftToFieldAligned::shiftZ(const Field3D &f, const Matrix< Array<dcomplex> > &phs, const REGION region) {
  ASSERT1(&mesh == f.getMesh());
  ASSERT1(region == RGN_NOX || region == RGN_NOBNDRY); // Never calculate x-guard cells here
  if(mesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&mesh);
  result.allocate();
  result.setLocation(f.getLocation());
  //result = 0.; // Set to value to avoid uninitialized value errors from Valgrind

  invalidateGuards(result); // Won't set x-guard cells, so allow checking to throw exception if they are used.

  // We only use methods in ShiftToFieldAligned to get fields for parallel operations
  // like interp_to or DDY.
  // Therefore we don't need x-guard cells, so do not set them.
  // (Note valgrind complains about corner guard cells if we try to loop over
  // the whole grid, because zShift is not initialized in the corner guard
  // cells.)
  for (const auto &i : mesh.getRegion2D(REGION_STRING(region))) {
    shiftZ(f(i.x(), i.y()), phs(i.x(), i.y()), result(i.x(), i.y()));
  }

  return result;

}

void ShiftToFieldAligned::shiftZ(const BoutReal *in, const Array<dcomplex> &phs, BoutReal *out) {
  // Take forward FFT
  rfft(in, mesh.LocalNz, cmplx.begin());

  //Following is an algorithm approach to write a = a*b where a and b are
  //vectors of dcomplex.
  //  std::transform(cmplxOneOff.begin(),cmplxOneOff.end(), ptr.begin(),
  //                 cmplxOneOff.begin(), std::multiplies<dcomplex>());

  const int nmodes = cmplx.size();
  for(int jz=1;jz<nmodes;jz++) {
    cmplx[jz] *= phs[jz];
  }

  irfft(cmplx.begin(), mesh.LocalNz, out); // Reverse FFT
}

//Old approach retained so we can still specify a general zShift
const Field3D ShiftToFieldAligned::shiftZ(const Field3D &f, const Field2D &zangle, const REGION region) {
  ASSERT1(&mesh == f.getMesh());
  ASSERT1(region == RGN_NOX || region == RGN_NOBNDRY); // Never calculate x-guard cells here
  ASSERT1(f.getLocation() == zangle.getLocation());
  if(mesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&mesh);
  result.allocate();
  invalidateGuards(result); // Won't set x-guard cells, so allow checking to throw exception if they are used.

  // We only use methods in ShiftToFieldAligned to get fields for parallel operations
  // like interp_to or DDY.
  // Therefore we don't need x-guard cells, so do not set them.
  // (Note valgrind complains about corner guard cells if we try to loop over
  // the whole grid, because zShift is not initialized in the corner guard
  // cells.)
  for (const auto &i : mesh.getRegion2D(REGION_STRING(region))) {
    shiftZ(f(i.x(), i.y()), mesh.LocalNz, zangle(i.x(),i.y()), result(i.x(), i.y()));
  }

  return result;
}

void ShiftToFieldAligned::shiftZ(const BoutReal *in, int len, BoutReal zangle,  BoutReal *out) {
  int nmodes = len/2 + 1;

  // Complex array used for FFTs
  cmplxLoc = Array<dcomplex>(nmodes);

  // Take forward FFT
  rfft(in, len, cmplxLoc.begin());

  // Apply phase shift
  BoutReal zlength = mesh.getCoordinates()->zlength();
  for(int jz=1;jz<nmodes;jz++) {
    BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
    cmplxLoc[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(cmplxLoc.begin(), len, out); // Reverse FFT
}
