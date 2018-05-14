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

#include <output.hxx>

ShiftedMetric::ShiftedMetric(Mesh &m) : mesh(m), zShift(&m) {
  // Read the zShift angle from the mesh
  
  if(mesh.get(zShift, "zShift")) {
    // No zShift variable. Try qinty in BOUT grid files
    mesh.get(zShift, "qinty");
  }

  if (mesh.xstart >=2) {
    // Can interpolate in x-direction
    // Calculate staggered field for zShift and apply boundary conditions
    Field2D zShift_XLOW = interp_to(zShift, CELL_XLOW, RGN_ALL);
    zShift_XLOW.applyBoundary("neumann"); // Set boundary guard cells to closest grid cell value
    zShift.set(zShift_XLOW);
  }
  if (mesh.ystart >=2) {
    // Can interpolate in y-direction
    // Calculate staggered field for zShift and apply boundary conditions
    Field2D zShift_YLOW = interp_to(zShift, CELL_YLOW, RGN_ALL);
    zShift_YLOW.applyBoundary("neumann"); // Set boundary guard cells to closest grid cell value
    zShift.set(zShift_YLOW);
  }

  int nmodes = mesh.LocalNz/2 + 1;
  //Allocate storage for complex intermediate
  cmplx = Array<dcomplex>(nmodes);
  std::fill(cmplx.begin(), cmplx.end(), 0.0);
}

//As we're attached to a mesh we can expect the z direction to not change
//once we've been created so cache the complex phases used in transformations
//the first time they are needed
Matrix< Array<dcomplex> > ShiftedMetric::getFromAlignedPhs(CELL_LOC location) {
  // bools so we only calculate the cached values the first time for each location
  static bool first_CENTRE = true, first_XLOW=true, first_YLOW=true;

  switch (location) {
  case CELL_CENTRE: {
    if (first_CENTRE) {
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      first_CENTRE = false;
      fromAlignedPhs_CENTRE = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : fromAlignedPhs_CENTRE) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            fromAlignedPhs_CENTRE(jx, jy)[jz] = dcomplex(cos(kwave*zShift(jx,jy)) , -sin(kwave*zShift(jx,jy)));
          }
        }
      }
    }
    return fromAlignedPhs_CENTRE;
    break;
  }
  case CELL_XLOW: {
    if (first_XLOW) {
      ASSERT1(mesh.xstart>=2); //otherwise we cannot interpolate in the x-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at XLOW
      Field2D zShift_XLOW = zShift.get(CELL_XLOW);

      first_XLOW = false;
      fromAlignedPhs_XLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : fromAlignedPhs_XLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
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
    if (first_YLOW) {
      ASSERT1(mesh.ystart>=2); //otherwise we cannot interpolate in the y-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at YLOW
      Field2D zShift_YLOW = zShift.get(CELL_YLOW);

      first_YLOW = false;
      fromAlignedPhs_YLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : fromAlignedPhs_YLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
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

Matrix< Array<dcomplex> > ShiftedMetric::getToAlignedPhs(CELL_LOC location) {
  // bools so we only calculate the cached values the first time for each location
  static bool first_CENTRE = true, first_XLOW=true, first_YLOW=true;

  switch (location) {
  case CELL_CENTRE: {
    if (first_CENTRE) {
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      first_CENTRE = false;
      toAlignedPhs_CENTRE = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : toAlignedPhs_CENTRE) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            toAlignedPhs_CENTRE(jx, jy)[jz] = dcomplex(cos(kwave*zShift(jx,jy)), sin(kwave*zShift(jx,jy)));
          }
        }
      }
    }
    return toAlignedPhs_CENTRE;
    break;
  }
  case CELL_XLOW: {
    if (first_XLOW) {
      ASSERT1(mesh.xstart>=2); //otherwise we cannot interpolate in the x-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at XLOW
      Field2D zShift_XLOW = zShift.get(CELL_XLOW);

      first_XLOW = false;
      toAlignedPhs_XLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : toAlignedPhs_XLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
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
    if (first_YLOW) {
      ASSERT1(mesh.ystart>=2); //otherwise we cannot interpolate in the y-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at YLOW
      Field2D zShift_YLOW = zShift.get(CELL_YLOW);

      first_YLOW = false;
      toAlignedPhs_YLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : toAlignedPhs_YLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
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

Matrix< Array<dcomplex> > ShiftedMetric::getYupPhs1(CELL_LOC location) {
  // bools so we only calculate the cached values the first time for each location
  static bool first_CENTRE = true, first_XLOW=true, first_YLOW=true;

  switch (location) {
  case CELL_CENTRE: {
    if (first_CENTRE) {
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      first_CENTRE = false;
      yupPhs1_CENTRE = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : yupPhs1_CENTRE) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal yupShift1 = zShift(jx,jy) - zShift(jx,jy+1);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            yupPhs1_CENTRE(jx, jy)[jz] = dcomplex(cos(kwave*yupShift1) , -sin(kwave*yupShift1));
          }
        }
      }
    }
    return yupPhs1_CENTRE;
    break;
  }
  case CELL_XLOW: {
    if (first_XLOW) {
      ASSERT1(mesh.xstart>=2); //otherwise we cannot interpolate in the x-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at XLOW
      Field2D zShift_XLOW = zShift.get(CELL_XLOW);

      first_XLOW = false;
      yupPhs1_XLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : yupPhs1_XLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal yupShift1 = zShift_XLOW(jx,jy) - zShift_XLOW(jx,jy+1);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            yupPhs1_XLOW(jx, jy)[jz] = dcomplex(cos(kwave*yupShift1) , -sin(kwave*yupShift1));
          }
        }
      }
    }
    return yupPhs1_XLOW;
    break;
  }
  case CELL_YLOW: {
    if (first_YLOW) {
      ASSERT1(mesh.ystart>=2); //otherwise we cannot interpolate in the y-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at YLOW
      Field2D zShift_YLOW = zShift.get(CELL_YLOW);

      first_YLOW = false;
      yupPhs1_YLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : yupPhs1_YLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal yupShift1 = zShift_YLOW(jx,jy) - zShift_YLOW(jx,jy+1);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            yupPhs1_YLOW(jx, jy)[jz] = dcomplex(cos(kwave*yupShift1) , -sin(kwave*yupShift1));
          }
        }
      }
    }
    return yupPhs1_YLOW;
    break;
  }
  case CELL_ZLOW: {
    // shifts don't depend on z, so are the same for CELL_ZLOW as for CELL_CENTRE
    return getYupPhs1(CELL_CENTRE);
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

Matrix< Array<dcomplex> > ShiftedMetric::getYupPhs2(CELL_LOC location) {
  // bools so we only calculate the cached values the first time for each location
  static bool first_CENTRE = true, first_XLOW=true, first_YLOW=true;

  switch (location) {
  case CELL_CENTRE: {
    if (first_CENTRE) {
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      first_CENTRE = false;
      yupPhs2_CENTRE = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : yupPhs2_CENTRE) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal yupShift2 = zShift(jx,jy) - zShift(jx,jy+2);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            yupPhs2_CENTRE(jx, jy)[jz] = dcomplex(cos(kwave*yupShift2) , -sin(kwave*yupShift2));
          }
        }
      }
    }
    return yupPhs2_CENTRE;
    break;
  }
  case CELL_XLOW: {
    if (first_XLOW) {
      ASSERT1(mesh.xstart>=2); //otherwise we cannot interpolate in the x-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at XLOW
      Field2D zShift_XLOW = zShift.get(CELL_XLOW);

      first_XLOW = false;
      yupPhs2_XLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : yupPhs2_XLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal yupShift2 = zShift_XLOW(jx,jy) - zShift_XLOW(jx,jy+2);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            yupPhs2_XLOW(jx, jy)[jz] = dcomplex(cos(kwave*yupShift2) , -sin(kwave*yupShift2));
          }
        }
      }
    }
    return yupPhs2_XLOW;
    break;
  }
  case CELL_YLOW: {
    if (first_YLOW) {
      ASSERT1(mesh.ystart>=2); //otherwise we cannot interpolate in the y-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at YLOW
      Field2D zShift_YLOW = zShift.get(CELL_YLOW);

      first_YLOW = false;
      yupPhs2_YLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : yupPhs2_YLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal yupShift2 = zShift_YLOW(jx,jy) - zShift_YLOW(jx,jy+2);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            yupPhs2_YLOW(jx, jy)[jz] = dcomplex(cos(kwave*yupShift2) , -sin(kwave*yupShift2));
          }
        }
      }
    }
    return yupPhs2_YLOW;
    break;
  }
  case CELL_ZLOW: {
    // shifts don't depend on z, so are the same for CELL_ZLOW as for CELL_CENTRE
    return getYupPhs1(CELL_CENTRE);
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

Matrix< Array<dcomplex> > ShiftedMetric::getYdownPhs1(CELL_LOC location) {
  // bools so we only calculate the cached values the first time for each location
  static bool first_CENTRE = true, first_XLOW=true, first_YLOW=true;

  switch (location) {
  case CELL_CENTRE: {
    if (first_CENTRE) {
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      first_CENTRE = false;
      ydownPhs1_CENTRE = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : ydownPhs1_CENTRE) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal ydownShift1 = zShift(jx,jy) - zShift(jx,jy-1);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            ydownPhs1_CENTRE(jx, jy)[jz] = dcomplex(cos(kwave*ydownShift1) , -sin(kwave*ydownShift1));
          }
        }
      }
    }
    return ydownPhs1_CENTRE;
    break;
  }
  case CELL_XLOW: {
    if (first_XLOW) {
      ASSERT1(mesh.xstart>=2); //otherwise we cannot interpolate in the x-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at XLOW
      Field2D zShift_XLOW = zShift.get(CELL_XLOW);

      first_XLOW = false;
      ydownPhs1_XLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : ydownPhs1_XLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal ydownShift1 = zShift_XLOW(jx,jy) - zShift_XLOW(jx,jy-1);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            ydownPhs1_XLOW(jx, jy)[jz] = dcomplex(cos(kwave*ydownShift1) , -sin(kwave*ydownShift1));
          }
        }
      }
    }
    return ydownPhs1_XLOW;
    break;
  }
  case CELL_YLOW: {
    if (first_YLOW) {
      ASSERT1(mesh.ystart>=2); //otherwise we cannot interpolate in the y-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at YLOW
      Field2D zShift_YLOW = zShift.get(CELL_YLOW);

      first_YLOW = false;
      ydownPhs1_YLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : ydownPhs1_YLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal ydownShift1 = zShift_YLOW(jx,jy) - zShift_YLOW(jx,jy-1);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            ydownPhs1_YLOW(jx, jy)[jz] = dcomplex(cos(kwave*ydownShift1) , -sin(kwave*ydownShift1));
          }
        }
      }
    }
    return ydownPhs1_YLOW;
    break;
  }
  case CELL_ZLOW: {
    // shifts don't depend on z, so are the same for CELL_ZLOW as for CELL_CENTRE
    return getYdownPhs1(CELL_CENTRE);
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

Matrix< Array<dcomplex> > ShiftedMetric::getYdownPhs2(CELL_LOC location) {
  // bools so we only calculate the cached values the first time for each location
  static bool first_CENTRE = true, first_XLOW=true, first_YLOW=true;

  switch (location) {
  case CELL_CENTRE: {
    if (first_CENTRE) {
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      first_CENTRE = false;
      ydownPhs2_CENTRE = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : ydownPhs2_CENTRE) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal ydownShift2 = zShift(jx,jy) - zShift(jx,jy-2);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            ydownPhs2_CENTRE(jx, jy)[jz] = dcomplex(cos(kwave*ydownShift2) , -sin(kwave*ydownShift2));
          }
        }
      }
    }
    return ydownPhs2_CENTRE;
    break;
  }
  case CELL_XLOW: {
    if (first_XLOW) {
      ASSERT1(mesh.xstart>=2); //otherwise we cannot interpolate in the x-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at XLOW
      Field2D zShift_XLOW = zShift.get(CELL_XLOW);

      first_XLOW = false;
      ydownPhs2_XLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : ydownPhs2_XLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal ydownShift2 = zShift_XLOW(jx,jy) - zShift_XLOW(jx,jy-2);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            ydownPhs2_XLOW(jx, jy)[jz] = dcomplex(cos(kwave*ydownShift2) , -sin(kwave*ydownShift2));
          }
        }
      }
    }
    return ydownPhs2_XLOW;
    break;
  }
  case CELL_YLOW: {
    if (first_YLOW) {
      ASSERT1(mesh.ystart>=2); //otherwise we cannot interpolate in the y-direction
      int nmodes = mesh.LocalNz/2 + 1;
      BoutReal zlength = mesh.coordinates()->zlength();

      // get zShift at YLOW
      Field2D zShift_YLOW = zShift.get(CELL_YLOW);

      first_YLOW = false;
      ydownPhs2_YLOW = Matrix< Array<dcomplex> >(mesh.LocalNx, mesh.LocalNy);
      for (auto &element : ydownPhs2_YLOW) {
        element = Array<dcomplex>(mesh.LocalNz);
      }

      //To/From field aligned phases
      for(int jx=0;jx<mesh.LocalNx;jx++){
        for(int jy=0;jy<mesh.LocalNy;jy++){
          BoutReal ydownShift2 = zShift_YLOW(jx,jy) - zShift_YLOW(jx,jy-2);
          for(int jz=0;jz<nmodes;jz++) {
            BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
            ydownPhs2_YLOW(jx, jy)[jz] = dcomplex(cos(kwave*ydownShift2) , -sin(kwave*ydownShift2));
          }
        }
      }
    }
    return ydownPhs2_YLOW;
    break;
  }
  case CELL_ZLOW: {
    // shifts don't depend on z, so are the same for CELL_ZLOW as for CELL_CENTRE
    return getYdownPhs1(CELL_CENTRE);
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
 * Calculate the Y up and down fields
 */
void ShiftedMetric::calcYUpDown(Field3D &f) {
  f.splitYupYdown();
  CELL_LOC location = f.getLocation();
  
  Field3D& yup1 = f.yup();
  yup1.allocate();
  Matrix< Array<dcomplex> > phases = getYupPhs1(location);
  for(int jx=0;jx<mesh.LocalNx;jx++) {
    for(int jy=mesh.ystart;jy<=mesh.yend;jy++) {
      shiftZ(&(f(jx,jy+1,0)), phases(jx, jy), &(yup1(jx,jy+1,0)));
    }
  }
  if (mesh.ystart>1) {
    Field3D& yup2 = f.yup(2);
    yup2.allocate();
    phases = getYupPhs2(location);
    for(int jx=0;jx<mesh.LocalNx;jx++) {
      for(int jy=mesh.ystart;jy<=mesh.yend;jy++) {
        shiftZ(&(f(jx,jy+2,0)), phases(jx, jy), &(yup2(jx,jy+2,0)));
      }
    }
  }

  Field3D& ydown1 = f.ydown();
  ydown1.allocate();
  phases = getYdownPhs1(location);
  for(int jx=0;jx<mesh.LocalNx;jx++) {
    for(int jy=mesh.ystart;jy<=mesh.yend;jy++) {
      shiftZ(&(f(jx,jy-1,0)), phases(jx, jy), &(ydown1(jx,jy-1,0)));
    }
  }
  if (mesh.ystart > 1) {
    Field3D& ydown2 = f.ydown(2);
    ydown2.allocate();
    phases = getYdownPhs2(location);
    for(int jx=0;jx<mesh.LocalNx;jx++) {
      for(int jy=mesh.ystart;jy<=mesh.yend;jy++) {
        shiftZ(&(f(jx,jy-2,0)), phases(jx, jy), &(ydown2(jx,jy-2,0)));
      }
    }
  }
}
  
/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftedMetric::toFieldAligned(const Field3D &f) {
  return shiftZ(f, getToAlignedPhs(f.getLocation()));
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftedMetric::fromFieldAligned(const Field3D &f) {
  return shiftZ(f, getFromAlignedPhs(f.getLocation()));
}

const Field3D ShiftedMetric::shiftZ(const Field3D &f, const Matrix< Array<dcomplex> > &phs) {
  ASSERT1(&mesh == f.getMesh());
  if(mesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(f); // Initialize from f, mostly so location get set correctly. (Does not copy data because of copy-on-change).
  
  for(int jx=0;jx<mesh.LocalNx;jx++) {
    for(int jy=0;jy<mesh.LocalNy;jy++) {
      shiftZ(f(jx,jy), phs(jx, jy), result(jx,jy));
    }
  }
  
  return result;

}

void ShiftedMetric::shiftZ(const BoutReal *in, const Array<dcomplex> &phs, BoutReal *out) {
  // Take forward FFT
  rfft(in, mesh.LocalNz, cmplx.begin());

  //Following is an algorithm approach to write a = a*b where a and b are
  //vectors of dcomplex.
  //  std::transform(cmplxOneOff.begin(),cmplxOneOff.end(), ptr.begin(), 
  //		 cmplxOneOff.begin(), std::multiplies<dcomplex>());

  const int nmodes = cmplx.size();
  for(int jz=1;jz<nmodes;jz++) {
    cmplx[jz] *= phs[jz];
  }

  irfft(cmplx.begin(), mesh.LocalNz, out); // Reverse FFT
}

//Old approach retained so we can still specify a general zShift
const Field3D ShiftedMetric::shiftZ(const Field3D &f, const Field2D &zangle) {
  ASSERT1(&mesh == f.getMesh());
  ASSERT1(f.getLocation() == zangle.getLocation());
  if(mesh.LocalNz == 1)
    return f; // Shifting makes no difference

  Field3D result(&mesh);
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
  cmplxLoc = Array<dcomplex>(nmodes);
  
  // Take forward FFT
  rfft(in, len, cmplxLoc.begin());
  
  // Apply phase shift
  BoutReal zlength = mesh.coordinates()->zlength();
  for(int jz=1;jz<nmodes;jz++) {
    BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
    cmplxLoc[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(cmplxLoc.begin(), len, out); // Reverse FFT
}
