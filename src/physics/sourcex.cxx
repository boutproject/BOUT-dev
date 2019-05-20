/**************************************************************
 * radial source and mask operators
 **************************************************************/

#include <globals.hxx>
#include <cmath>

#include <bout/mesh.hxx>
#include <field2d.hxx>
#include <sourcex.hxx>
#include <msg_stack.hxx>

#include "unused.hxx"

BoutReal TanH(BoutReal a) {
  BoutReal temp = exp(a);
  return ((temp - 1.0 / temp) / (temp + 1.0 / temp));
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field2D source_tanhx(const Field2D &f, BoutReal swidth, BoutReal slength) {
  Mesh* localmesh = f.getMesh();

  Field2D result{emptyFrom(f)};
  result.allocate();

  // create a radial buffer zone to set jpar zero near radial boundary
  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    BoutReal lx = localmesh->GlobalX(i.x()) - slength;
    BoutReal dampl = TanH(lx / swidth);
    result[i] = 0.5 * (1.0 - dampl);
  }

  // Need to communicate boundaries
  localmesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field2D source_expx2(const Field2D &f, BoutReal swidth, BoutReal slength) {
  Mesh* localmesh = f.getMesh();

  Field2D result{emptyFrom(f)};

  result.allocate();

  // create a radial buffer zone to set jpar zero near radial boundary
  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    BoutReal lx = localmesh->GlobalX(i.x()) - slength;
    BoutReal dampl = exp(-lx * lx / swidth / swidth);
    result[i] = dampl;
  }

  // Need to communicate boundaries
  localmesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhx(const Field2D &UNUSED(f0), const Field3D &f, BoutReal swidth,
                         BoutReal slength, bool UNUSED(BoutRealspace)) {
  Mesh* localmesh = f.getMesh();

  Field3D result{emptyFrom(f)};
  result.allocate();

  // create a radial buffer zone to set jpar zero near radial boundary
  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    BoutReal rlx = 1. - localmesh->GlobalX(i.x()) - slength;
    BoutReal dampr = TanH(rlx / swidth);
    result[i] = 0.5 * (1.0 - dampr) * f[i];
  }

  // Need to communicate boundaries
  localmesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D mask_x(const Field3D &f, bool UNUSED(BoutRealspace)) {
  TRACE("mask_x");

  Mesh* localmesh = f.getMesh();

  Field3D result{emptyFrom(f)};
  result.allocate();

  // create a radial buffer zone to set jpar zero near radial boundary
  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    BoutReal lx = localmesh->GlobalX(i.x());
    BoutReal dampl = TanH(lx / 40.0);
    BoutReal dampr = TanH((1. - lx) / 40.0);

    result[i] = (1.0 - dampl * dampr) * f[i];
  }

  // Need to communicate boundaries
  localmesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhxl(const Field2D &UNUSED(f0), const Field3D &f, BoutReal swidth,
                          BoutReal slength, bool UNUSED(BoutRealspace)) {
  TRACE("sink_tanhx");

  Mesh* localmesh = f.getMesh();

  Field3D result{emptyFrom(f)};

  result.allocate();

  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    BoutReal lx = localmesh->GlobalX(i.x()) - slength;
    BoutReal dampl = TanH(lx / swidth);

    result[i] = 0.5 * (1.0 - dampl) * f[i];
  }

  // Need to communicate boundaries
  localmesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhxr(const Field2D &UNUSED(f0), const Field3D &f, BoutReal swidth,
                          BoutReal slength, bool UNUSED(BoutRealspace)) {
  TRACE("sink_tanhxr");

  Mesh* localmesh = f.getMesh();

  Field3D result{emptyFrom(f)};
  result.allocate();

  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    BoutReal rlx = 1. - localmesh->GlobalX(i.x()) - slength;
    BoutReal dampr = TanH(rlx / swidth);

    result[i] = 0.5 * (1.0 - dampr) * f[i];
  }

  // Need to communicate boundaries
  localmesh->communicate(result);

  return result;
}

// create radial buffer zones to damp Psi to zero near radial boundaries
const Field3D buff_x(const Field3D &f, bool UNUSED(BoutRealspace)) {
  TRACE("buff_x");

  Mesh* localmesh = f.getMesh();

  Field3D result{emptyFrom(f)};
  result.allocate();

  const BoutReal dampl = 1.e0;
  const BoutReal dampr = 1.e0;
  const BoutReal deltal = 0.05;
  const BoutReal deltar = 0.05;

  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    BoutReal lx = localmesh->GlobalX(i.x());
    BoutReal rlx = 1. - lx;

    result[i] = (dampl * exp(-(lx * lx) / (deltal * deltal)) +
                 dampr * exp(-(rlx * rlx) / (deltar * deltar))) *
                f[i];
  }

  // Need to communicate boundaries
  localmesh->communicate(result);

  return result;
}
