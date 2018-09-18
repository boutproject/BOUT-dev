/**************************************************************
 * radial source and mask operators
 **************************************************************/

#include <globals.hxx>
#include <math.h>

#include <sourcex.hxx>
#include <msg_stack.hxx>

#include "unused.hxx"

BoutReal TanH(BoutReal a) {
  BoutReal temp = exp(a);
  return ((temp - 1.0 / temp) / (temp + 1.0 / temp));
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field2D source_tanhx(const Field2D &UNUSED(f), BoutReal swidth, BoutReal slength) {
  Field2D result;
  result.allocate();

  // create a radial buffer zone to set jpar zero near radial boundary
  for (const auto &i : result) {
    BoutReal lx = mesh->GlobalX(i.x) - slength;
    BoutReal dampl = TanH(lx / swidth);
    result[i] = 0.5 * (1.0 - dampl);
  }

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field2D source_expx2(const Field2D &UNUSED(f), BoutReal swidth, BoutReal slength) {
  Field2D result;

  result.allocate();

  // create a radial buffer zone to set jpar zero near radial boundary

  for (const auto &i : result) {
    BoutReal lx = mesh->GlobalX(i.x) - slength;
    BoutReal dampl = exp(-lx * lx / swidth / swidth);
    result[i] = dampl;
  }

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhx(const Field2D &UNUSED(f0), const Field3D &f, BoutReal swidth,
                         BoutReal slength, bool UNUSED(BoutRealspace)) {
  Field3D result;
  result.allocate();

  // create a radial buffer zone to set jpar zero near radial boundary
  for (const auto &i : result) {
    BoutReal rlx = 1. - mesh->GlobalX(i.x) - slength;
    BoutReal dampr = TanH(rlx / swidth);
    result[i] = 0.5 * (1.0 - dampr) * f[i];
  }

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D mask_x(const Field3D &f, bool UNUSED(BoutRealspace)) {
  TRACE("mask_x");

  Field3D result;
  result.allocate();

  // create a radial buffer zone to set jpar zero near radial boundary
  for (const auto &i : result) {
    BoutReal lx = mesh->GlobalX(i.x);
    BoutReal dampl = TanH(lx / 40.0);
    BoutReal dampr = TanH((1. - lx) / 40.0);

    result[i] = (1.0 - dampl * dampr) * f[i];
  }

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhxl(const Field2D &UNUSED(f0), const Field3D &f, BoutReal swidth,
                          BoutReal slength, bool UNUSED(BoutRealspace)) {
  TRACE("sink_tanhx");

  Field3D result;

  result.allocate();

  for (const auto &i : result) {
    BoutReal lx = mesh->GlobalX(i.x) - slength;
    BoutReal dampl = TanH(lx / swidth);

    result[i] = 0.5 * (1.0 - dampl) * f[i];
  }

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhxr(const Field2D &UNUSED(f0), const Field3D &f, BoutReal swidth,
                          BoutReal slength, bool UNUSED(BoutRealspace)) {
  TRACE("sink_tanhxr");
  Field3D result;
  result.allocate();

  for (const auto &i : result) {
    BoutReal rlx = 1. - mesh->GlobalX(i.x) - slength;
    BoutReal dampr = TanH(rlx / swidth);

    result[i] = 0.5 * (1.0 - dampr) * f[i];
  }

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to damp Psi to zero near radial boundaries
const Field3D buff_x(const Field3D &f, bool UNUSED(BoutRealspace)) {
  TRACE("buff_x");

  Field3D result;
  result.allocate();

  for (const auto &i : result) {
    BoutReal lx = mesh->GlobalX(i.x);
    BoutReal rlx = 1. - lx;
    BoutReal dampl = 1.e0;
    BoutReal dampr = 1.e0;
    BoutReal deltal = 0.05;
    BoutReal deltar = 0.05;

    result[i] = (dampl * exp(-(lx * lx) / (deltal * deltal)) +
                 dampr * exp(-(rlx * rlx) / (deltar * deltar))) *
                f[i];
  }

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}
