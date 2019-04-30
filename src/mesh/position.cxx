#include "bout/sys/position.hxx"
#include "boundary_region.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"

Position::Position(const BoundaryRegion* bndry, int iz, CELL_LOC loc, BoutReal _t)
    : iz(iz), _t(_t), sx(loc == CELL_XLOW), sy(loc == CELL_YLOW), sz(loc == CELL_ZLOW) {
  ix = bndry->x;
  iy = bndry->y;
  if (bndry->bx) {
    if (sx && bndry->bx > 0) {
      sx = false;
      ix -= 1;
    } else if (sx && bndry->bx < 0) {
      sx = false;
    } else if (bndry->bx > 0) {
      sx = true;
    } else {
      sx = true;
      ix += 1;
    }
  }
  if (bndry->by) {
    if (sy && bndry->by > 0) {
      sy = false;
      iy -= 1;
    } else if (sy && bndry->by < 0) {
      sy = false;
    } else if (bndry->by > 0) {
      sy = true;
    } else {
      sy = true;
      iy += 1;
    }
  }
};

BoutReal Position::x() {
  BoutReal fac = 1;
  switch (fx) {
  case TWOPI:
    fac = 2 * PI;
    BOUT_FALL_THROUGH;
  case DEFAULT:
  case ONE:
    if (sx) {
      return fac * 0.5 * (msh->GlobalX(ix) + msh->GlobalX(ix - 1));
    } else {
      return fac * msh->GlobalX(ix);
    }
  case VALUE:
    return _x;
  case REAL:
  default:
    throw BoutException("Not implemented");
  }
}

BoutReal Position::y() {
  BoutReal fac = 2 * PI;
  switch (fy) {
  case ONE:
    fac = 1;
    BOUT_FALL_THROUGH;
  case TWOPI:
  case DEFAULT:
    if (sy) {
      return fac * 0.5 * (msh->GlobalY(iy) + msh->GlobalY(iy - 1));
    } else {
      return fac * msh->GlobalY(iy);
    }
  case VALUE:
    return _y;
  case REAL:
  default:
    throw BoutException("Not implemented");
  }
}
BoutReal Position::z() {
  BoutReal fac = 2 * PI;
  switch (fz) {
  case ONE:
    fac = 1;
    BOUT_FALL_THROUGH;
  case TWOPI:
  case DEFAULT:
    if (sz) {
      return fac * (iz - 0.5) / static_cast<BoutReal>(msh->LocalNz);
    } else {
      return fac * iz / static_cast<BoutReal>(msh->LocalNz);
    }
  case VALUE:
    return _z;
  case REAL:
  default:
    throw BoutException("Not implemented");
  }
}
