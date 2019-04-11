#include "boundary_region.hxx"
#include "bout/sys/position.hxx"
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


BoutReal Position::x(){
  switch (fx) {
  case DEFAULT:
    if (sx) {
      return 0.5 * (msh->GlobalX(ix) + msh->GlobalX(ix-1));
    } else {
      return msh->GlobalX(ix);
    }
    break;
  case VALUE:
    return _x;
    }

}
