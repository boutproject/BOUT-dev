#pragma once

#include "bout/region.hxx"

#if defined(__GNUC__)
# define BOUT_FALL_THROUGH __attribute__ ((fallthrough));
#else
# define BOUT_FALL_THROUGH
#endif


class BoundaryRegion;
class Mesh;

class Position{
public:
  template <IND_TYPE N>
  Position(SpecificInd<N> i, CELL_LOC loc, Mesh * msh, BoutReal _t)
    : Position(i.x(), i.y(), i.z(), loc, msh, _t){};
  Position(int ix, int iy, int iz, CELL_LOC loc, Mesh * msh, BoutReal _t)
    : ix(ix), iy(iy), iz(iz), _t(_t), msh(msh), sx(loc == CELL_XLOW), sy(loc == CELL_YLOW),
      sz(loc == CELL_ZLOW){};
  Position(): valid(false) {};
  Position(const BoundaryRegion* bndry, int iz, CELL_LOC loc, BoutReal _t);
  Position(const BoundaryRegion* bndry, CELL_LOC loc, BoutReal _t): Position(bndry, 0, loc, _t) {};
  virtual BoutReal x();
  virtual BoutReal y();
  virtual BoutReal z();
  virtual BoutReal t(){
    return _t;
  };
  virtual int getIx() {
    return ix;
  };
  virtual int getIy() {
    return iy;
  };
  virtual int getIz() {
    return iz;
  }
  void setX(BoutReal val){
    _x = val;
    fx = VALUE;
  }
  void setY(BoutReal val){
    _y = val;
    fy = VALUE;
  }
  void setZ(BoutReal val){
    _z = val;
    fz = VALUE;
  }
private:
  enum flags { DEFAULT, VALUE, ONE, TWOPI, REAL};
  int ix;
  int iy;
  int iz;
  BoutReal _t;
  Mesh * msh;
  bool valid = true; // ignored
  bool sx = false;
  bool sy = false;
  bool sz = false;
  BoutReal _x = BoutNaN, _y = BoutNaN, _z = BoutNaN;
  flags fx = DEFAULT, fy = DEFAULT, fz = DEFAULT;
};
