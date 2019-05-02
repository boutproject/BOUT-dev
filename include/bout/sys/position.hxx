#pragma once

#include "bout/region.hxx"

#if defined(__cplusplus) && __cplusplus >= 201703L
# define BOUT_FALL_THROUGH [[fallthrough]];
#elif defined(__GNUC__) && __GNUC__ >= 7
# define BOUT_FALL_THROUGH __attribute__ ((fallthrough));
#else
# define BOUT_FALL_THROUGH
#endif


class BoundaryRegion;
class Mesh;

class Position{
public:
  Position(Ind2D i, CELL_LOC loc, Mesh* msh, BoutReal _t)
    : Position(i.x(), i.y(), i.z(), loc, msh, _t){
    fz = INVALID;
  };
  Position(Ind3D i, CELL_LOC loc, Mesh * msh, BoutReal _t)
    : Position(i.x(), i.y(), i.z(), loc, msh, _t){};
  Position(int ix, int iy, int iz, CELL_LOC loc, Mesh * msh, BoutReal _t)
    : ix(ix), iy(iy), iz(iz), _t(_t), msh(msh), sx(loc == CELL_XLOW), sy(loc == CELL_YLOW),
      sz(loc == CELL_ZLOW){};
  Position(): ix(-1), iy(-1), iz(-1), _t(0), msh(nullptr), _x(0), _y(0), _z(0), fx(INVALID), fy(INVALID), fz(INVALID){};
  Position(const BoundaryRegion* bndry, int iz, CELL_LOC loc, BoutReal t, Mesh * msh);
  Position(const BoundaryRegion* bndry, CELL_LOC loc, BoutReal t, Mesh * msh): Position(bndry, 0, loc, t, msh) {};
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
  enum flags { DEFAULT, VALUE, ONE, TWOPI, REAL, INVALID};
  int ix;
  int iy;
  int iz;
  BoutReal _t;
  Mesh * msh;
  bool sx = false;
  bool sy = false;
  bool sz = false;
  BoutReal _x = BoutNaN, _y = BoutNaN, _z = BoutNaN;
  flags fx = DEFAULT, fy = DEFAULT, fz = DEFAULT;
};
