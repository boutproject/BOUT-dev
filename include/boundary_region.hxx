
class BoundaryRegion;

#ifndef __BNDRY_REGION_H__
#define __BNDRY_REGION_H__

#include <string>
using std::string;

/// Location of boundary
enum BndryLoc {BNDRY_XIN=1,
               BNDRY_XOUT=2,
               BNDRY_YDOWN=4,
               BNDRY_YUP=8,
               BNDRY_ALL=15,
               BNDRY_PAR_FWD=16,   // Don't include parallel boundaries
               BNDRY_PAR_BKWD=32};

class BoundaryRegionBase {
public:
  BoundaryRegionBase() {}
  BoundaryRegionBase(const string &name) : label(name) {}
  BoundaryRegionBase(const string &name, BndryLoc loc) : label(name), location(loc) {}
  virtual ~BoundaryRegionBase() {}

  string label; ///< Label for this boundary region

  BndryLoc location;         ///< Which side of the domain is it on?
  bool isParallel = false;   ///< Is this a parallel boundary?

  virtual void first() = 0;  ///< Move the region iterator to the start
  virtual void next() = 0;   ///< Get the next element in the loop
                             ///  over every element from inside out (in
                             ///  X or Y first)
  virtual bool isDone() = 0; ///< Returns true if outside domain. Can use this with nested nextX, nextY
};

/// Describes a region of the boundary, and a means of iterating over it
class BoundaryRegion : public BoundaryRegionBase {
public:
  BoundaryRegion() {}
  BoundaryRegion(const string &name, BndryLoc loc) : BoundaryRegionBase(name, loc) {}
  BoundaryRegion(const string &name, int xd, int yd) : BoundaryRegionBase(name), bx(xd), by(yd), width(2) {}
  virtual ~BoundaryRegion() {}

  int x,y; ///< Indices of the point in the boundary
  int bx, by; ///< Direction of the boundary [x+dx][y+dy] is going outwards

  int width; ///< Width of the boundary

  virtual void next1d() = 0; ///< Loop over the innermost elements
  virtual void nextX() = 0;  ///< Just loop over X
  virtual void nextY() = 0;  ///< Just loop over Y
};

class BoundaryRegionXIn : public BoundaryRegion {
public:
  BoundaryRegionXIn(const string &name, int ymin, int ymax);

  void first();
  void next();
  void next1d();
  void nextX();
  void nextY();
  bool isDone();
private:
  int ys, ye;
};

class BoundaryRegionXOut : public BoundaryRegion {
public:
  BoundaryRegionXOut(const string &name, int ymin, int ymax);

  void first();
  void next();
  void next1d();
  void nextX();
  void nextY();
  bool isDone();
private:
  int ys, ye;
};

class BoundaryRegionYDown : public BoundaryRegion {
public:
  BoundaryRegionYDown(const string &name, int xmin, int xmax);

  void first();
  void next();
  void next1d();
  void nextX();
  void nextY();
  bool isDone();
private:
  int xs, xe;
};

class BoundaryRegionYUp : public BoundaryRegion {
public:
  BoundaryRegionYUp(const string &name, int xmin, int xmax);

  void first();
  void next();
  void next1d();
  void nextX();
  void nextY();
  bool isDone();
private:
  int xs, xe;
};

#endif // __BNDRY_REGION_H__
