
class BoundaryRegion;

#ifndef __BNDRY_REGION_H__
#define __BNDRY_REGION_H__

#include <bout/assert.hxx>

#include <string>
#include <utility>
using std::string;

class Mesh;
extern Mesh* mesh;

/// Location of boundary
enum BndryLoc {BNDRY_XIN=1,
               BNDRY_XOUT=2,
               BNDRY_YDOWN=4,
               BNDRY_YUP=8,
               BNDRY_ALL=15,
               BNDRY_PAR_FWD=16,   // Don't include parallel boundaries
               BNDRY_PAR_BKWD=32};

/// Describes a region of the boundary, and a means of iterating over it
class BoundaryRegion {
public:
  BoundaryRegion() = delete;
  BoundaryRegion(std::string name, int xd, int yd, Mesh *passmesh = nullptr)
      : bx(xd), by(yd), width(2), localmesh(passmesh ? passmesh : mesh), label(std::move(name)) {
        ASSERT1(!(bx == 0 && by == 0) && !(bx != 0 && by != 0));
      }
  virtual ~BoundaryRegion() {}

  int x,y; ///< Indices of the point in the boundary
  const int bx, by; ///< Direction of the boundary [x+dx][y+dy] is going outwards

  int width; ///< Width of the boundary

  virtual void next1d() = 0; ///< Loop over the innermost elements
  virtual void nextX() = 0;  ///< Just loop over X
  virtual void nextY() = 0;  ///< Just loop over Y

  Mesh* localmesh; ///< Mesh does this boundary region belongs to

  string label; ///< Label for this boundary region

  BndryLoc location;         ///< Which side of the domain is it on?

  virtual void first() = 0;  ///< Move the region iterator to the start
  virtual void next() = 0;   ///< Get the next element in the loop
                             ///  over every element from inside out (in
                             ///  X or Y first)
  virtual bool isDone() = 0; ///< Returns true if outside domain. Can use this with nested nextX, nextY
};

class BoundaryRegionXIn : public BoundaryRegion {
public:
  BoundaryRegionXIn(std::string name, int ymin, int ymax, Mesh* passmesh = nullptr);

  void first() override;
  void next() override;
  void next1d() override;
  void nextX() override;
  void nextY() override;
  bool isDone() override;

private:
  int ys, ye;
};

class BoundaryRegionXOut : public BoundaryRegion {
public:
  BoundaryRegionXOut(std::string name, int ymin, int ymax, Mesh* passmesh = nullptr);

  void first() override;
  void next() override;
  void next1d() override;
  void nextX() override;
  void nextY() override;
  bool isDone() override;

private:
  int ys, ye;
};

class BoundaryRegionYDown : public BoundaryRegion {
public:
  BoundaryRegionYDown(std::string name, int xmin, int xmax, Mesh* passmesh = nullptr);

  void first() override;
  void next() override;
  void next1d() override;
  void nextX() override;
  void nextY() override;
  bool isDone() override;

private:
  int xs, xe;
};

class BoundaryRegionYUp : public BoundaryRegion {
public:
  BoundaryRegionYUp(std::string name, int xmin, int xmax, Mesh* passmesh = nullptr);

  void first() override;
  void next() override;
  void next1d() override;
  void nextX() override;
  void nextY() override;
  bool isDone() override;

private:
  int xs, xe;
};

#endif // __BNDRY_REGION_H__
