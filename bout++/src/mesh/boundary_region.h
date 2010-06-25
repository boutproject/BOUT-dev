
class BoundaryRegion;

#ifndef __BNDRY_REGION_H__
#define __BNDRY_REGION_H__

enum BndryLoc {BNDRY_XIN=1, BNDRY_XOUT=2, BNDRY_YDOWN=4, BNDRY_YUP=8, BNDRY_ALL=15};

/// Describes a region of the boundary, and a means of iterating over it
class BoundaryRegion {
 public:
  virtual ~BoundaryRegion() {}
  
  BndryLoc location; // Which side of the domain is it on?
  
  int x,y; // Indices of the point in the boundary
  int bx, by; // Direction of the boundary [x+dx][y+dy] is going outwards

  virtual void first() = 0;
  virtual void next() = 0; // Loop over every element from inside out (in X or Y first)
  virtual void nextX() = 0; // Just loop over X
  virtual void nextY() = 0; // Just loop over Y
  virtual bool isDone() = 0; // Returns true if outside domain. Can use this with nested nextX, nextY
};

class BoundaryRegionXIn : public BoundaryRegion {
 public:
  BoundaryRegionXIn(int ymin, int ymax);
  
  void first();
  void next();
  void nextX();
  void nextY();
  bool isDone();
 private:
  int ys, ye;
};

class BoundaryRegionXOut : public BoundaryRegion {
 public:
  BoundaryRegionXOut(int ymin, int ymax);
  
  void first();
  void next();
  void nextX();
  void nextY();
  bool isDone();
 private:
  int ys, ye;
};

class BoundaryRegionYDown : public BoundaryRegion {
 public:
  BoundaryRegionYDown(int xmin, int xmax);
  
  void first();
  void next();
  void nextX();
  void nextY();
  bool isDone();
 private:
  int xs, xe;
};

class BoundaryRegionYUp : public BoundaryRegion {
 public:
  BoundaryRegionYUp(int xmin, int xmax);
  
  void first();
  void next();
  void nextX();
  void nextY();
  bool isDone();
 private:
  int xs, xe;
};

#endif // __BNDRY_REGION_H__
