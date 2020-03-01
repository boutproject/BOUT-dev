
class BoundaryRegion;

#ifndef __BNDRY_REGION_H__
#define __BNDRY_REGION_H__

#include <string>
#include <utility>

class Mesh;
namespace bout {
namespace globals {
  extern Mesh* mesh; ///< Global mesh
} // namespace bout
} // namespace globals

/// Location of boundary
enum class BndryLoc {xin,
                     xout,
                     ydown,
                     yup,
                     all,
                     par_fwd,   // Don't include parallel boundaries
                     par_bkwd};
constexpr BndryLoc BNDRY_XIN = BndryLoc::xin;
constexpr BndryLoc BNDRY_XOUT = BndryLoc::xout;
constexpr BndryLoc BNDRY_YDOWN = BndryLoc::ydown;
constexpr BndryLoc BNDRY_YUP = BndryLoc::yup;
constexpr BndryLoc BNDRY_ALL = BndryLoc::all;
constexpr BndryLoc BNDRY_PAR_FWD = BndryLoc::par_fwd;
constexpr BndryLoc BNDRY_PAR_BKWD = BndryLoc::par_bkwd;

class BoundaryRegionBase {
public:

  BoundaryRegionBase() = delete;
  BoundaryRegionBase(std::string name, Mesh *passmesh = nullptr)
      : localmesh(passmesh ? passmesh : bout::globals::mesh), label(std::move(name)) {}
  BoundaryRegionBase(std::string name, BndryLoc loc, Mesh *passmesh = nullptr)
      : localmesh(passmesh ? passmesh : bout::globals::mesh), label(std::move(name)), location(loc) {}

  virtual ~BoundaryRegionBase() = default;

  Mesh* localmesh; ///< Mesh does this boundary region belongs to

  std::string label; ///< Label for this boundary region

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
  BoundaryRegion() = delete;
  BoundaryRegion(std::string name, BndryLoc loc, Mesh *passmesh = nullptr)
      : BoundaryRegionBase(name, loc, passmesh) {}
  BoundaryRegion(std::string name, int xd, int yd, Mesh *passmesh = nullptr)
      : BoundaryRegionBase(name, passmesh), bx(xd), by(yd), width(2) {}
  ~BoundaryRegion() override = default;

  int x,y; ///< Indices of the point in the boundary
  int bx, by; ///< Direction of the boundary [x+bx][y+by] is going outwards

  int width; ///< Width of the boundary

  virtual void next1d() = 0; ///< Loop over the innermost elements
  virtual void nextX() = 0;  ///< Just loop over X
  virtual void nextY() = 0;  ///< Just loop over Y
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
