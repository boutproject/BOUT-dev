
class BoundaryOp;
class BoundaryModifier;

#ifndef __BNDRY_OP__
#define __BNDRY_OP__

#include "boundary_region.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"
#include "unused.hxx"

#include <cmath>
#include <string>
#include <list>
#include <map>
using std::string;
using std::list;

/// An operation on a boundary
class BoundaryOp {
public:
  BoundaryOp()  : bndry(nullptr), apply_to_ddt(false), gen(nullptr) {}
  BoundaryOp(BoundaryRegion *region)
    : bndry(region), apply_to_ddt(false), gen(nullptr) {}
  BoundaryOp(BoundaryRegion *region, std::shared_ptr<FieldGenerator> g)
    : bndry(region), apply_to_ddt(false), gen(std::move(g)) {}
  virtual ~BoundaryOp() {}

  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOp* clone(BoundaryRegion *UNUSED(region), const list<string> &UNUSED(args)) {
    throw BoutException("BoundaryOp::clone not implemented");
  }

  /// Clone using positional args and keywords
  /// If not implemented, check if keywords are passed, then call two-argument version
  virtual BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                            const std::map<std::string, std::string> &keywords) {
    if (!keywords.empty()) {
      // Given keywords, but not using
      throw BoutException("Keywords ignored in boundary : %s", keywords.begin()->first.c_str());
    }
    
    return clone(region, args);
  }


  /// Apply a boundary condition on field f
  virtual void apply(Field2D &f,BoutReal t = 0.) {
    applyTemplate(f, t);
  }
  virtual void apply(Field3D &f,BoutReal t = 0.) {
    applyTemplate(f, t);
  }

  virtual void apply(Vector2D &f) {
    apply(f.x);
    apply(f.y);
    apply(f.z);
  }

  virtual void apply(Vector3D &f) {
    apply(f.x);
    apply(f.y);
    apply(f.z);
  }

  /// Apply a boundary condition on ddt(f)
  virtual void apply_ddt(Field2D &f);
  virtual void apply_ddt(Field3D &f);
  virtual void apply_ddt(Vector2D &f) {
    apply(ddt(f));
  }
  virtual void apply_ddt(Vector3D &f) {
    apply(ddt(f));
  }

  BoundaryRegion *bndry;
  bool apply_to_ddt; // True if this boundary condition should be applied on the time derivatives, false if it should be applied to the field values

protected:
  std::shared_ptr<FieldGenerator> gen; // Generator

  // Apply boundary condition at a point
  virtual void applyAtPoint(Field2D &UNUSED(f), BoutReal UNUSED(val), int
      UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z),
      Coordinates* UNUSED(metric)) {
    throw BoutException("BoundaryOp::applyAtPoint() should never be called. A "
        "subclass should either override BoundaryOp::apply() or override "
        "applyAtPoint() and applyAtPointStaggered().");
  }
  virtual void applyAtPoint(Field3D &UNUSED(f), BoutReal UNUSED(val), int
      UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z),
      Coordinates* UNUSED(metric)) {
    throw BoutException("BoundaryOp::applyAtPoint() should never be called. A "
        "subclass should either override BoundaryOp::apply() or override "
        "applyAtPoint() and applyAtPointStaggered().");
  }

  // Apply to staggered grid
  virtual void applyAtPointStaggered(Field2D &UNUSED(f), BoutReal UNUSED(val),
      int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int
      UNUSED(z), Coordinates* UNUSED(metric)) {
    throw BoutException("BoundaryOp::applyAtPointStaggered() should never be "
        "called. A subclass should either override BoundaryOp::apply() or "
        "override applyAtPoint() and applyAtPointStaggered().");
  }
  virtual void applyAtPointStaggered(Field3D &UNUSED(f), BoutReal UNUSED(val),
      int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int
      UNUSED(z), Coordinates* UNUSED(metric)) {
    throw BoutException("BoundaryOp::applyAtPointStaggered() should never be "
        "called. A subclass should either override BoundaryOp::apply() or "
        "override applyAtPoint() and applyAtPointStaggered().");
  }

  // extrapolate to further guard cells
  virtual void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  virtual void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);

private:
  template<typename T>
  void applyTemplate(T &f, BoutReal t);
};

class BoundaryModifier : public BoundaryOp {
public:
  BoundaryModifier() : op(nullptr) {}
  BoundaryModifier(BoundaryOp *operation) : BoundaryOp(operation->bndry), op(operation) {}
  virtual BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args) = 0;
  virtual BoundaryOpPar* cloneMod(BoundaryOpPar *UNUSED(op), const list<string> &UNUSED(args)) {
    throw BoutException("BoundaryModifier should not be called on a BoundaryOpPar.");
  }
protected:
  BoundaryOp *op;
};

#endif // __BNDRY_OP__
