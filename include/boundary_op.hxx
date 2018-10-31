
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
  BoundaryOp(bool apply_ddt = false) : bndry(nullptr), apply_to_ddt(apply_ddt),
                                       val(0.), gen(nullptr), width(0) {}
  BoundaryOp(BoundaryRegion *region, int width_in = 0, bool apply_ddt = false)
    : bndry(region), apply_to_ddt(apply_ddt), val(0.), gen(nullptr),
      width(width_in ? width_in : region->width) {}
  BoundaryOp(BoundaryRegion *region, BoutReal val_in, std::shared_ptr<FieldGenerator> g, int width_in = 0)
    : bndry(region), apply_to_ddt(false), val(val_in), gen(std::move(g)),
      width(width_in ? width_in : region->width) {}
  virtual ~BoundaryOp() {}

  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOp *clone(BoundaryRegion *UNUSED(region), const list<string> &UNUSED(args),
                            const std::map<std::string, std::string> &UNUSED(keywords)) {
    throw BoutException("BoundaryOp::clone not implemented");
    return nullptr;
  }


  /// Apply a boundary condition on field f
  virtual void apply(Field2D &f, BoutReal t = 0.) = 0;
  virtual void apply(Field3D &f, BoutReal t = 0.) = 0;

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
  const bool apply_to_ddt; // True if this boundary condition should be applied on the time derivatives, false if it should be applied to the field values
protected:
  const BoutReal val; // constant value for boundary condition
  std::shared_ptr<FieldGenerator> gen; // Generator
  const int width; // boundary width, stored in case we change it from the default
};

/// An operation on a boundary
template<typename Derived, bool needs_delta = false>
class BoundaryOpWithApply : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp;

  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOp* clone(BoundaryRegion *UNUSED(region), const list<string> &UNUSED(args)) {
    ASSERT1(false); // this implementation should never get called
    return nullptr;
  }

  /// Apply a boundary condition on field f
  void apply(Field2D &f, BoutReal t = 0.) override {
    applyTemplate(f, t);
  }
  void apply(Field3D &f, BoutReal t = 0.) override {
    applyTemplate(f, t);
  }

protected:
  //// Apply boundary condition at a point
  //virtual void applyAtPoint(Field2D &UNUSED(f), BoutReal UNUSED(val), int
  //    UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z),
  //    Coordinates* UNUSED(metric)) = 0;
  //virtual void applyAtPoint(Field3D &UNUSED(f), BoutReal UNUSED(val), int
  //    UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int UNUSED(z),
  //    Coordinates* UNUSED(metric)) = 0;

  //// Apply to staggered grid
  //virtual void applyAtPointStaggered(Field2D &UNUSED(f), BoutReal UNUSED(val),
  //    int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int
  //    UNUSED(z), Coordinates* UNUSED(metric)) = 0;
  //virtual void applyAtPointStaggered(Field3D &UNUSED(f), BoutReal UNUSED(val),
  //    int UNUSED(x), int UNUSED(bx), int UNUSED(y), int UNUSED(by), int
  //    UNUSED(z), Coordinates* UNUSED(metric)) = 0;

  //// extrapolate to further guard cells
  //virtual void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z) = 0;
  //virtual void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z) = 0;

private:
  template<typename T>
  void applyTemplate(T &f, BoutReal t);
};

class BoundaryModifier : public BoundaryOp {
public:
  BoundaryModifier(bool apply_ddt = false) : BoundaryOp(apply_ddt), op(nullptr) {}
  BoundaryModifier(BoundaryOp *operation, bool apply_ddt = false)
    : BoundaryOp(operation->bndry, 0, apply_ddt), op(operation) {}
  virtual BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args) = 0;
  virtual BoundaryOpPar* cloneMod(BoundaryOpPar *UNUSED(op), const list<string> &UNUSED(args)) {
    throw BoutException("BoundaryModifier should not be called on a BoundaryOpPar.");
  }
protected:
  BoundaryOp *op;
};

#endif // __BNDRY_OP__
