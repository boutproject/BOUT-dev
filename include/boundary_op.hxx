
class BoundaryOp;
class BoundaryModifier;

#ifndef __BNDRY_OP__
#define __BNDRY_OP__

#include "boundary_region.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "unused.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

#include <cmath>
#include <list>
#include <map>
#include <string>
using std::string;
using std::list;

/// An operation on a boundary
class BoundaryOp {
public:
  BoundaryOp(bool apply_ddt = false)
      : bndry(nullptr), apply_to_ddt(apply_ddt), val(0.), gen(nullptr) {}
  BoundaryOp(BoundaryRegion *region, bool apply_ddt = false)
      : bndry(region), apply_to_ddt(apply_ddt), val(0.), gen(nullptr) {}
  BoundaryOp(BoundaryRegion *region, BoutReal val_in, std::shared_ptr<FieldGenerator> g)
      : bndry(region), apply_to_ddt(false), val(val_in), gen(std::move(g)) {}
  virtual ~BoundaryOp() {}

  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOp *clone(BoundaryRegion *UNUSED(region),
                            const list<string> &UNUSED(args),
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
  virtual void apply_ddt(Vector2D &f) { apply(ddt(f)); }
  virtual void apply_ddt(Vector3D &f) { apply(ddt(f)); }

  BoundaryRegion *bndry;
  const bool apply_to_ddt; // True if this boundary condition should be applied on the
                           // time derivatives, false if it should be applied to the field
                           // values
protected:
  const BoutReal val;                  // constant value for boundary condition
  std::shared_ptr<FieldGenerator> gen; // Generator
};

/// An operation on a boundary
template <typename Derived, bool needs_delta = false>
class BoundaryOpWithApply : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp;

  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOp *clone(BoundaryRegion *UNUSED(region),
                            const list<string> &UNUSED(args)) {
    ASSERT1(false); // this implementation should never get called
    return nullptr;
  }

  /// Apply a boundary condition on field f
  void apply(Field2D &f, BoutReal t = 0.) override { applyTemplate(f, t); }
  void apply(Field3D &f, BoutReal t = 0.) override { applyTemplate(f, t); }

private:
  template <typename T> void applyTemplate(T &f, BoutReal t);
};

class BoundaryModifier : public BoundaryOp {
public:
  BoundaryModifier(bool apply_ddt = false) : BoundaryOp(apply_ddt), op(nullptr) {}
  BoundaryModifier(BoundaryOp *operation, bool apply_ddt = false)
      : BoundaryOp(operation->bndry, apply_ddt), op(operation) {}
  virtual BoundaryOp *cloneMod(BoundaryOp *op, const list<string> &args) = 0;
  virtual BoundaryOpPar *cloneMod(BoundaryOpPar *UNUSED(op),
                                  const list<string> &UNUSED(args)) {
    throw BoutException("BoundaryModifier should not be called on a BoundaryOpPar.");
  }

protected:
  BoundaryOp *op;
};

#endif // __BNDRY_OP__
