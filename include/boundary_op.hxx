
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

class BoundaryOpBase {
public:
  BoundaryOpBase() {}
  virtual ~BoundaryOpBase() {}

  /// Apply a boundary condition on field f
  virtual void apply(Field2D &f) = 0;
  virtual void apply(Field2D &f,BoutReal UNUSED(t)){return apply(f);}//JMAD
  virtual void apply(Field3D &f) = 0;
  virtual void apply(Field3D &f,BoutReal UNUSED(t)){return apply(f);}//JMAD

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
};

/// An operation on a boundary
class BoundaryOp : public BoundaryOpBase {
public:
  BoundaryOp() {
    bndry = nullptr;
    apply_to_ddt = false;
  }
  BoundaryOp(BoundaryRegion *region) {bndry = region; apply_to_ddt=false;}
  ~BoundaryOp() override {}

  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOp* clone(BoundaryRegion *UNUSED(region), const list<string> &UNUSED(args)) {
    return nullptr;
  }
  
  virtual BoundaryOp* clone(BoundaryRegion *region, const list<string> &args, const std::map<std::string, std::string>& UNUSED(keywords)) {
    // If not implemented, call two-argument version
    return clone(region, args);
  }

  /// Apply a boundary condition on ddt(f)
  virtual void apply_ddt(Field2D &f) {
    apply(ddt(f));
  }
  virtual void apply_ddt(Field3D &f) {
    apply(ddt(f));
  }
  virtual void apply_ddt(Vector2D &f) {
    apply(ddt(f));
  }
  virtual void apply_ddt(Vector3D &f) {
    apply(ddt(f));
  }

  BoundaryRegion *bndry;
  bool apply_to_ddt; // True if this boundary condition should be applied on the time derivatives, false if it should be applied to the field values
};

class BoundaryModifier : public BoundaryOp {
public:
  BoundaryModifier() : op(nullptr) {}
  BoundaryModifier(BoundaryOp *operation) : BoundaryOp(operation->bndry), op(operation) {}
  virtual BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args) = 0;
protected:
  BoundaryOp *op;
};

#endif // __BNDRY_OP__
