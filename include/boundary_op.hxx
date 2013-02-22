
class BoundaryOp;
class BoundaryModifier;

#ifndef __BNDRY_OP__
#define __BNDRY_OP__

#include "boundary_region.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

#include <cmath>
#include <string>
#include <list>
using std::string;
using std::list;

/// An operation on a boundary
class BoundaryOp {
 public:
  BoundaryOp() {bndry = NULL;}
  BoundaryOp(BoundaryRegion *region) {bndry = region;}
  virtual ~BoundaryOp() {}
  
  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOp* clone(BoundaryRegion *region, const list<string> &args) {return NULL; } 
  
  /// Apply a boundary condition on field f
  virtual void apply(Field2D &f) = 0;
  virtual void apply(Field3D &f) = 0;
  
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
};

class BoundaryModifier : public BoundaryOp {
 public:
  BoundaryModifier() : op(NULL) {}
  BoundaryModifier(BoundaryOp *operation) : BoundaryOp(operation->bndry), op(operation) {}
  virtual BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args) = 0;
 protected:
  BoundaryOp *op;
};

#endif // __BNDRY_OP__
