
class BoundaryOp;
class BoundaryModifier;

#ifndef __BNDRY_OP__
#define __BNDRY_OP__

#include "boundary_region.h"
#include "field2d.h"
#include "field3d.h"
#include "vector2d.h"
#include "vector3d.h"
#include "bout.h"

#include <cmath>

/// An operation on a boundary
class BoundaryOp {
 public:
  BoundaryOp() {bndry = NULL;}
  BoundaryOp(BoundaryRegion *region) {bndry = region;}
  
  // Note: All methods must implement clone, except for modifiers (see below)
  virtual BoundaryOp* clone(BoundaryRegion *region) {return NULL;}
  
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
    apply_ddt(f.x);
    apply_ddt(f.y);
    apply_ddt(f.z);
  }
  virtual void apply_ddt(Vector3D &f) {
    apply_ddt(f.x);
    apply_ddt(f.y);
    apply_ddt(f.z);
  }
 protected:
  BoundaryRegion *bndry;
};

class BoundaryModifier : public BoundaryOp {
 public:
  virtual BoundaryOp* clone(BoundaryOp *op) = 0;
 protected:
  BoundaryOp *op;
};

#endif // __BNDRY_OP__
