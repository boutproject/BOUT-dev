/// Some standard boundary conditions

#ifndef __BNDRY_STD_H__
#define __BNDRY_STD_H__

#include "boundary_op.h"
#include "bout_types.h"

/// Dirichlet (set to zero) boundary condition
class BoundaryDirichlet : public BoundaryOp {
 public:
  BoundaryDirichlet() {}
 BoundaryDirichlet(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region);
  void apply(Field2D &f);
  void apply(Field3D &f);
};

/// Neumann (zero-gradient) boundary condition
class BoundaryNeumann : public BoundaryOp {
 public:
  BoundaryNeumann() {}
 BoundaryNeumann(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region);
  void apply(Field2D &f);
  void apply(Field3D &f);
};

/// Convert a boundary condition to a relaxing one
class BoundaryRelax : public BoundaryModifier {
 public:
  BoundaryRelax(BoutReal rate) {r = fabs(rate);}
  BoundaryOp* clone(BoundaryOp *op);
  
  void apply(Field2D &f);
  void apply(Field3D &f);
  
  void apply_ddt(Field2D &f);
  void apply_ddt(Field3D &f);
 private:
  BoundaryRelax() {} // Must be initialised with a rate
  BoutReal r;
};

#endif // __BNDRY_STD_H__
