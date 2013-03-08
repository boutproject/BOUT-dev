/// Some standard boundary conditions

#ifndef __BNDRY_STD_H__
#define __BNDRY_STD_H__

#include "boundary_op.hxx"
#include "bout_types.hxx"

/// Dirichlet (set to zero) boundary condition
class BoundaryDirichlet : public BoundaryOp {
 public:
  BoundaryDirichlet() : val(0.) {}
  BoundaryDirichlet(const BoutReal setval): val(setval) {}
  BoundaryDirichlet(BoundaryRegion *region, BoutReal setval=0.):BoundaryOp(region),val(setval) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f);
  void apply(Field3D &f);
  
  void apply_ddt(Field2D &f);
  void apply_ddt(Field3D &f);
 private:
  BoutReal val;
};

/// Neumann (zero-gradient) boundary condition
class BoundaryNeumann : public BoundaryOp {
 public:
  BoundaryNeumann(): val(0.) {}
  BoundaryNeumann(const BoutReal setval): val(setval) {}
  BoundaryNeumann(BoundaryRegion *region, BoutReal setval=0.):BoundaryOp(region),val(setval) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f);
  void apply(Field3D &f);
 private:
  BoutReal val;
};

/// Neumann (zero-gradient) boundary condition, using 2nd order on boundary
class BoundaryNeumann2 : public BoundaryOp {
 public:
  BoundaryNeumann2() {}
  BoundaryNeumann2(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f);
  void apply(Field3D &f);
};

/// NeumannPar (zero-gradient) boundary condition on
/// the variable / sqrt(g_22)
class BoundaryNeumannPar : public BoundaryOp {
 public:
  BoundaryNeumannPar() {}
  BoundaryNeumannPar(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f);
  void apply(Field3D &f);
};

/// Robin (mix of Dirichlet and Neumann)
class BoundaryRobin : public BoundaryOp {
 public:
  BoundaryRobin() : aval(0.), bval(0.), gval(0.) {}
  BoundaryRobin(BoundaryRegion *region, BoutReal a, BoutReal b, BoutReal g):BoundaryOp(region), aval(a), bval(b), gval(g) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f);
  void apply(Field3D &f);
private:
  BoutReal aval, bval, gval;
};

/// Constant gradient (zero second derivative)
class BoundaryConstGradient : public BoundaryOp {
 public:
  BoundaryConstGradient() {}
  BoundaryConstGradient(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f);
  void apply(Field3D &f);
};

/// Zero Laplacian, decaying solution
class BoundaryZeroLaplace : public BoundaryOp {
 public:
  BoundaryZeroLaplace() {}
  BoundaryZeroLaplace(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f);
  void apply(Field3D &f);
};

/// Zero Laplacian
class BoundaryZeroLaplace2 : public BoundaryOp {
 public:
  BoundaryZeroLaplace2() {}
  BoundaryZeroLaplace2(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f);
  void apply(Field3D &f);
};

/// Constant Laplacian, decaying solution
class BoundaryConstLaplace : public BoundaryOp {
 public:
  BoundaryConstLaplace() {}
  BoundaryConstLaplace(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f);
  void apply(Field3D &f);
};

/// Vector boundary condition Div(B) = 0, Curl(B) = 0
class BoundaryDivCurl : public BoundaryOp {
 public:
  BoundaryDivCurl() {}
  BoundaryDivCurl(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
  void apply(Field2D &f) { bout_error("ERROR: DivCurl boundary only for vectors"); }
  void apply(Field3D &f) { bout_error("ERROR: DivCurl boundary only for vectors"); }
  void apply(Vector2D &f);
  void apply(Vector3D &f);
};

/////////////////////////////////////////////////////////

/// Convert a boundary condition to a relaxing one
class BoundaryRelax : public BoundaryModifier {
 public:
  BoundaryRelax() : r(10.) {}  // Set default rate
  BoundaryRelax(BoundaryOp *operation, BoutReal rate) : BoundaryModifier(operation) {r = fabs(rate);}
  BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args);
  
  void apply(Field2D &f);
  void apply(Field3D &f);
  
  void apply_ddt(Field2D &f);
  void apply_ddt(Field3D &f);
 private:
  BoutReal r;
};

/// Apply boundary condition in shifted coordinates
class BoundaryShifted : public BoundaryModifier {
public:
  BoundaryShifted() {}
  BoundaryShifted(BoundaryOp *operation) : BoundaryModifier(operation) {}
  BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args);
  
  void apply(Field2D &f);
  void apply(Field3D &f);
  
  void apply_ddt(Field2D &f);
  void apply_ddt(Field3D &f);
private:
  
};

/// Increase the width of a boundary
class BoundaryWidth : public BoundaryModifier {
public:
  BoundaryWidth() : width(2) {}
  BoundaryWidth(BoundaryOp *operation, int wid) : BoundaryModifier(operation), width(wid) {}
  BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args);
  
  void apply(Field2D &f);
  void apply(Field3D &f);
  
  void apply_ddt(Field2D &f);
  void apply_ddt(Field3D &f);
private:
  int width;
};

#endif // __BNDRY_STD_H__
