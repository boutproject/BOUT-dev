/// Some standard boundary conditions

#ifndef __BNDRY_STD_H__
#define __BNDRY_STD_H__

#include "boundary_op.hxx"
#include "bout_types.hxx"
#include <field_factory.hxx>
#include "unused.hxx"

#include <utility>

/// Dirichlet boundary condition set half way between guard cell and grid cell at 2nd order accuracy
class BoundaryDirichlet_2ndOrder : public BoundaryOp {
 public:
  BoundaryDirichlet_2ndOrder() : val(0.) {}
  BoundaryDirichlet_2ndOrder(BoutReal setval ): val(setval) {}
  BoundaryDirichlet_2ndOrder(BoundaryRegion *region, BoutReal setval=0.):BoundaryOp(region),val(setval) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  BoutReal val;
};

/// Dirichlet (set to zero) boundary condition
class BoundaryDirichlet : public BoundaryOp {
 public:
  BoundaryDirichlet() : gen(nullptr) {}
  BoundaryDirichlet(BoundaryRegion *region, std::shared_ptr<FieldGenerator> g)
      : BoundaryOp(region), gen(std::move(g)) {}

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field2D &f,BoutReal t) override;
  void apply(Field3D &f) override;
  void apply(Field3D &f,BoutReal t) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  std::shared_ptr<FieldGenerator>  gen; // Generator
};

BoutReal default_func(BoutReal t, int x, int y, int z);

/// 3nd-order boundary condition
class BoundaryDirichlet_O3 : public BoundaryOp {
 public:
  BoundaryDirichlet_O3() : gen(nullptr) {}
  BoundaryDirichlet_O3(BoundaryRegion *region, std::shared_ptr<FieldGenerator> g)
      : BoundaryOp(region), gen(std::move(g)) {}

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field2D &f,BoutReal t) override;
  void apply(Field3D &f) override;
  void apply(Field3D &f,BoutReal t) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  std::shared_ptr<FieldGenerator>  gen; // Generator
};

/// 4th-order boundary condition
class BoundaryDirichlet_O4 : public BoundaryOp {
 public:
  BoundaryDirichlet_O4() : gen(nullptr) {}
  BoundaryDirichlet_O4(BoundaryRegion *region, std::shared_ptr<FieldGenerator> g)
      : BoundaryOp(region), gen(std::move(g)) {}

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field2D &f,BoutReal t) override;
  void apply(Field3D &f) override;
  void apply(Field3D &f,BoutReal t) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  std::shared_ptr<FieldGenerator>  gen; // Generator
};

/// Dirichlet boundary condition set half way between guard cell and grid cell at 4th order accuracy
class BoundaryDirichlet_4thOrder : public BoundaryOp {
 public:
  BoundaryDirichlet_4thOrder() : val(0.) {}
  BoundaryDirichlet_4thOrder(BoutReal setval ): val(setval) {}
  BoundaryDirichlet_4thOrder(BoundaryRegion *region, BoutReal setval=0.):BoundaryOp(region),val(setval) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  BoutReal val;
};

/// Neumann (zero-gradient) boundary condition for non-orthogonal meshes
class BoundaryNeumann_NonOrthogonal : public BoundaryOp {
 public:
  BoundaryNeumann_NonOrthogonal(): val(0.) {}
  BoundaryNeumann_NonOrthogonal(BoutReal setval ): val(setval) {}
  BoundaryNeumann_NonOrthogonal(BoundaryRegion *region, BoutReal setval=0.):BoundaryOp(region),val(setval) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;
 private:
  BoutReal val;
};

/// Neumann (zero-gradient) boundary condition, using 2nd order on boundary
class BoundaryNeumann2 : public BoundaryOp {
 public:
  BoundaryNeumann2() {}
  BoundaryNeumann2(BoundaryRegion *region):BoundaryOp(region) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;
};

/// Neumann boundary condition set half way between guard cell and grid cell at 2nd order accuracy
class BoundaryNeumann_2ndOrder : public BoundaryOp {
 public:
  BoundaryNeumann_2ndOrder() : val(0.) {}
  BoundaryNeumann_2ndOrder(BoutReal setval ): val(setval) {}
  BoundaryNeumann_2ndOrder(BoundaryRegion *region, BoutReal setval=0.):BoundaryOp(region),val(setval) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  BoutReal val;
};

// Neumann boundary condition set half way between guard cell and grid cell at 2nd order accuracy
class BoundaryNeumann : public BoundaryOp {
 public:
  BoundaryNeumann() : gen(nullptr) {}
  BoundaryNeumann(BoundaryRegion *region, std::shared_ptr<FieldGenerator> g)
      : BoundaryOp(region), gen(std::move(g)) {}

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field2D &f, BoutReal t) override;
  void apply(Field3D &f) override;
  void apply(Field3D &f,BoutReal t) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Neumann boundary condition set half way between guard cell and grid cell at 4th order accuracy
class BoundaryNeumann_4thOrder : public BoundaryOp {
 public:
  BoundaryNeumann_4thOrder() : val(0.) {}
  BoundaryNeumann_4thOrder(BoutReal setval ): val(setval) {}
  BoundaryNeumann_4thOrder(BoundaryRegion *region, BoutReal setval=0.):BoundaryOp(region),val(setval) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  BoutReal val;
};

/// Neumann boundary condition set half way between guard cell and grid cell at 4th order accuracy
class BoundaryNeumann_O4 : public BoundaryOp {
 public:
  BoundaryNeumann_O4() : gen(nullptr) {}
  BoundaryNeumann_O4(BoundaryRegion *region, std::shared_ptr<FieldGenerator> g)
      : BoundaryOp(region), gen(std::move(g)) {}

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field2D &f, BoutReal t) override;
  void apply(Field3D &f) override;
  void apply(Field3D &f,BoutReal t) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  std::shared_ptr<FieldGenerator> gen;
};

/// NeumannPar (zero-gradient) boundary condition on
/// the variable / sqrt(g_22)
class BoundaryNeumannPar : public BoundaryOp {
 public:
  BoundaryNeumannPar() {}
  BoundaryNeumannPar(BoundaryRegion *region):BoundaryOp(region) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;
};

/// Robin (mix of Dirichlet and Neumann)
class BoundaryRobin : public BoundaryOp {
 public:
  BoundaryRobin() : aval(0.), bval(0.), gval(0.) {}
  BoundaryRobin(BoundaryRegion *region, BoutReal a, BoutReal b, BoutReal g):BoundaryOp(region), aval(a), bval(b), gval(g) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;
private:
  BoutReal aval, bval, gval;
};

/// Constant gradient (zero second derivative)
class BoundaryConstGradient : public BoundaryOp {
 public:
  BoundaryConstGradient() {}
  BoundaryConstGradient(BoundaryRegion *region):BoundaryOp(region) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;
};

/// Zero Laplacian, decaying solution
class BoundaryZeroLaplace : public BoundaryOp {
 public:
  BoundaryZeroLaplace() {}
  BoundaryZeroLaplace(BoundaryRegion *region):BoundaryOp(region) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;
};

/// Zero Laplacian
class BoundaryZeroLaplace2 : public BoundaryOp {
 public:
  BoundaryZeroLaplace2() {}
  BoundaryZeroLaplace2(BoundaryRegion *region):BoundaryOp(region) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;
};

/// Constant Laplacian, decaying solution
class BoundaryConstLaplace : public BoundaryOp {
 public:
  BoundaryConstLaplace() {}
  BoundaryConstLaplace(BoundaryRegion *region):BoundaryOp(region) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;
};

/// Vector boundary condition Div(B) = 0, Curl(B) = 0
class BoundaryDivCurl : public BoundaryOp {
 public:
  BoundaryDivCurl() {}
  BoundaryDivCurl(BoundaryRegion *region):BoundaryOp(region) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &UNUSED(f)) override { throw BoutException("ERROR: DivCurl boundary only for vectors"); }
  void apply(Field3D &UNUSED(f)) override { throw BoutException("ERROR: DivCurl boundary only for vectors"); }
  void apply(Vector2D &f) override;
  void apply(Vector3D &f) override;
};

/// Free boundary condition (evolve the field in the guard cells, using non-centred derivatives to calculate the ddt)
class BoundaryFree : public BoundaryOp {
 public:
  BoundaryFree() : val(0.) {apply_to_ddt = true;}
  BoundaryFree(BoutReal setval): val(setval) {}
  BoundaryFree(BoundaryRegion *region, BoutReal setval=0.):BoundaryOp(region),val(setval) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  BoutReal val;
};

// L. Easy
/// Alternative free boundary condition (evolve the field in the guard cells, using non-centred derivatives to calculate the ddt)
class BoundaryFree_O2 : public BoundaryOp {
public:
  BoundaryFree_O2()  {}
  BoundaryFree_O2(BoundaryRegion *region):BoundaryOp(region) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
};

class BoundaryFree_O3 : public BoundaryOp {
public:
  BoundaryFree_O3() {}
  BoundaryFree_O3(BoundaryRegion *region):BoundaryOp(region) { }

  using BoundaryOp::clone;
  BoundaryOp* clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override;
  void apply(Field3D &f) override;

  using BoundaryOp::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;

};
// End L.Easy


/////////////////////////////////////////////////////////

/// Convert a boundary condition to a relaxing one
class BoundaryRelax : public BoundaryModifier {
 public:
  BoundaryRelax() : r(10.) {apply_to_ddt = true;}  // Set default rate
  BoundaryRelax(BoundaryOp *operation, BoutReal rate) : BoundaryModifier(operation) {r = fabs(rate); apply_to_ddt = true;}
  BoundaryOp* cloneMod(BoundaryOp *op, const std::list<std::string> &args) override;

  using BoundaryModifier::apply;
  void apply(Field2D &f) override {apply(f, 0.);};
  void apply(Field2D &f, BoutReal t) override;
  void apply(Field3D &f) override {apply(f, 0.);};
  void apply(Field3D &f, BoutReal t) override;

  using BoundaryModifier::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
 private:
  BoutReal r;
};

/// Increase the width of a boundary
class BoundaryWidth : public BoundaryModifier {
public:
  BoundaryWidth() : width(2) {}
  BoundaryWidth(BoundaryOp *operation, int wid) : BoundaryModifier(operation), width(wid) {}
  BoundaryOp* cloneMod(BoundaryOp *op, const std::list<std::string> &args) override;

  using BoundaryModifier::apply;
  void apply(Field2D &f) override {apply(f, 0.);};
  void apply(Field2D &f, BoutReal t) override;
  void apply(Field3D &f) override {apply(f, 0.);};
  void apply(Field3D &f, BoutReal t) override;

  using BoundaryModifier::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
private:
  int width;
};

/// Convert input field fromFieldAligned, apply boundary and then convert back toFieldAligned
/// Equivalent to converting the boundary condition to "Field Aligned" from "orthogonal"
class BoundaryToFieldAligned : public BoundaryModifier {
public:
  BoundaryToFieldAligned(){}
  BoundaryToFieldAligned(BoundaryOp *operation) : BoundaryModifier(operation){}
  BoundaryOp* cloneMod(BoundaryOp *op, const std::list<std::string> &args) override;

  using BoundaryModifier::apply;
  void apply(Field2D &f) override {apply(f, 0.);};
  void apply(Field2D &f, BoutReal t) override;
  void apply(Field3D &f) override {apply(f, 0.);};
  void apply(Field3D &f, BoutReal t) override;

  using BoundaryModifier::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
private:
};

/// Convert input field toFieldAligned, apply boundary and then convert back fromFieldAligned
/// Equivalent to converting the boundary condition from "Field Aligned" to "orthogonal"
class BoundaryFromFieldAligned : public BoundaryModifier {
public:
  BoundaryFromFieldAligned(){}
  BoundaryFromFieldAligned(BoundaryOp *operation) : BoundaryModifier(operation){}
  BoundaryOp* cloneMod(BoundaryOp *op, const std::list<std::string> &args) override;

  using BoundaryModifier::apply;
  void apply(Field2D &f) override {apply(f, 0.);};
  void apply(Field2D &f, BoutReal t) override;
  void apply(Field3D &f) override {apply(f, 0.);};
  void apply(Field3D &f, BoutReal t) override;

  using BoundaryModifier::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
private:
};

#endif // __BNDRY_STD_H__
