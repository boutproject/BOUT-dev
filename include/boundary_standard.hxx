/// Some standard boundary conditions

#ifndef __BNDRY_STD_H__
#define __BNDRY_STD_H__

#include "boundary_op.hxx"
#include "bout_types.hxx"
#include <field_factory.hxx>
#include "unused.hxx"

#include <utility>

/// Dirichlet (set to zero) boundary condition
class BoundaryDirichlet : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  using BoundaryOp::apply;
  void apply(Field2D &f,BoutReal t = 0.) override {
    applyTemplate<Field2D>(f, t);
  }
  void apply(Field3D &f,BoutReal t = 0.) override {
    applyTemplate<Field3D>(f, t);
  }

 private:
  template<typename T>
  void applyTemplate(T &f, BoutReal t);
};

BoutReal default_func(BoutReal t, int x, int y, int z);

/// 3nd-order boundary condition
class BoundaryDirichlet_O3 : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
  void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z) override;
  void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z) override;
};

/// 4th-order boundary condition
class BoundaryDirichlet_O4 : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
  void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z) override;
  void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z) override;
};

/// Dirichlet boundary condition, tries to smooth out grid-scale oscillations at the boundary
class BoundaryDirichlet_smooth : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
};

/// Dirichlet boundary condition set half way between guard cell and grid cell at 2nd order accuracy
class BoundaryDirichlet_2ndOrder : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
};

/// Dirichlet boundary condition set half way between guard cell and grid cell at 4th order accuracy
class BoundaryDirichlet_O5 : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int UNUSED(bx), int y, int UNUSED(by), int z, Coordinates* UNUSED(metric)) override;
  void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z) override;
  void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z) override;
};

/// Neumann (zero-gradient) boundary condition for non-orthogonal meshes
class BoundaryNeumann_NonOrthogonal : public BoundaryOp {
 public:
  BoundaryNeumann_NonOrthogonal(): val(0.) {}
  BoundaryNeumann_NonOrthogonal(BoutReal setval ): val(setval) {}
  BoundaryNeumann_NonOrthogonal(BoundaryRegion *region, BoutReal setval=0.):BoundaryOp(region),val(setval) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  using BoundaryOp::apply;
  void apply(Field2D &f, BoutReal t = 0.) override {
    applyTemplate(f, t);
  }
  void apply(Field3D &f, BoutReal t = 0.) override {
    applyTemplate(f, t);
  }
 private:
  BoutReal val;

  template<typename T>
  void applyTemplate(T &f, BoutReal t);
};

/// Neumann (zero-gradient) boundary condition, using 2nd order on boundary
class BoundaryNeumann2 : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
};

/// Neumann boundary condition set half way between guard cell and grid cell at 2nd order accuracy
class BoundaryNeumann_2ndOrder : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
};

// Neumann boundary condition set half way between guard cell and grid cell at 2nd order accuracy
class BoundaryNeumann : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
};

/// Neumann boundary condition set half way between guard cell and grid cell at 4th order accuracy
class BoundaryNeumann_O4 : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z) override;
  void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z) override;
};

/// Neumann boundary condition set half way between guard cell and grid cell at 4th order accuracy
class BoundaryNeumann_4thOrder : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z) override;
  void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z) override;
};

/// NeumannPar (zero-gradient) boundary condition on
/// the variable / sqrt(g_22)
class BoundaryNeumannPar : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
};

/// Robin (mix of Dirichlet and Neumann)
class BoundaryRobin : public BoundaryOp {
 public:
  BoundaryRobin() : aval(0.), bval(0.), gval(0.) {}
  BoundaryRobin(BoundaryRegion *region, BoutReal a, BoutReal b, BoutReal g)
    : BoundaryOp(region), aval(a), bval(b), gval(g) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  using BoundaryOp::apply;
  void apply(Field2D &f, BoutReal t = 0.) override {
    applyTemplate(f, t);
  }
  void apply(Field3D &f, BoutReal t = 0.) override {
    applyTemplate(f, t);
  }
private:
  BoutReal aval, bval, gval;

  template<typename T>
  void applyTemplate(T &f, BoutReal t);
};

/// Constant gradient (zero second derivative)
class BoundaryConstGradient : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
};

/// Zero Laplacian, decaying solution
class BoundaryZeroLaplace : public BoundaryOp {
 public:
  BoundaryZeroLaplace() {}
  BoundaryZeroLaplace(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  using BoundaryOp::apply;
  void apply(Field2D &f, BoutReal UNUSED(t)) override;
  void apply(Field3D &f, BoutReal UNUSED(t)) override;
};

/// Zero Laplacian
class BoundaryZeroLaplace2 : public BoundaryOp {
 public:
  BoundaryZeroLaplace2() {}
  BoundaryZeroLaplace2(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  using BoundaryOp::apply;
  void apply(Field2D &f, BoutReal t) override;
  void apply(Field3D &f, BoutReal t) override;
};

/// Constant Laplacian, decaying solution
class BoundaryConstLaplace : public BoundaryOp {
 public:
  BoundaryConstLaplace() {}
  BoundaryConstLaplace(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  using BoundaryOp::apply;
  void apply(Field2D &f, BoutReal t) override;
  void apply(Field3D &f, BoutReal t) override;
};

/// Vector boundary condition Div(B) = 0, Curl(B) = 0
class BoundaryDivCurl : public BoundaryOp {
 public:
  BoundaryDivCurl() {}
  BoundaryDivCurl(BoundaryRegion *region):BoundaryOp(region) { }
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  using BoundaryOp::apply;
  void apply(Field2D &UNUSED(f), BoutReal UNUSED(t)) override { throw BoutException("ERROR: DivCurl boundary only for vectors"); }
  void apply(Field3D &UNUSED(f), BoutReal UNUSED(t)) override { throw BoutException("ERROR: DivCurl boundary only for vectors"); }
  void apply(Vector2D &f) override;
  void apply(Vector3D &f) override;
};

/// Free boundary condition (evolve the field in the guard cells, using non-centred derivatives to calculate the ddt)
class BoundaryFree : public BoundaryOp {
 public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  using BoundaryOp::apply;
  void apply(Field2D &f, BoutReal UNUSED(t)) override;
  void apply(Field3D &f, BoutReal UNUSED(t)) override;

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
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
};

class BoundaryFree_O3 : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z) override;
  void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z) override;
};
// End L.Easy

class BoundaryFree_O4 : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z) override;
  void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z) override;
};

class BoundaryFree_O5 : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp* clone(BoundaryRegion *region, const list<string> &args,
      const std::map<std::string, std::string> &keywords) override;

  void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* metric) override;
  void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z, Coordinates* UNUSED(metric)) override;
  void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z) override;
  void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z) override;
};

/////////////////////////////////////////////////////////

/// Convert a boundary condition to a relaxing one
class BoundaryRelax : public BoundaryModifier {
 public:
  BoundaryRelax() : BoundaryModifier(true), r(10.) {}  // Set default rate
  BoundaryRelax(BoundaryOp *operation, BoutReal rate) : BoundaryModifier(operation, true) { r = fabs(rate); }
  BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args) override;

  using BoundaryModifier::apply;
  void apply(Field2D &f, BoutReal t) override;
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
  BoundaryOp* cloneMod(BoundaryOp *UNUSED(op), const list<string> &UNUSED(args)) override {
    throw BoutException("WARNING: BoundaryWidth modifier is deprecated, use 'width' keyword to boundary conditions instead");
    return new BoundaryWidth(nullptr, 0);
  }

  void apply(Field2D &UNUSED(f), BoutReal UNUSED(t)) override {};
  void apply(Field3D &UNUSED(f), BoutReal UNUSED(t)) override {};

  void apply_ddt(Field2D &UNUSED(f)) override {};
  void apply_ddt(Field3D &UNUSED(f)) override {};
private:
  int width;
};

/// Convert input field fromFieldAligned, apply boundary and then convert back toFieldAligned
/// Equivalent to converting the boundary condition to "Field Aligned" from "orthogonal"
class BoundaryToFieldAligned : public BoundaryModifier {
public:
  BoundaryToFieldAligned(){}
  BoundaryToFieldAligned(BoundaryOp *operation) : BoundaryModifier(operation){}
  BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args) override;

  using BoundaryModifier::apply;
  void apply(Field2D &f, BoutReal t) override;
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
  BoundaryOp* cloneMod(BoundaryOp *op, const list<string> &args) override;

  using BoundaryModifier::apply;
  void apply(Field2D &f, BoutReal t) override;
  void apply(Field3D &f, BoutReal t) override;

  using BoundaryModifier::apply_ddt;
  void apply_ddt(Field2D &f) override;
  void apply_ddt(Field3D &f) override;
private:
};

#endif // __BNDRY_STD_H__
