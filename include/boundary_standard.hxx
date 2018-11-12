/// Some standard boundary conditions

#ifndef __BNDRY_STD_H__
#define __BNDRY_STD_H__

#include "boundary_op.hxx"
#include "bout_types.hxx"
#include "unused.hxx"
#include <field_factory.hxx>

#include <utility>

/// Dirichlet (set to zero) boundary condition
class BoundaryDirichlet : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  void apply(Field2D &f, BoutReal t = 0.) final { applyTemplate<Field2D>(f, t); }
  void apply(Field3D &f, BoutReal t = 0.) final { applyTemplate<Field3D>(f, t); }

private:
  template <typename T> void applyTemplate(T &f, BoutReal t);
};

BoutReal default_func(BoutReal t, int x, int y, int z);

/// 3nd-order boundary condition
class BoundaryDirichlet_O3 : public BoundaryOpWithApply<BoundaryDirichlet_O3> {
public:
  using BoundaryOpWithApply<
      BoundaryDirichlet_O3>::BoundaryOpWithApply; // inherit BoundaryOp constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// 4th-order boundary condition
class BoundaryDirichlet_O4 : public BoundaryOpWithApply<BoundaryDirichlet_O4> {
public:
  using BoundaryOpWithApply<BoundaryDirichlet_O4>::
      BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// Dirichlet boundary condition, tries to smooth out grid-scale oscillations at the
/// boundary
class BoundaryDirichlet_smooth : public BoundaryOpWithApply<BoundaryDirichlet_smooth> {
public:
  using BoundaryOpWithApply<BoundaryDirichlet_smooth>::
      BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// Dirichlet boundary condition set half way between guard cell and grid cell at 2nd
/// order accuracy
class BoundaryDirichlet_2ndOrder
    : public BoundaryOpWithApply<BoundaryDirichlet_2ndOrder> {
public:
  using BoundaryOpWithApply<BoundaryDirichlet_2ndOrder>::
      BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// Dirichlet boundary condition set half way between guard cell and grid cell at 4th
/// order accuracy
class BoundaryDirichlet_O5 : public BoundaryOpWithApply<BoundaryDirichlet_O5> {
public:
  using BoundaryOpWithApply<BoundaryDirichlet_O5>::
      BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// Neumann (zero-gradient) boundary condition for non-orthogonal meshes
class BoundaryNeumann_NonOrthogonal : public BoundaryOp {
public:
  BoundaryNeumann_NonOrthogonal() {}
  BoundaryNeumann_NonOrthogonal(BoutReal setval)
      : BoundaryOp(nullptr, setval, nullptr) {}
  BoundaryNeumann_NonOrthogonal(BoundaryRegion *region, BoutReal setval = 0.)
      : BoundaryOp(region, setval, nullptr) {}
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  void apply(Field2D &f, BoutReal t = 0.) final { applyTemplate(f, t); }
  void apply(Field3D &f, BoutReal t = 0.) final { applyTemplate(f, t); }

private:
  template <typename T> void applyTemplate(T &f, BoutReal t);
};

/// Neumann (zero-gradient) boundary condition, using 2nd order on boundary
class BoundaryNeumann2 : public BoundaryOpWithApply<BoundaryNeumann2> {
public:
  using BoundaryOpWithApply<
      BoundaryNeumann2>::BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// Neumann boundary condition set half way between guard cell and grid cell at 2nd order
/// accuracy
class BoundaryNeumann_2ndOrder
    : public BoundaryOpWithApply<BoundaryNeumann_2ndOrder, true> {
public:
  using BoundaryOpWithApply<BoundaryNeumann_2ndOrder, true>::
      BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

// Neumann boundary condition set half way between guard cell and grid cell at 2nd order
// accuracy
class BoundaryNeumann : public BoundaryOpWithApply<BoundaryNeumann, true> {
public:
  using BoundaryOpWithApply<BoundaryNeumann, true>::
      BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// Neumann boundary condition set half way between guard cell and grid cell at 4th order
/// accuracy
class BoundaryNeumann_O4 : public BoundaryOpWithApply<BoundaryNeumann_O4, true> {
public:
  using BoundaryOpWithApply<BoundaryNeumann_O4, true>::
      BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// Neumann boundary condition set half way between guard cell and grid cell at 4th order
/// accuracy
class BoundaryNeumann_4thOrder
    : public BoundaryOpWithApply<BoundaryNeumann_4thOrder, true> {
public:
  using BoundaryOpWithApply<BoundaryNeumann_4thOrder, true>::
      BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// NeumannPar (zero-gradient) boundary condition on
/// the variable / sqrt(g_22)
class BoundaryNeumannPar : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOp constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  void apply(Field2D &f, BoutReal t = 0.) final { applyTemplate(f, t); }
  void apply(Field3D &f, BoutReal t = 0.) final { applyTemplate(f, t); }

private:
  template <typename T> void applyTemplate(T &f, BoutReal t);
};

/// Robin (mix of Dirichlet and Neumann)
class BoundaryRobin : public BoundaryOp {
public:
  BoundaryRobin() : aval(0.), bval(0.), gval(0.) {}
  BoundaryRobin(BoundaryRegion *region, BoutReal a, BoutReal b, BoutReal g)
      : BoundaryOp(region), aval(a), bval(b), gval(g) {}
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  void apply(Field2D &f, BoutReal t = 0.) final { applyTemplate(f, t); }
  void apply(Field3D &f, BoutReal t = 0.) final { applyTemplate(f, t); }

private:
  const BoutReal aval, bval, gval;

  template <typename T> void applyTemplate(T &f, BoutReal t);
};

/// Constant gradient (zero second derivative)
class BoundaryConstGradient : public BoundaryOpWithApply<BoundaryConstGradient> {
public:
  using BoundaryOpWithApply<BoundaryConstGradient>::
      BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/// Zero Laplacian, decaying solution
class BoundaryZeroLaplace : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp;
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  void apply(Field2D &f, BoutReal t) final;
  void apply(Field3D &f, BoutReal t) final;
};

/// Zero Laplacian
class BoundaryZeroLaplace2 : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp;
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  void apply(Field2D &f, BoutReal t) final;
  void apply(Field3D &f, BoutReal t) final;
};

/// Constant Laplacian, decaying solution
class BoundaryConstLaplace : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp;
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  void apply(Field2D &f, BoutReal t) final;
  void apply(Field3D &f, BoutReal t) final;
};

/// Vector boundary condition Div(B) = 0, Curl(B) = 0
class BoundaryDivCurl : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp;
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  void apply(Field2D &UNUSED(f), BoutReal UNUSED(t)) final {
    throw BoutException("ERROR: DivCurl boundary only for vectors");
  }
  void apply(Field3D &UNUSED(f), BoutReal UNUSED(t)) final {
    throw BoutException("ERROR: DivCurl boundary only for vectors");
  }
  void apply(Vector2D &f) final;
  void apply(Vector3D &f) final;
};

/// Free boundary condition (evolve the field in the guard cells, using non-centred
/// derivatives to calculate the ddt)
class BoundaryFree : public BoundaryOp {
public:
  using BoundaryOp::BoundaryOp; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  void apply(Field2D &f, BoutReal t) final;
  void apply(Field3D &f, BoutReal t) final;

  void apply_ddt(Field2D &f) final;
  void apply_ddt(Field3D &f) final;
};

// L. Easy
/// Alternative free boundary condition (evolve the field in the guard cells, using
/// non-centred derivatives to calculate the ddt)
class BoundaryFree_O2 : public BoundaryOpWithApply<BoundaryFree_O2> {
public:
  using BoundaryOpWithApply<
      BoundaryFree_O2>::BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

class BoundaryFree_O3 : public BoundaryOpWithApply<BoundaryFree_O3> {
public:
  using BoundaryOpWithApply<
      BoundaryFree_O3>::BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};
// End L.Easy

class BoundaryFree_O4 : public BoundaryOpWithApply<BoundaryFree_O4> {
public:
  using BoundaryOpWithApply<
      BoundaryFree_O4>::BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

class BoundaryFree_O5 : public BoundaryOpWithApply<BoundaryFree_O5> {
public:
  using BoundaryOpWithApply<
      BoundaryFree_O5>::BoundaryOpWithApply; // inherit BoundaryOpWithApply constructors
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args,
                    const std::map<std::string, std::string> &keywords) override;

  static void applyAtPoint(Field2D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPoint(Field3D &f, BoutReal val, int x, int bx, int y, int by, int z,
                           BoutReal delta);
  static void applyAtPointStaggered(Field2D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void applyAtPointStaggered(Field3D &f, BoutReal val, int x, int bx, int y,
                                    int by, int z, BoutReal delta);
  static void extrapolateFurther(Field2D &f, int x, int bx, int y, int by, int z);
  static void extrapolateFurther(Field3D &f, int x, int bx, int y, int by, int z);
};

/////////////////////////////////////////////////////////

/// Convert a boundary condition to a relaxing one
class BoundaryRelax : public BoundaryModifier {
public:
  BoundaryRelax() : BoundaryModifier(true), r(10.) {} // Set default rate
  BoundaryRelax(BoundaryOp *operation, BoutReal rate)
      : BoundaryModifier(operation, true) {
    r = fabs(rate);
  }
  BoundaryOp *cloneMod(BoundaryOp *op, const list<string> &args) final;

  void apply(Field2D &f, BoutReal t) final;
  void apply(Field3D &f, BoutReal t) final;

  void apply_ddt(Field2D &f) final;
  void apply_ddt(Field3D &f) final;

private:
  BoutReal r;
};

/// Increase the width of a boundary
class BoundaryWidth : public BoundaryModifier {
public:
  BoundaryWidth() : bndry(nullptr) {}
  BoundaryWidth(BoundaryOp *operation, int wid);
  BoundaryOp *cloneMod(BoundaryOp *op,
                       const list<string> &args) final;

  void apply(Field2D &f, BoutReal t) final {
    op->apply(f, t);
  }
  void apply(Field3D &f, BoutReal t) final {
    op->apply(f, t);
  }

  void apply_ddt(Field2D &f) final {
    op->apply_ddt(f);
  }
  void apply_ddt(Field3D &f) final {
    op->apply_ddt(f);
  }

private:
  std::unique_ptr<BoundaryRegion> bndry;
};

/// Convert input field fromFieldAligned, apply boundary and then convert back
/// toFieldAligned
/// Equivalent to converting the boundary condition to "Field Aligned" from "orthogonal"
class BoundaryToFieldAligned : public BoundaryModifier {
public:
  BoundaryToFieldAligned() {}
  BoundaryToFieldAligned(BoundaryOp *operation) : BoundaryModifier(operation) {}
  BoundaryOp *cloneMod(BoundaryOp *op, const list<string> &args) final;

  void apply(Field2D &f, BoutReal t) final;
  void apply(Field3D &f, BoutReal t) final;

  void apply_ddt(Field2D &f) final;
  void apply_ddt(Field3D &f) final;

private:
};

/// Convert input field toFieldAligned, apply boundary and then convert back
/// fromFieldAligned
/// Equivalent to converting the boundary condition from "Field Aligned" to "orthogonal"
class BoundaryFromFieldAligned : public BoundaryModifier {
public:
  BoundaryFromFieldAligned() {}
  BoundaryFromFieldAligned(BoundaryOp *operation) : BoundaryModifier(operation) {}
  BoundaryOp *cloneMod(BoundaryOp *op, const list<string> &args) final;

  void apply(Field2D &f, BoutReal t) final;
  void apply(Field3D &f, BoutReal t) final;

  void apply_ddt(Field2D &f) final;
  void apply_ddt(Field3D &f) final;

private:
};

#endif // __BNDRY_STD_H__
