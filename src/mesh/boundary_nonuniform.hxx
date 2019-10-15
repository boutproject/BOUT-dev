
#include <utility>

#include "boundary_op.hxx"

class BoundaryDirichletNonUniform_O2 : public BoundaryOp {
public:
  BoundaryDirichletNonUniform_O2() {}
  BoundaryDirichletNonUniform_O2(BoundaryRegion* region,
                                 std::shared_ptr<FieldGenerator> gen = nullptr)
      : BoundaryOp(region), gen(gen) {}
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override;

  using BoundaryOp::apply;
  void apply(Field2D& UNUSED(f)) override { throw BoutException("Not Implemented"); };

  void apply(Field3D& f) override { apply(f, 0.0); };
  void apply(Field3D& f, BoutReal t) override;

private:
  std::shared_ptr<FieldGenerator> gen; // Generator
  void calc_interp_to_stencil(BoutReal x0, BoutReal x1, BoutReal& fac0,
                              BoutReal& fac1) const;
};

class BoundaryNeumannNonUniform_O2 : public BoundaryOp {
public:
  BoundaryNeumannNonUniform_O2() {}
  BoundaryNeumannNonUniform_O2(BoundaryRegion* region,
                               std::shared_ptr<FieldGenerator> gen = nullptr)
      : BoundaryOp(region), gen(gen) {}
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override;

  using BoundaryOp::apply;
  void apply(Field2D& UNUSED(f)) override { throw BoutException("Not Implemented"); };

  void apply(Field3D& f) override { apply(f, 0.0); };
  void apply(Field3D& f, BoutReal t) override;

private:
  std::shared_ptr<FieldGenerator> gen; // Generator
  void calc_interp_to_stencil(BoutReal x0, BoutReal x1, BoutReal& fac0,
                              BoutReal& fac1) const;
};

class BoundaryFreeNonUniform_O2 : public BoundaryOp {
public:
  BoundaryFreeNonUniform_O2() {}
  BoundaryFreeNonUniform_O2(BoundaryRegion* region,
                            std::shared_ptr<FieldGenerator> gen = nullptr)
      : BoundaryOp(region), gen(gen) {}
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override;

  using BoundaryOp::apply;
  void apply(Field2D& UNUSED(f)) override { throw BoutException("Not Implemented"); };

  void apply(Field3D& f) override { apply(f, 0.0); };
  void apply(Field3D& f, BoutReal t) override;

private:
  std::shared_ptr<FieldGenerator> gen; // Generator
  void calc_interp_to_stencil(BoutReal x0, BoutReal x1, BoutReal& fac0,
                              BoutReal& fac1) const;
};

class BoundaryDirichletNonUniform_O3 : public BoundaryOp {
public:
  BoundaryDirichletNonUniform_O3() {}
  BoundaryDirichletNonUniform_O3(BoundaryRegion* region,
                                 std::shared_ptr<FieldGenerator> gen = nullptr)
      : BoundaryOp(region), gen(gen) {}
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override;

  using BoundaryOp::apply;
  void apply(Field2D& UNUSED(f)) override { throw BoutException("Not Implemented"); };

  void apply(Field3D& f) override { apply(f, 0.0); };
  void apply(Field3D& f, BoutReal t) override;

private:
  std::shared_ptr<FieldGenerator> gen; // Generator
  void calc_interp_to_stencil(BoutReal x0, BoutReal x1, BoutReal x2, BoutReal& fac0,
                              BoutReal& fac1, BoutReal& fac2) const;
};

class BoundaryNeumannNonUniform_O3 : public BoundaryOp {
public:
  BoundaryNeumannNonUniform_O3() {}
  BoundaryNeumannNonUniform_O3(BoundaryRegion* region,
                               std::shared_ptr<FieldGenerator> gen = nullptr)
      : BoundaryOp(region), gen(gen) {}
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override;

  using BoundaryOp::apply;
  void apply(Field2D& UNUSED(f)) override { throw BoutException("Not Implemented"); };

  void apply(Field3D& f) override { apply(f, 0.0); };
  void apply(Field3D& f, BoutReal t) override;

private:
  std::shared_ptr<FieldGenerator> gen; // Generator
  void calc_interp_to_stencil(BoutReal x0, BoutReal x1, BoutReal x2, BoutReal& fac0,
                              BoutReal& fac1, BoutReal& fac2) const;
};

class BoundaryFreeNonUniform_O3 : public BoundaryOp {
public:
  BoundaryFreeNonUniform_O3() {}
  BoundaryFreeNonUniform_O3(BoundaryRegion* region,
                            std::shared_ptr<FieldGenerator> gen = nullptr)
      : BoundaryOp(region), gen(gen) {}
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override;

  using BoundaryOp::apply;
  void apply(Field2D& UNUSED(f)) override { throw BoutException("Not Implemented"); };

  void apply(Field3D& f) override { apply(f, 0.0); };
  void apply(Field3D& f, BoutReal t) override;

private:
  std::shared_ptr<FieldGenerator> gen; // Generator
  void calc_interp_to_stencil(BoutReal x0, BoutReal x1, BoutReal x2, BoutReal& fac0,
                              BoutReal& fac1, BoutReal& fac2) const;
};

class BoundaryDirichletNonUniform_O4 : public BoundaryOp {
public:
  BoundaryDirichletNonUniform_O4() {}
  BoundaryDirichletNonUniform_O4(BoundaryRegion* region,
                                 std::shared_ptr<FieldGenerator> gen = nullptr)
      : BoundaryOp(region), gen(gen) {}
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override;

  using BoundaryOp::apply;
  void apply(Field2D& UNUSED(f)) override { throw BoutException("Not Implemented"); };

  void apply(Field3D& f) override { apply(f, 0.0); };
  void apply(Field3D& f, BoutReal t) override;

private:
  std::shared_ptr<FieldGenerator> gen; // Generator
  void calc_interp_to_stencil(BoutReal x0, BoutReal x1, BoutReal x2, BoutReal x3,
                              BoutReal& fac0, BoutReal& fac1, BoutReal& fac2,
                              BoutReal& fac3) const;
};

class BoundaryNeumannNonUniform_O4 : public BoundaryOp {
public:
  BoundaryNeumannNonUniform_O4() {}
  BoundaryNeumannNonUniform_O4(BoundaryRegion* region,
                               std::shared_ptr<FieldGenerator> gen = nullptr)
      : BoundaryOp(region), gen(gen) {}
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override;

  using BoundaryOp::apply;
  void apply(Field2D& UNUSED(f)) override { throw BoutException("Not Implemented"); };

  void apply(Field3D& f) override { apply(f, 0.0); };
  void apply(Field3D& f, BoutReal t) override;

private:
  std::shared_ptr<FieldGenerator> gen; // Generator
  void calc_interp_to_stencil(BoutReal x0, BoutReal x1, BoutReal x2, BoutReal x3,
                              BoutReal& fac0, BoutReal& fac1, BoutReal& fac2,
                              BoutReal& fac3) const;
};

class BoundaryFreeNonUniform_O4 : public BoundaryOp {
public:
  BoundaryFreeNonUniform_O4() {}
  BoundaryFreeNonUniform_O4(BoundaryRegion* region,
                            std::shared_ptr<FieldGenerator> gen = nullptr)
      : BoundaryOp(region), gen(gen) {}
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override;

  using BoundaryOp::apply;
  void apply(Field2D& UNUSED(f)) override { throw BoutException("Not Implemented"); };

  void apply(Field3D& f) override { apply(f, 0.0); };
  void apply(Field3D& f, BoutReal t) override;

private:
  std::shared_ptr<FieldGenerator> gen; // Generator
  void calc_interp_to_stencil(BoutReal x0, BoutReal x1, BoutReal x2, BoutReal x3,
                              BoutReal& fac0, BoutReal& fac1, BoutReal& fac2,
                              BoutReal& fac3) const;
};
