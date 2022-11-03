
#include <utility>

#include "boundary_op.hxx"

struct vec2 {
  BoutReal f0;
  BoutReal f1;
  void operator+=(BoutReal v) {
    f0 += v;
    f1 += v;
  }
};

struct vec3 {
  BoutReal f0;
  BoutReal f1;
  BoutReal f2;
  void operator+=(BoutReal v) {
    f0 += v;
    f1 += v;
    f2 += v;
  }
};

struct vec4 {
  BoutReal f0;
  BoutReal f1;
  BoutReal f2;
  BoutReal f3;
  void operator+=(BoutReal v) {
    f0 += v;
    f1 += v;
    f2 += v;
    f3 += v;
  }
};

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
  static vec2 calc_interp_to_stencil(const vec2& spacing);
  void apply_no_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_co_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_anti_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                          std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                          const int x_boundary_offset, const int y_boundary_offset);
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
  static vec2 calc_interp_to_stencil(const vec2& spacing);
  void apply_no_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_co_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_anti_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                          std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                          const int x_boundary_offset, const int y_boundary_offset);
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
  static vec2 calc_interp_to_stencil(const vec2& spacing);
  void apply_no_stagger(Field3D& f, Mesh* mesh

  );
  void apply_co_stagger(Field3D& f, Mesh* mesh

  );
  void apply_anti_stagger(Field3D& f, Mesh* mesh

  );
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
  static vec3 calc_interp_to_stencil(const vec3& spacing);
  void apply_no_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_co_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_anti_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                          std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                          const int x_boundary_offset, const int y_boundary_offset);
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
  static vec3 calc_interp_to_stencil(const vec3& spacing);
  void apply_no_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_co_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_anti_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                          std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                          const int x_boundary_offset, const int y_boundary_offset);
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
  static vec3 calc_interp_to_stencil(const vec3& spacing);
  void apply_no_stagger(Field3D& f, Mesh* mesh

  );
  void apply_co_stagger(Field3D& f, Mesh* mesh

  );
  void apply_anti_stagger(Field3D& f, Mesh* mesh

  );
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
  static vec4 calc_interp_to_stencil(const vec4& spacing);
  void apply_no_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_co_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_anti_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                          std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                          const int x_boundary_offset, const int y_boundary_offset);
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
  static vec4 calc_interp_to_stencil(const vec4& spacing);
  void apply_no_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_co_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                        std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                        const int x_boundary_offset, const int y_boundary_offset);
  void apply_anti_stagger(Field3D& f, Mesh* mesh, BoutReal t,
                          std::shared_ptr<FieldGenerator> fg, std::vector<BoutReal> vals,
                          const int x_boundary_offset, const int y_boundary_offset);
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
  static vec4 calc_interp_to_stencil(const vec4& spacing);
  void apply_no_stagger(Field3D& f, Mesh* mesh

  );
  void apply_co_stagger(Field3D& f, Mesh* mesh

  );
  void apply_anti_stagger(Field3D& f, Mesh* mesh

  );
};
