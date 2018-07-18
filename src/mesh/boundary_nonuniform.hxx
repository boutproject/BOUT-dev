
class BoundaryDirichletNonUniform_O4 : public BoundaryOp {
public:
  BoundaryDirichletNonUniform_O4() {}
  BoundaryDirichletNonUniform_O4(BoundaryRegion *region, std::shared_ptr<FieldGenerator> gen = nullptr) : BoundaryOp(region), gen(gen) {}
  BoundaryOp *clone(BoundaryRegion *region, const list<string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &f) override {
    throw BoutException("Not Implemented");
  };
 
  void apply(Field3D &f) override {
    apply(f,0.0);
  };
  void apply(Field3D &f, BoutReal t) override;
  
private:
  std::shared_ptr<FieldGenerator>  gen; // Generator
  void calc_interp_to_stencil(BoutReal a, BoutReal b, BoutReal c, BoutReal d, BoutReal &e,
                              BoutReal &f, BoutReal &g, BoutReal &h) const;
};
