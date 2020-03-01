#!/usr/bin/env python3

from jinja2 import Environment

env=Environment(trim_blocks=True);


orders=range(2,5)
boundaries=["Dirichlet","Neumann","Free"]

header="""
#include <utility>

#include "boundary_op.hxx"

struct fac2{
BoutReal f0, f1;};
struct fac3{
BoutReal f0, f1, f2;};
struct fac4{
BoutReal f0;
BoutReal f1;
BoutReal f2;
BoutReal f3;};
"""

class_str="""

class {{class}} : public BoundaryOp {
public:
  {{class}}() {}
  {{class}}(BoundaryRegion *region, std::shared_ptr<FieldGenerator> gen = nullptr) : BoundaryOp(region), gen(gen) {}
  BoundaryOp *clone(BoundaryRegion *region, const std::list<std::string> &args) override;

  using BoundaryOp::apply;
  void apply(Field2D &UNUSED(f)) override {
    throw BoutException("Not Implemented");
  };
 
  void apply(Field3D &f) override {
    apply(f,0.0);
  };
  void apply(Field3D &f, BoutReal t) override;
  
private:
  std::shared_ptr<FieldGenerator>  gen; // Generator
  fac{{order}} calc_interp_to_stencil(
{% for i in range(order) %}BoutReal x{{i}}{% if loop.last %}){% else %}, {% endif %}{% endfor %} const ;
};
"""


if __name__ == "__main__":
    print(header)
    for order in orders:
        for boundary in boundaries:
            args={
                'order':order,
                'boundary':boundary,
                'class':"Boundary%sNonUniform_O%d"%(boundary,order),
            }
            print(env.from_string(class_str).render(**args))
