#!/usr/bin/env python3

from jinja2 import Environment

env=Environment(trim_blocks=True);


orders=range(2,5)
whats=["Dirichlet","Neumann","Free"]

header="""
#include <utility>

#include "boundary_op.hxx"
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
  void calc_interp_to_stencil(
{% for i in range(order) %}BoutReal x{{i}}, {% endfor %}
{% for i in range(order) %}BoutReal &fac{{i}}{% if loop.last %}){% else %}, {% endif %}{% endfor %} const ;
};
"""


if __name__ == "__main__":
    print(header)
    for order in orders:
        for what in whats:
            args={
                'order':order,
                'what':what,
                'class':"Boundary%sNonUniform_O%d"%(what,order),
            }
            print(env.from_string(class_str).render(**args))
