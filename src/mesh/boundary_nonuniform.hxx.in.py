#!/usr/bin/env python3

import sys

from jinja2 import Environment
from boundary_nonuniform_common import orders, boundaries, maybeopen

env = Environment(trim_blocks=True)

header = """
#include <utility>

#include "boundary_op.hxx"

"""
vecs = """
struct vec{{order}}{
{% for i in range(order) %}
  BoutReal f{{i}};
{% endfor %}
  void operator+=(BoutReal v){
{% for i in range(order) %}
    f{{i}} += v;
{% endfor %}
  }
};

"""


class_str = """

class {{class}} : public BoundaryOp {
public:
  {{class}}() {}
  {{class}}(BoundaryRegion *region, std::shared_ptr<FieldGenerator> gen = nullptr) : BoundaryOp(region), gen(std::move(gen)) {}
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
  static vec{{order}} calc_interp_to_stencil(const vec{{order}}& spacing);
{% for stagger in "no", "co", "anti" %}
  void apply_ {{- stagger -}} _stagger (Field3D &f, Mesh * mesh
{{ ", BoutReal t, const std::shared_ptr<FieldGenerator>& fg, std::vector<BoutReal>& vals, int x_boundary_offset, int y_boundary_offset" if boundary != "Free"  }}
);
{% endfor %}
};
"""


if __name__ == "__main__":
    with maybeopen(sys.argv) as print:
        print(header)
        for order in orders:
            print(env.from_string(vecs).render(order=order))
        for order in orders:
            for boundary in boundaries:
                args = {
                    "order": order,
                    "boundary": boundary,
                    "class": "Boundary%sNonUniform_O%d" % (boundary, order),
                }
                print(env.from_string(class_str).render(**args))
