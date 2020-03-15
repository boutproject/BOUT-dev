#!/usr/bin/env python3

from jinja2 import Environment
import stencils_sympy as sten

header="""\
#include <boundary_standard.hxx>
#include <bout/constants.hxx>
#include <boutexception.hxx>
#include <derivs.hxx>
#include <fft.hxx>
#include <globals.hxx>
#include <invert_laplace.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <utils.hxx>

#include "boundary_nonuniform.hxx"

static void update_stagger_offsets(int& x_boundary_offset, int& y_boundary_offset, int& stagger, CELL_LOC loc){
  // NB: bx is going outwards
  // NB: XLOW means shifted in -x direction
  // `stagger` stagger direction with respect to direction of boundary
  //   0 : no stagger or orthogonal to boundary direction
  //   1 : staggerd in direction of boundary
  //  -1 : staggerd in oposite direction of boundary
  // Also note that all offsets are basically half a cell
  if (loc == CELL_XLOW) {
    if (x_boundary_offset == 0) {
      x_boundary_offset = -1;
    } else if (x_boundary_offset < 0) {
      stagger = -1;
    } else {
      stagger = 1;
    }
  }
  if (loc == CELL_YLOW) {
    if (y_boundary_offset == 0) {
      y_boundary_offset = -1;
    } else if (y_boundary_offset < 0) {
      stagger = -1;
    } else {
      stagger = 1;
    }
  }
}
"""

env=Environment(trim_blocks=True);

apply_str="""
void Boundary{{type}}NonUniform_O{{order}}::apply(Field3D &f, MAYBE_UNUSED(BoutReal t)) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  Mesh *mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

{% if type != "Free" %}
  std::vector<BoutReal> vals;
  vals.reserve(mesh->LocalNz);
{% endif %}

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
{% if type != "Free" %}
    if (fg) {
      // Calculate the X and Y normalised values half-way between the guard cell and
      // grid cell
      const BoutReal xnorm = 0.5 *
                     (mesh->GlobalX(bndry->x)          // In the guard cell
                      + mesh->GlobalX(bndry->x - x_boundary_offset)); // the grid cell
      const BoutReal ynorm = TWOPI * 0.5 *
                     (mesh->GlobalY(bndry->y)          // In the guard cell
                      + mesh->GlobalY(bndry->y - y_boundary_offset)); // the grid cell
      const BoutReal zfac =  TWOPI / mesh->LocalNz;
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        vals[zk] = fg->generate(bout::generator::Context().set("x", xnorm, "y", ynorm, "z", zfac * zk, "t" , t));
      }
    }
{% endif %}{# type != Free #}


    vec{{order}} spacing;
    vec{{order}} facs;

    const Field2D &coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
{% if type != "Free" %}
{% for i in range(1,order) %}
    Indices i{{i}}{bndry->x - {{i}} * bndry->bx, bndry->y - {{i}} * bndry->by, 0};
{% endfor %}
    if (stagger == 0) {
      BoutReal t;
      spacing.f0 = 0;
      BoutReal st=0;
{% for i in range(1,order) %}
      t = coords_field(i{{i}}.x, i{{i}}.y);
      spacing.f{{i}} = st + t / 2;
      st += t;
{% endfor %}
    } else {
      spacing.f0 = 0;
{% if type == "Neumann" %}
      // Check if we are staggered and also boundary in low
      //  direction
      // In the case of Neumann we have in this case two values
      //  defined at the same point
      if (stagger == -1 && 
              (    (bndry->bx && x_boundary_offset == -1)
                || (bndry->by && y_boundary_offset == -1))){
        spacing.f1 = spacing.f0;
{% for i in range(2,order) %}
        spacing.f{{i}} = spacing.f{{i-1}} + coords_field(i{{i-1}}.x, i{{i-1}}.y);
{% endfor %}
      } else {
{% for i in range(1,order) %}
        spacing.f{{i}} = spacing.f{{i-1}} + coords_field(i{{i}}.x, i{{i}}.y);
{% endfor %}
      }
{% else %}
{% for i in range(1,order) %}
        spacing.f{{i}} = spacing.f{{i-1}} + coords_field(i{{i}}.x, i{{i}}.y);
{% endfor %}
{% endif %} {# end of special case for neuman #}

    }
{% if type == "Dirichlet" %}
    if (stagger == -1) {
{% for i in range(1,order) %}
      i{{i}} = {bndry->x - {{i+1}} * bndry->bx, bndry->y - {{i+1}} * bndry->by, 0};
{% endfor %}
    }
{% endif %}
{% else %} {# type != Free #}
{% for i in range(order) %}
    const Indices i{{i}}{bndry->x - {{i+1}} * bndry->bx, bndry->y - {{i+1}} * bndry->by, 0};
{% endfor %}
    if (stagger == 0) {
      BoutReal st=0;
      BoutReal t;
{% for i in range(order) %}
      t = coords_field(i{{i}}.x, i{{i}}.y);
      spacing.f{{i}} = st + t / 2;
      st += t;
{% endfor %}
    } else {
      spacing.f0 = coords_field(i0.x, i0.y);
{% for i in range(1,order) %}
      spacing.f{{i}} = spacing.f{{i-1}} + coords_field(i{{i}}.x, i{{i}}.y);
{% endfor %}
    }

{% endif %} {# type != Free #}
{% if type == "Dirichlet" %}
    // with dirichlet, we specify the value on the boundary, even if
    // the value is part of the evolving system.
    for (int i = ((stagger == -1) ? -1 : 0); i < bndry->width; i++) {
{% else %}
    // With free and neumann the value is not set if the point is
    // evolved and it is on the boundary.
    for (int i = 0; i < bndry->width; i++) {
{% endif %}
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        BoutReal to_add = coords_field(ic.x, ic.y) / 2;
        spacing += to_add;
        facs = calc_interp_to_stencil(spacing);
        spacing += to_add;
      } else {
        BoutReal to_add = coords_field(ic.x, ic.y);
{% if type == "Free" %}
{% elif type == "Dirichlet" %}
        if (stagger == -1
              && i != -1) {
          spacing += to_add;
        }
{% else %}
        if (stagger == -1) {
          spacing += to_add;
        }
{% endif %}
        facs = calc_interp_to_stencil(spacing);
{% if type != "Free" %}
        if (stagger == 1) {
          spacing += to_add;
        }
{% else %}
        spacing += to_add;
{% endif %}
      }
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
{% if type != "Free" %}
        const BoutReal val = (fg) ? vals[iz] : 0.0;
        const BoutReal set = facs.f0 * val
{% else %}
        const BoutReal set = facs.f0 * f(i0.x, i0.y, iz)
{% endif %}
{% for i in range(1,order) %}
           + facs.f{{i}} *f(i{{i}}.x, i{{i}}.y, iz)
{% endfor %}
        ;
        
        f(ic.x, ic.y, iz) = set;
      }
    }
  }
}
"""
clone_str="""
BoundaryOp * {{class}}::clone(BoundaryRegion *region,
   const std::list<std::string> &args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new {{class}}(region, newgen);
}
"""
stencil_str="""
vec{{order}} {{class}}::calc_interp_to_stencil(const vec{{order}} & spacing) const {
vec{{order}} facs;
// Stencil Code
{{stencil_code}}
return facs;
}
"""

orders=range(2,5)
boundaries=["Dirichlet","Neumann","Free"]

if __name__ == "__main__":
    print(header)

    for order in orders:
        for boundary in boundaries:
            if boundary == "Neumann":
                mat=sten.neumann
            else:
                mat=sten.dirichlet
            try:
                code=sten.gen_code(order,mat)
            except:
                import sys
                print("Order:",order,"boundary:",boundary,file=sys.stderr)
                raise
            args={
                'order':order,
                'type':boundary,
                'class':"Boundary%sNonUniform_O%d"%(boundary,order),
                'stencil_code': code,
            }

            print(env.from_string(apply_str).render(**args))
            print(env.from_string(clone_str).render(**args))
            print(env.from_string(stencil_str).render(**args))
