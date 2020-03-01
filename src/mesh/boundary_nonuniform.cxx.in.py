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
"""

env=Environment(trim_blocks=True);

apply_str="""
void Boundary{{type}}NonUniform_O{{order}}::apply(Field3D &f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh *mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  // NB: bx is going outwards
  int bx = bndry->bx;
  int by = bndry->by;
  // NB: XLOW means shifted in -x direction
  // `stagger` stagger direction with respect to direction of boundary
  //   0 : no stagger or orthogonal to boundary direction
  //   1 : staggerd in direction of boundary
  //  -1 : staggerd in oposite direction of boundary
  // Also note that all offsets are basically half a cell
  int stagger = 0;
  if (loc == CELL_XLOW) {
    if (bx == 0) {
      bx = -1;
    } else if (bx < 0) {
      stagger = -1;
    } else {
      stagger = 1;
    }
  }
  if (loc == CELL_YLOW) {
    if (by == 0) {
      by = -1;
    } else if (by < 0) {
      stagger = -1;
    } else {
      stagger = 1;
    }
  }
{% if type == "Dirichlet" %}
  int istart = (stagger == -1) ? -1 : 0;
{% endif %}

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5 * (mesh->GlobalX(bndry->x)          // In the guard cell
                                + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5 * (mesh->GlobalY(bndry->y)          // In the guard cell
                                + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }
{% for i in range(order) %}
    BoutReal x{{i}};
{% endfor %}
    fac{{order}} facs;

    const Field2D &spacing =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
{% if type != "Free" %}
{% for i in range(1,order) %}
    Indices i{{i}}{bndry->x - {{i}} * bndry->bx, bndry->y - {{i}} * bndry->by, 0};
{% endfor %}
    BoutReal t;
    if (stagger == 0) {
      x0 = 0;
      BoutReal st=0;
{% for i in range(1,order) %}
      t = spacing(i{{i}}.x, i{{i}}.y);
      x{{i}} = st + t / 2;
      st += t;
{% endfor %}
    } else {
      x0 = 0; // spacing(bndry->x, bndry->y) / 2;
{% for i in range(1,order) %}
      x{{i}} = x{{i-1}} + spacing(i{{i}}.x, i{{i}}.y);
{% endfor %}
    }
{% if type == "Dirichlet" %}
    if (stagger == -1) {
{% for i in range(1,order) %}
      i{{i}} = {bndry->x - {{i+1}} * bndry->bx, bndry->y - {{i+1}} * bndry->by, 0};
{% endfor %}
    }
{% endif %}
{% else %}
{% for i in range(order) %}
    Indices i{{i}}{bndry->x - {{i+1}} * bndry->bx, bndry->y - {{i+1}} * bndry->by, 0};
{% endfor %}
    if (stagger == 0) {
      BoutReal st=0;
{% for i in range(order) %}
      t = spacing(i{{i}}.x, i{{i}}.y);
      x{{i}} = st + t / 2;
      st += t;
{% endfor %}
    } else {
      x0 = 0;
{% for i in range(1,order) %}
      x{{i}} = x{{i-1}} + spacing(i{{i}}.x, i{{i}}.y);
{% endfor %}
    }

{% endif %}
{% if type == "Dirichlet" %}
    for (int i = istart; i < bndry->width; i++) {
{% else %}
    for (int i = 0; i < bndry->width; i++) {
{% endif %}
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = spacing(ic.x, ic.y) / 2;
{% for i in range(order) %}
        x{{i}} += t;
{% endfor %}
        // printf("%+2d: %d %d %g %g %g %g\\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        facs = calc_interp_to_stencil(
{% for i in range(order) %}x{{i}}{% if not loop.last %}, {% endif %}{% endfor %});
{% for i in range(order) %}
        x{{i}} += t;
{% endfor %}
      } else {
        t = spacing(ic.x, ic.y);
        if (stagger == -1
{% if type == "Dirichlet" %}
              && i != -1
{% endif %}
                     ) {
{% for i in range(order) %}
          x{{i}} += t;
{% endfor %}
        }
        facs = calc_interp_to_stencil(
{% for i in range(order) %}x{{i}}{% if not loop.last %}, {% endif %}{% endfor %});
        if (stagger == 1) {
{% for i in range(order) %}
          x{{i}} += t;
{% endfor %}
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
{% for i in range(1,order) %}
        i{{i}}.z = ic.z;
{% endfor %}
{% if type != "Free" %}
        val = (fg) ? vals[ic.z] : 0.0;
        t = facs.f0 * val 
{% else %}
        t = facs.f0 * f(i0.x, i0.y, i0.z)
{% endif %}
{% for i in range(1,order) %}
           + facs.f{{i}} *f(i{{i}}.x, i{{i}}.y, i{{i}}.z)
{% endfor %}
        ;
        
        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}
"""
clone_str="""
BoundaryOp * {{class}}::clone(BoundaryRegion *region,
   const std::list<std::string> &args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new {{class}}(region, newgen);
}
"""
stencil_str="""
fac{{order}} {{class}}::calc_interp_to_stencil(
{% for i in range(order) %}BoutReal x{{i}}{% if loop.last %}){% else %}, {% endif %}{% endfor %} const {
fac{{order}} facs;
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
