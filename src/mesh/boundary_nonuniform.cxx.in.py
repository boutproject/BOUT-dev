#!/usr/bin/env python3

from jinja2 import Environment
import sys
import stencils_sympy as sten

from boundary_nonuniform_common import orders, boundaries, maybeopen

header = """\
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

#if ! BOUT_USE_METRIC_3D
using IndMetric = Ind2D;
#define IND(var, x, y, z) IndMetric var(x*localNy+y, localNy, 1)
#define IND3D(var, x, y, z) Ind3D var((x*localNy+y)*localNz + z, localNy, localNz)
#else
using IndMetric = Ind3D;
#define IND(var, x, y, z) IndMetric var((x*localNy+y)*localNz + z, localNy, localNz)
#define IND3D(var, x, y, z)
#endif

"""

env = Environment(trim_blocks=True)

apply_str = """
{% macro setindicies(order, start=1, offset=0, const='') -%}
    const int localNy = mesh->LocalNy;
    const int localNz = mesh->LocalNz;
    const int offset = (bndry->bx * localNy + bndry->by)
# if BOUT_USE_METRIC_3D
         * localNz
# endif
                ;
    const IND(temp, bndry->x, bndry->y, 0);
    IND3D(temp3d,  bndry->x, bndry->y, 0);
{% for i in range(start, order) %}
    {{ const }} IndMetric i{{i}}{temp - {{i + offset}} * offset};
# if ! BOUT_USE_METRIC_3D
    {{ const }} Ind3D i{{i}}3d{temp3d - {{i + offset}} * offset * localNz};
# endif
{% endfor %}
{%- endmacro %}
void Boundary{{type}}NonUniform_O{{order}}::apply(Field3D &f, MAYBE_UNUSED(BoutReal t)) {
  bndry->first();
  Mesh *mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

{% if with_fg %}
  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg) {
    fg = f.getBndryGenerator(bndry->location);
  }

  std::vector<BoutReal> vals;
  vals.reserve(mesh->LocalNz);
{% endif %}

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
{% if with_fg %}
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
{% endif %}

    vec{{order}} spacing;

    const Coordinates::FieldMetric &coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;

{% if type == "Dirichlet" %}{# DIRICHLET  DIRICHLET  DIRICHLET  DIRICHLET  DIRICHLET  DIRICHLET #}
# if BOUT_USE_METRIC_3D
    for (int iz{0}; iz < mesh->LocalNz; iz++) {
# else
    const int iz = 0;
# endif
{{    setindicies(order) }}
      if (stagger == 0) {
        spacing.f0 = 0;
        BoutReal total_offset = 0;
{% for i in range(1,order) %}
        {% if loop.first %}BoutReal {% endif %}offset = coords_field[i{{i}}+iz];
        spacing.f{{i}} = total_offset + offset / 2;
 {% if not loop.last %}
        total_offset += offset;
 {% endif %}
{% endfor %}
      } else {
        spacing.f0 = 0;
{% for i in range(1,order) %}
          spacing.f{{i}} = spacing.f{{i-1}} + coords_field[i{{i}} + iz];
{% endfor %}
      }
      if (stagger == -1) {
{% for i in range(1,order) %}
# if BOUT_USE_METRIC_3D
        i{{i}} = temp - {{i+1}} * offset;
# else
        i{{i}}3d = temp3d - {{i+1}} * offset;
# endif
{% endfor %}
      }
      
      // with dirichlet, we specify the value on the boundary, even if
      // the value is part of the evolving system.
      for (int i = ((stagger == -1) ? -1 : 0); i < bndry->width; i++) {
        IndMetric icm{temp + IndMetric(i*offset)};
# if BOUT_USE_METRIC_3D
        icm + iz;
        IndMetric ic = icm;
# else
        Ind3D ic{temp3d + Ind3D(i*offset*localNz)};
# endif
        vec{{order}} facs;
        if (stagger == 0) {
          BoutReal to_add = coords_field[icm] / 2;
          spacing += to_add;
          facs = calc_interp_to_stencil(spacing);
          spacing += to_add;
        } else {
          if (stagger == -1
                && i != -1) {
            spacing += coords_field[icm];
          }
          facs = calc_interp_to_stencil(spacing);
          if (stagger == 1) {
            spacing += coords_field[icm];
          }
        }
{% elif type == "Neumann" %}{# NEUMANN  NEUMANN  NEUMANN  NEUMANN  NEUMANN  NEUMANN  NEUMANN #}
{{  setindicies(order) }}
# if BOUT_USE_METRIC_3D
    for (int iz{0}; iz < mesh->LocalNz; iz++) {
# else
    const int iz = 0;
# endif
      if (stagger == 0) {
        spacing.f0 = 0;
        BoutReal total_offset=0;
{% for i in range(1,order) %}
        {% if loop.first %}BoutReal {% endif %}offset = coords_field[i{{i}} + iz];
        spacing.f{{i}} = total_offset + offset / 2;
 {% if not loop.last %}
        total_offset += offset;
 {% endif %}
{% endfor %}
      } else { // stagger != 0
        spacing.f0 = 0;
        // Check if we are staggered and also boundary in low
        //  direction
        // In the case of Neumann we have in this case two values
        //  defined at the same point
        if (stagger == -1 &&
                (    (bndry->bx && x_boundary_offset == -1)
                  || (bndry->by && y_boundary_offset == -1))){
          spacing.f1 = spacing.f0;
{% for i in range(2,order) %}
          spacing.f{{i}} = spacing.f{{i-1}} + coords_field[i{{i-1}} + iz];
{% endfor %}
        } else {
{% for i in range(1,order) %}
          spacing.f{{i}} = spacing.f{{i-1}} + coords_field[i{{i}} + iz];
{% endfor %}
        }
      } // stagger != 0
      // With neumann (and free) the value is not set if the point is
      // evolved and it is on the boundary.
      for (int i = 0; i < bndry->width; i++) {
        IndMetric icm{temp + IndMetric(i*offset)};
# if BOUT_USE_METRIC_3D
        icm += iz;
        IndMetric ic = icm;
# else
        Ind3D ic{temp3d + Ind3D(i*offset*localNz)};
# endif
        vec{{order}} facs;
        if (stagger == 0) {
          BoutReal to_add = coords_field[icm] / 2;
          spacing += to_add;
          facs = calc_interp_to_stencil(spacing);
          spacing += to_add;
        } else {
          if (stagger == -1) {
            spacing += coords_field[icm];
          }
          facs = calc_interp_to_stencil(spacing);
          if (stagger == 1) {
            spacing += coords_field[icm];
          }
        }
{% elif type == "Free" %}{# FREE  FREE  FREE  FREE  FREE  FREE  FREE  FREE  FREE  FREE  FREE  FREE #}
{{  setindicies(order, 0, 1, 'const ') }}
# if BOUT_USE_METRIC_3D
    for (int iz{0}; iz < mesh->LocalNz; iz++) {
# else
    const int iz = 0;
# endif
      if (stagger == 0) {
        BoutReal total_offset = 0;
{% for i in range(order) %}
        {% if loop.first %}BoutReal {% endif %}offset = coords_field[i{{i}} + iz];
        spacing.f{{i}} = total_offset + offset / 2;
 {% if not loop.last %}
        total_offset += offset;
 {% endif %}
{% endfor %}
      } else {
        spacing.f0 = coords_field[i0 + iz];
{% for i in range(1,order) %}
        spacing.f{{i}} = spacing.f{{i-1}} + coords_field[i{{i}} + iz];
{% endfor %}
      }

      // With free (and neumann) the value is not set if the point is
      // evolved and it is on the boundary.
      for (int i = 0; i < bndry->width; i++) {
        IndMetric icm{temp + IndMetric(i*offset)};
# if BOUT_USE_METRIC_3D
        icm += iz;
        IndMetric ic = icm;
# else
        Ind3D ic{temp3d + Ind3D(i*offset*localNz)};
# endif
        vec{{order}} facs;
        if (stagger == 0) {
          BoutReal to_add = coords_field[icm] / 2;
          spacing += to_add;
          facs = calc_interp_to_stencil(spacing);
          spacing += to_add;
        } else {
          facs = calc_interp_to_stencil(spacing);
          spacing += coords_field[icm];
        }
{% endif %}{# END  END  END  END  END  END  END  END  END  END  END  END #}
# if ! BOUT_USE_METRIC_3D
    for (int iz{0}; iz < mesh->LocalNz; iz++) {
# endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
        const BoutReal val =
{% if type != "Free" %}
           (fg) ? vals[iz] : 0.0;
{% else %}
           f[MAKE3D(i0) + iz];
{% endif %}
        f[ic] = facs.f0 * val
{% for i in range(1,order) %}
           + facs.f{{i}} *f[MAKE3D(i{{i}}) + iz]
{% endfor %}
        ;
      }
    }
  }
}

"""
clone_str = """
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
stencil_str = """
vec{{order}} {{class}}::calc_interp_to_stencil(const vec{{order}} & spacing) const {
vec{{order}} facs;
// Stencil Code
{{stencil_code}}
return facs;
}
"""


if __name__ == "__main__":
    with maybeopen(sys.argv) as print:
        print(header)

        for order in orders:
            for boundary in boundaries:
                if boundary == "Neumann":
                    mat = sten.neumann
                else:
                    mat = sten.dirichlet
                try:
                    code = sten.gen_code(order, mat)
                except:
                    import sys

                    print("Order:", order, "boundary:", boundary, file=sys.stderr)
                    raise
                args = {
                    "order": order,
                    "type": boundary,
                    "class": "Boundary%sNonUniform_O%d" % (boundary, order),
                    "stencil_code": code,
                    "with_fg": boundary != "Free",
                }

                print(env.from_string(apply_str).render(**args))
                print(env.from_string(clone_str).render(**args))
                print(env.from_string(stencil_str).render(**args))
