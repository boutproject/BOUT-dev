#include <bout/constants.hxx>
#include <boundary_standard.hxx>
#include <boutexception.hxx>
#include <derivs.hxx>
#include <fft.hxx>
#include <globals.hxx>
#include <invert_laplace.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <utils.hxx>

#include "boundary_nonuniform.hxx"

static void update_bx_by_stagger(int& bx, int& by, int& stagger, CELL_LOC loc) {
  // NB: bx is going outwards
  // NB: XLOW means shifted in -x direction
  // `stagger` stagger direction with respect to direction of boundary
  //   0 : no stagger or orthogonal to boundary direction
  //   1 : staggerd in direction of boundary
  //  -1 : staggerd in oposite direction of boundary
  // Also note that all offsets are basically half a cell
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
}

void BoundaryDirichletNonUniform_O2::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int by = bndry->by;
  int stagger = 0;
  update_bx_by_stagger(bx, by, stagger, loc);
  int istart = (stagger == -1) ? -1 : 0;

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5
                         * (mesh->GlobalX(bndry->x)          // In the guard cell
                            + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5
                         * (mesh->GlobalY(bndry->y)          // In the guard cell
                            + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }

    vec2 spacing;
    vec2 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      spacing.f0 = 0;
      BoutReal st = 0;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = 0;
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
    }
    if (stagger == -1) {
      i1 = {bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    }
    for (int i = istart; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing.f0 += t;
        spacing.f1 += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        if (stagger == -1 && i != -1) {
          spacing += t;
        }
        facs = calc_interp_to_stencil(spacing);
        if (stagger == 1) {
          spacing += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, i1.z);

        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}

BoundaryOp* BoundaryDirichletNonUniform_O2::clone(BoundaryRegion* region,
                                                  const std::list<std::string>& args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichletNonUniform_O2(region, newgen);
}

vec2 BoundaryDirichletNonUniform_O2::calc_interp_to_stencil(const vec2& spacing) const {
  vec2 facs;
  // Stencil Code
  facs.f0 = -spacing.f1 / (spacing.f0 - spacing.f1);
  facs.f1 = spacing.f0 / (spacing.f0 - spacing.f1);

  return facs;
}

void BoundaryNeumannNonUniform_O2::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int by = bndry->by;
  int stagger = 0;
  update_bx_by_stagger(bx, by, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5
                         * (mesh->GlobalX(bndry->x)          // In the guard cell
                            + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5
                         * (mesh->GlobalY(bndry->y)          // In the guard cell
                            + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }

    vec2 spacing;
    vec2 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      spacing.f0 = 0;
      BoutReal st = 0;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = 0;
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
    }
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing.f0 += t;
        spacing.f1 += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        if (stagger == -1) {
          spacing += t;
        }
        facs = calc_interp_to_stencil(spacing);
        if (stagger == 1) {
          spacing += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, i1.z);

        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}

BoundaryOp* BoundaryNeumannNonUniform_O2::clone(BoundaryRegion* region,
                                                const std::list<std::string>& args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumannNonUniform_O2(region, newgen);
}

vec2 BoundaryNeumannNonUniform_O2::calc_interp_to_stencil(const vec2& spacing) const {
  vec2 facs;
  // Stencil Code
  facs.f0 = -spacing.f1;
  facs.f1 = 1;

  return facs;
}

void BoundaryFreeNonUniform_O2::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int by = bndry->by;
  int stagger = 0;
  update_bx_by_stagger(bx, by, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5
                         * (mesh->GlobalX(bndry->x)          // In the guard cell
                            + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5
                         * (mesh->GlobalY(bndry->y)          // In the guard cell
                            + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }

    vec2 spacing;
    vec2 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i0{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i1{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    if (stagger == 0) {
      BoutReal st = 0;
      t = coords_field(i0.x, i0.y);
      spacing.f0 = st + t / 2;
      st += t;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = 0;
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
    }

    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing.f0 += t;
        spacing.f1 += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        if (stagger == -1) {
          spacing += t;
        }
        facs = calc_interp_to_stencil(spacing);
        if (stagger == 1) {
          spacing += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        t = facs.f0 * f(i0.x, i0.y, i0.z) + facs.f1 * f(i1.x, i1.y, i1.z);

        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}

BoundaryOp* BoundaryFreeNonUniform_O2::clone(BoundaryRegion* region,
                                             const std::list<std::string>& args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryFreeNonUniform_O2(region, newgen);
}

vec2 BoundaryFreeNonUniform_O2::calc_interp_to_stencil(const vec2& spacing) const {
  vec2 facs;
  // Stencil Code
  facs.f0 = -spacing.f1 / (spacing.f0 - spacing.f1);
  facs.f1 = spacing.f0 / (spacing.f0 - spacing.f1);

  return facs;
}

void BoundaryDirichletNonUniform_O3::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int by = bndry->by;
  int stagger = 0;
  update_bx_by_stagger(bx, by, stagger, loc);
  int istart = (stagger == -1) ? -1 : 0;

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5
                         * (mesh->GlobalX(bndry->x)          // In the guard cell
                            + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5
                         * (mesh->GlobalY(bndry->y)          // In the guard cell
                            + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }

    vec3 spacing;
    vec3 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i2{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      spacing.f0 = 0;
      BoutReal st = 0;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
      t = coords_field(i2.x, i2.y);
      spacing.f2 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = 0;
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
      spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
    }
    if (stagger == -1) {
      i1 = {bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
      i2 = {bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    }
    for (int i = istart; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing.f0 += t;
        spacing.f1 += t;
        spacing.f2 += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        if (stagger == -1 && i != -1) {
          spacing += t;
        }
        facs = calc_interp_to_stencil(spacing);
        if (stagger == 1) {
          spacing += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, i1.z) + facs.f2 * f(i2.x, i2.y, i2.z);

        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}

BoundaryOp* BoundaryDirichletNonUniform_O3::clone(BoundaryRegion* region,
                                                  const std::list<std::string>& args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichletNonUniform_O3(region, newgen);
}

vec3 BoundaryDirichletNonUniform_O3::calc_interp_to_stencil(const vec3& spacing) const {
  vec3 facs;
  // Stencil Code
  facs.f0 =
      spacing.f1 * spacing.f2 / ((spacing.f0 - spacing.f1) * (spacing.f0 - spacing.f2));
  facs.f1 =
      -spacing.f0 * spacing.f2 / ((spacing.f0 - spacing.f1) * (spacing.f1 - spacing.f2));
  facs.f2 =
      spacing.f0 * spacing.f1 / ((spacing.f0 - spacing.f2) * (spacing.f1 - spacing.f2));

  return facs;
}

void BoundaryNeumannNonUniform_O3::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int by = bndry->by;
  int stagger = 0;
  update_bx_by_stagger(bx, by, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5
                         * (mesh->GlobalX(bndry->x)          // In the guard cell
                            + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5
                         * (mesh->GlobalY(bndry->y)          // In the guard cell
                            + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }

    vec3 spacing;
    vec3 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i2{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      spacing.f0 = 0;
      BoutReal st = 0;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
      t = coords_field(i2.x, i2.y);
      spacing.f2 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = 0;
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
      spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
    }
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing.f0 += t;
        spacing.f1 += t;
        spacing.f2 += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        if (stagger == -1) {
          spacing += t;
        }
        facs = calc_interp_to_stencil(spacing);
        if (stagger == 1) {
          spacing += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, i1.z) + facs.f2 * f(i2.x, i2.y, i2.z);

        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}

BoundaryOp* BoundaryNeumannNonUniform_O3::clone(BoundaryRegion* region,
                                                const std::list<std::string>& args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumannNonUniform_O3(region, newgen);
}

vec3 BoundaryNeumannNonUniform_O3::calc_interp_to_stencil(const vec3& spacing) const {
  vec3 facs;
  // Stencil Code
  facs.f0 = spacing.f1 * spacing.f2 / (2 * spacing.f0 - spacing.f1 - spacing.f2);
  facs.f1 = -spacing.f2 * (2 * spacing.f0 - spacing.f2)
            / ((spacing.f1 - spacing.f2) * (2 * spacing.f0 - spacing.f1 - spacing.f2));
  facs.f2 = spacing.f1 * (2 * spacing.f0 - spacing.f1)
            / ((spacing.f1 - spacing.f2) * (2 * spacing.f0 - spacing.f1 - spacing.f2));

  return facs;
}

void BoundaryFreeNonUniform_O3::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int by = bndry->by;
  int stagger = 0;
  update_bx_by_stagger(bx, by, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5
                         * (mesh->GlobalX(bndry->x)          // In the guard cell
                            + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5
                         * (mesh->GlobalY(bndry->y)          // In the guard cell
                            + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }

    vec3 spacing;
    vec3 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i0{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i1{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    Indices i2{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    if (stagger == 0) {
      BoutReal st = 0;
      t = coords_field(i0.x, i0.y);
      spacing.f0 = st + t / 2;
      st += t;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
      t = coords_field(i2.x, i2.y);
      spacing.f2 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = 0;
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
      spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
    }

    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing.f0 += t;
        spacing.f1 += t;
        spacing.f2 += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        if (stagger == -1) {
          spacing += t;
        }
        facs = calc_interp_to_stencil(spacing);
        if (stagger == 1) {
          spacing += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        t = facs.f0 * f(i0.x, i0.y, i0.z) + facs.f1 * f(i1.x, i1.y, i1.z)
            + facs.f2 * f(i2.x, i2.y, i2.z);

        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}

BoundaryOp* BoundaryFreeNonUniform_O3::clone(BoundaryRegion* region,
                                             const std::list<std::string>& args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryFreeNonUniform_O3(region, newgen);
}

vec3 BoundaryFreeNonUniform_O3::calc_interp_to_stencil(const vec3& spacing) const {
  vec3 facs;
  // Stencil Code
  facs.f0 =
      spacing.f1 * spacing.f2 / ((spacing.f0 - spacing.f1) * (spacing.f0 - spacing.f2));
  facs.f1 =
      -spacing.f0 * spacing.f2 / ((spacing.f0 - spacing.f1) * (spacing.f1 - spacing.f2));
  facs.f2 =
      spacing.f0 * spacing.f1 / ((spacing.f0 - spacing.f2) * (spacing.f1 - spacing.f2));

  return facs;
}

void BoundaryDirichletNonUniform_O4::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int by = bndry->by;
  int stagger = 0;
  update_bx_by_stagger(bx, by, stagger, loc);
  int istart = (stagger == -1) ? -1 : 0;

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5
                         * (mesh->GlobalX(bndry->x)          // In the guard cell
                            + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5
                         * (mesh->GlobalY(bndry->y)          // In the guard cell
                            + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }

    vec4 spacing;
    vec4 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i2{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    Indices i3{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      spacing.f0 = 0;
      BoutReal st = 0;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
      t = coords_field(i2.x, i2.y);
      spacing.f2 = st + t / 2;
      st += t;
      t = coords_field(i3.x, i3.y);
      spacing.f3 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = 0;
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
      spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
      spacing.f3 = spacing.f2 + coords_field(i3.x, i3.y);
    }
    if (stagger == -1) {
      i1 = {bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
      i2 = {bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
      i3 = {bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by, 0};
    }
    for (int i = istart; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing.f0 += t;
        spacing.f1 += t;
        spacing.f2 += t;
        spacing.f3 += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        if (stagger == -1 && i != -1) {
          spacing += t;
        }
        facs = calc_interp_to_stencil(spacing);
        if (stagger == 1) {
          spacing += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        i3.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, i1.z) + facs.f2 * f(i2.x, i2.y, i2.z)
            + facs.f3 * f(i3.x, i3.y, i3.z);

        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}

BoundaryOp* BoundaryDirichletNonUniform_O4::clone(BoundaryRegion* region,
                                                  const std::list<std::string>& args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichletNonUniform_O4(region, newgen);
}

vec4 BoundaryDirichletNonUniform_O4::calc_interp_to_stencil(const vec4& spacing) const {
  vec4 facs;
  // Stencil Code
  facs.f0 = -spacing.f1 * spacing.f2 * spacing.f3
            / ((spacing.f0 - spacing.f1) * (spacing.f0 - spacing.f2)
               * (spacing.f0 - spacing.f3));
  facs.f1 = spacing.f0 * spacing.f2 * spacing.f3
            / ((spacing.f0 - spacing.f1) * (spacing.f1 - spacing.f2)
               * (spacing.f1 - spacing.f3));
  facs.f2 = -spacing.f0 * spacing.f1 * spacing.f3
            / ((spacing.f0 - spacing.f2) * (spacing.f1 - spacing.f2)
               * (spacing.f2 - spacing.f3));
  facs.f3 = spacing.f0 * spacing.f1 * spacing.f2
            / ((spacing.f0 - spacing.f3) * (spacing.f1 - spacing.f3)
               * (spacing.f2 - spacing.f3));

  return facs;
}

void BoundaryNeumannNonUniform_O4::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int by = bndry->by;
  int stagger = 0;
  update_bx_by_stagger(bx, by, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5
                         * (mesh->GlobalX(bndry->x)          // In the guard cell
                            + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5
                         * (mesh->GlobalY(bndry->y)          // In the guard cell
                            + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }

    vec4 spacing;
    vec4 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i2{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    Indices i3{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      spacing.f0 = 0;
      BoutReal st = 0;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
      t = coords_field(i2.x, i2.y);
      spacing.f2 = st + t / 2;
      st += t;
      t = coords_field(i3.x, i3.y);
      spacing.f3 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = 0;
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
      spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
      spacing.f3 = spacing.f2 + coords_field(i3.x, i3.y);
    }
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing.f0 += t;
        spacing.f1 += t;
        spacing.f2 += t;
        spacing.f3 += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        if (stagger == -1) {
          spacing += t;
        }
        facs = calc_interp_to_stencil(spacing);
        if (stagger == 1) {
          spacing += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        i3.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, i1.z) + facs.f2 * f(i2.x, i2.y, i2.z)
            + facs.f3 * f(i3.x, i3.y, i3.z);

        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}

BoundaryOp* BoundaryNeumannNonUniform_O4::clone(BoundaryRegion* region,
                                                const std::list<std::string>& args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumannNonUniform_O4(region, newgen);
}

vec4 BoundaryNeumannNonUniform_O4::calc_interp_to_stencil(const vec4& spacing) const {
  vec4 facs;
  // Stencil Code
  facs.f0 =
      -spacing.f1 * spacing.f2 * spacing.f3
      / (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
         - 2 * spacing.f0 * spacing.f2 - 2 * spacing.f0 * spacing.f3
         + spacing.f1 * spacing.f2 + spacing.f1 * spacing.f3 + spacing.f2 * spacing.f3);
  facs.f1 = spacing.f2 * spacing.f3
            * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f2
               - 2 * spacing.f0 * spacing.f3 + spacing.f2 * spacing.f3)
            / ((spacing.f1 - spacing.f2) * (spacing.f1 - spacing.f3)
               * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
                  - 2 * spacing.f0 * spacing.f2 - 2 * spacing.f0 * spacing.f3
                  + spacing.f1 * spacing.f2 + spacing.f1 * spacing.f3
                  + spacing.f2 * spacing.f3));
  facs.f2 = -spacing.f1 * spacing.f3
            * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
               - 2 * spacing.f0 * spacing.f3 + spacing.f1 * spacing.f3)
            / ((spacing.f1 - spacing.f2) * (spacing.f2 - spacing.f3)
               * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
                  - 2 * spacing.f0 * spacing.f2 - 2 * spacing.f0 * spacing.f3
                  + spacing.f1 * spacing.f2 + spacing.f1 * spacing.f3
                  + spacing.f2 * spacing.f3));
  facs.f3 = spacing.f1 * spacing.f2
            * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
               - 2 * spacing.f0 * spacing.f2 + spacing.f1 * spacing.f2)
            / ((spacing.f1 - spacing.f3) * (spacing.f2 - spacing.f3)
               * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
                  - 2 * spacing.f0 * spacing.f2 - 2 * spacing.f0 * spacing.f3
                  + spacing.f1 * spacing.f2 + spacing.f1 * spacing.f3
                  + spacing.f2 * spacing.f3));

  return facs;
}

void BoundaryFreeNonUniform_O4::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  BoutReal val = 0.0;

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int by = bndry->by;
  int stagger = 0;
  update_bx_by_stagger(bx, by, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        // Calculate the X and Y normalised values half-way between the guard cell and
        // grid cell
        BoutReal xnorm = 0.5
                         * (mesh->GlobalX(bndry->x)          // In the guard cell
                            + mesh->GlobalX(bndry->x - bx)); // the grid cell

        BoutReal ynorm = 0.5
                         * (mesh->GlobalY(bndry->y)          // In the guard cell
                            + mesh->GlobalY(bndry->y - by)); // the grid cell

        vals[zk] = fg->generate(xnorm, TWOPI * ynorm, TWOPI * zk / (mesh->LocalNz), t);
      }
    }

    vec4 spacing;
    vec4 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i0{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i1{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    Indices i2{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    Indices i3{bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by, 0};
    if (stagger == 0) {
      BoutReal st = 0;
      t = coords_field(i0.x, i0.y);
      spacing.f0 = st + t / 2;
      st += t;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
      t = coords_field(i2.x, i2.y);
      spacing.f2 = st + t / 2;
      st += t;
      t = coords_field(i3.x, i3.y);
      spacing.f3 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = 0;
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
      spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
      spacing.f3 = spacing.f2 + coords_field(i3.x, i3.y);
    }

    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing.f0 += t;
        spacing.f1 += t;
        spacing.f2 += t;
        spacing.f3 += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        if (stagger == -1) {
          spacing += t;
        }
        facs = calc_interp_to_stencil(spacing);
        if (stagger == 1) {
          spacing += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        i3.z = ic.z;
        t = facs.f0 * f(i0.x, i0.y, i0.z) + facs.f1 * f(i1.x, i1.y, i1.z)
            + facs.f2 * f(i2.x, i2.y, i2.z) + facs.f3 * f(i3.x, i3.y, i3.z);

        f(ic.x, ic.y, ic.z) = t;
      }
    }
  }
}

BoundaryOp* BoundaryFreeNonUniform_O4::clone(BoundaryRegion* region,
                                             const std::list<std::string>& args) {

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryFreeNonUniform_O4(region, newgen);
}

vec4 BoundaryFreeNonUniform_O4::calc_interp_to_stencil(const vec4& spacing) const {
  vec4 facs;
  // Stencil Code
  facs.f0 = -spacing.f1 * spacing.f2 * spacing.f3
            / ((spacing.f0 - spacing.f1) * (spacing.f0 - spacing.f2)
               * (spacing.f0 - spacing.f3));
  facs.f1 = spacing.f0 * spacing.f2 * spacing.f3
            / ((spacing.f0 - spacing.f1) * (spacing.f1 - spacing.f2)
               * (spacing.f1 - spacing.f3));
  facs.f2 = -spacing.f0 * spacing.f1 * spacing.f3
            / ((spacing.f0 - spacing.f2) * (spacing.f1 - spacing.f2)
               * (spacing.f2 - spacing.f3));
  facs.f3 = spacing.f0 * spacing.f1 * spacing.f2
            / ((spacing.f0 - spacing.f3) * (spacing.f1 - spacing.f3)
               * (spacing.f2 - spacing.f3));

  return facs;
}
