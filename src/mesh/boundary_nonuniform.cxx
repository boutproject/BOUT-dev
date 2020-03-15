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

static void update_stagger_offsets(int& x_boundary_offset, int& y_boundary_offset,
                                   int& stagger, CELL_LOC loc) {
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

void BoundaryDirichletNonUniform_O2::apply(Field3D& f, BoutReal t) {
  bndry->first();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  std::vector<BoutReal> vals;
  vals.reserve(mesh->LocalNz);

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      // Calculate the X and Y normalised values half-way between the guard cell and
      // grid cell
      const BoutReal xnorm =
          0.5
          * (mesh->GlobalX(bndry->x)                         // In the guard cell
             + mesh->GlobalX(bndry->x - x_boundary_offset)); // the grid cell
      const BoutReal ynorm =
          TWOPI * 0.5
          * (mesh->GlobalY(bndry->y)                         // In the guard cell
             + mesh->GlobalY(bndry->y - y_boundary_offset)); // the grid cell
      const BoutReal zfac = TWOPI / mesh->LocalNz;
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        vals[zk] = fg->generate(bout::generator::Context().set("x", xnorm, "y", ynorm,
                                                               "z", zfac * zk, "t", t));
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
    // with dirichlet, we specify the value on the boundary, even if
    // the value is part of the evolving system.
    for (int i = ((stagger == -1) ? -1 : 0); i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing += t;
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
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        const BoutReal val = (fg) ? vals[iz] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, iz);

        f(ic.x, ic.y, iz) = t;
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

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  std::vector<BoutReal> vals;
  vals.reserve(mesh->LocalNz);

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      // Calculate the X and Y normalised values half-way between the guard cell and
      // grid cell
      const BoutReal xnorm =
          0.5
          * (mesh->GlobalX(bndry->x)                         // In the guard cell
             + mesh->GlobalX(bndry->x - x_boundary_offset)); // the grid cell
      const BoutReal ynorm =
          TWOPI * 0.5
          * (mesh->GlobalY(bndry->y)                         // In the guard cell
             + mesh->GlobalY(bndry->y - y_boundary_offset)); // the grid cell
      const BoutReal zfac = TWOPI / mesh->LocalNz;
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        vals[zk] = fg->generate(bout::generator::Context().set("x", xnorm, "y", ynorm,
                                                               "z", zfac * zk, "t", t));
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
      // Check if we are staggered and also boundary in low
      //  direction
      // In the case of Neumann we have in this case two values
      //  defined at the same point
      if (stagger == -1
          && ((bndry->bx && x_boundary_offset == -1)
              || (bndry->by && y_boundary_offset == -1))) {
        spacing.f1 = spacing.f0;
      } else {
        spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
      }
    }
    // With free and neumann the value is not set if the point is
    // evolved and it is on the boundary.
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing += t;
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
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        const BoutReal val = (fg) ? vals[iz] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, iz);

        f(ic.x, ic.y, iz) = t;
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

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {

    vec2 spacing;
    vec2 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    const Indices i0{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    const Indices i1{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    if (stagger == 0) {
      BoutReal st = 0;
      t = coords_field(i0.x, i0.y);
      spacing.f0 = st + t / 2;
      st += t;
      t = coords_field(i1.x, i1.y);
      spacing.f1 = st + t / 2;
      st += t;
    } else {
      spacing.f0 = coords_field(i0.x, i0.y);
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
    }

    // With free and neumann the value is not set if the point is
    // evolved and it is on the boundary.
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      }
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        t = facs.f0 * f(i0.x, i0.y, iz) + facs.f1 * f(i1.x, i1.y, iz);

        f(ic.x, ic.y, iz) = t;
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

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  std::vector<BoutReal> vals;
  vals.reserve(mesh->LocalNz);

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      // Calculate the X and Y normalised values half-way between the guard cell and
      // grid cell
      const BoutReal xnorm =
          0.5
          * (mesh->GlobalX(bndry->x)                         // In the guard cell
             + mesh->GlobalX(bndry->x - x_boundary_offset)); // the grid cell
      const BoutReal ynorm =
          TWOPI * 0.5
          * (mesh->GlobalY(bndry->y)                         // In the guard cell
             + mesh->GlobalY(bndry->y - y_boundary_offset)); // the grid cell
      const BoutReal zfac = TWOPI / mesh->LocalNz;
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        vals[zk] = fg->generate(bout::generator::Context().set("x", xnorm, "y", ynorm,
                                                               "z", zfac * zk, "t", t));
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
    // with dirichlet, we specify the value on the boundary, even if
    // the value is part of the evolving system.
    for (int i = ((stagger == -1) ? -1 : 0); i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing += t;
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
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        const BoutReal val = (fg) ? vals[iz] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, iz) + facs.f2 * f(i2.x, i2.y, iz);

        f(ic.x, ic.y, iz) = t;
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

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  std::vector<BoutReal> vals;
  vals.reserve(mesh->LocalNz);

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      // Calculate the X and Y normalised values half-way between the guard cell and
      // grid cell
      const BoutReal xnorm =
          0.5
          * (mesh->GlobalX(bndry->x)                         // In the guard cell
             + mesh->GlobalX(bndry->x - x_boundary_offset)); // the grid cell
      const BoutReal ynorm =
          TWOPI * 0.5
          * (mesh->GlobalY(bndry->y)                         // In the guard cell
             + mesh->GlobalY(bndry->y - y_boundary_offset)); // the grid cell
      const BoutReal zfac = TWOPI / mesh->LocalNz;
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        vals[zk] = fg->generate(bout::generator::Context().set("x", xnorm, "y", ynorm,
                                                               "z", zfac * zk, "t", t));
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
      // Check if we are staggered and also boundary in low
      //  direction
      // In the case of Neumann we have in this case two values
      //  defined at the same point
      if (stagger == -1
          && ((bndry->bx && x_boundary_offset == -1)
              || (bndry->by && y_boundary_offset == -1))) {
        spacing.f1 = spacing.f0;
        spacing.f2 = spacing.f1 + coords_field(i1.x, i1.y);
      } else {
        spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
        spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
      }
    }
    // With free and neumann the value is not set if the point is
    // evolved and it is on the boundary.
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing += t;
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
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        const BoutReal val = (fg) ? vals[iz] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, iz) + facs.f2 * f(i2.x, i2.y, iz);

        f(ic.x, ic.y, iz) = t;
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

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {

    vec3 spacing;
    vec3 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    const Indices i0{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    const Indices i1{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    const Indices i2{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
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
      spacing.f0 = coords_field(i0.x, i0.y);
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
      spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
    }

    // With free and neumann the value is not set if the point is
    // evolved and it is on the boundary.
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      }
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        t = facs.f0 * f(i0.x, i0.y, iz) + facs.f1 * f(i1.x, i1.y, iz)
            + facs.f2 * f(i2.x, i2.y, iz);

        f(ic.x, ic.y, iz) = t;
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

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  std::vector<BoutReal> vals;
  vals.reserve(mesh->LocalNz);

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      // Calculate the X and Y normalised values half-way between the guard cell and
      // grid cell
      const BoutReal xnorm =
          0.5
          * (mesh->GlobalX(bndry->x)                         // In the guard cell
             + mesh->GlobalX(bndry->x - x_boundary_offset)); // the grid cell
      const BoutReal ynorm =
          TWOPI * 0.5
          * (mesh->GlobalY(bndry->y)                         // In the guard cell
             + mesh->GlobalY(bndry->y - y_boundary_offset)); // the grid cell
      const BoutReal zfac = TWOPI / mesh->LocalNz;
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        vals[zk] = fg->generate(bout::generator::Context().set("x", xnorm, "y", ynorm,
                                                               "z", zfac * zk, "t", t));
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
    // with dirichlet, we specify the value on the boundary, even if
    // the value is part of the evolving system.
    for (int i = ((stagger == -1) ? -1 : 0); i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing += t;
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
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        const BoutReal val = (fg) ? vals[iz] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, iz) + facs.f2 * f(i2.x, i2.y, iz)
            + facs.f3 * f(i3.x, i3.y, iz);

        f(ic.x, ic.y, iz) = t;
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

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  std::vector<BoutReal> vals;
  vals.reserve(mesh->LocalNz);

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {
    if (fg) {
      // Calculate the X and Y normalised values half-way between the guard cell and
      // grid cell
      const BoutReal xnorm =
          0.5
          * (mesh->GlobalX(bndry->x)                         // In the guard cell
             + mesh->GlobalX(bndry->x - x_boundary_offset)); // the grid cell
      const BoutReal ynorm =
          TWOPI * 0.5
          * (mesh->GlobalY(bndry->y)                         // In the guard cell
             + mesh->GlobalY(bndry->y - y_boundary_offset)); // the grid cell
      const BoutReal zfac = TWOPI / mesh->LocalNz;
      for (int zk = 0; zk < mesh->LocalNz; zk++) {
        vals[zk] = fg->generate(bout::generator::Context().set("x", xnorm, "y", ynorm,
                                                               "z", zfac * zk, "t", t));
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
      // Check if we are staggered and also boundary in low
      //  direction
      // In the case of Neumann we have in this case two values
      //  defined at the same point
      if (stagger == -1
          && ((bndry->bx && x_boundary_offset == -1)
              || (bndry->by && y_boundary_offset == -1))) {
        spacing.f1 = spacing.f0;
        spacing.f2 = spacing.f1 + coords_field(i1.x, i1.y);
        spacing.f3 = spacing.f2 + coords_field(i2.x, i2.y);
      } else {
        spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
        spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
        spacing.f3 = spacing.f2 + coords_field(i3.x, i3.y);
      }
    }
    // With free and neumann the value is not set if the point is
    // evolved and it is on the boundary.
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing += t;
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
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        const BoutReal val = (fg) ? vals[iz] : 0.0;
        t = facs.f0 * val + facs.f1 * f(i1.x, i1.y, iz) + facs.f2 * f(i2.x, i2.y, iz)
            + facs.f3 * f(i3.x, i3.y, iz);

        f(ic.x, ic.y, iz) = t;
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

  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  int x_boundary_offset = bndry->bx;
  int y_boundary_offset = bndry->by;
  int stagger = 0;
  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

  for (; !bndry->isDone(); bndry->next1d()) {

    vec4 spacing;
    vec4 facs;

    const Field2D& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    const Indices i0{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    const Indices i1{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    const Indices i2{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    const Indices i3{bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by, 0};
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
      spacing.f0 = coords_field(i0.x, i0.y);
      spacing.f1 = spacing.f0 + coords_field(i1.x, i1.y);
      spacing.f2 = spacing.f1 + coords_field(i2.x, i2.y);
      spacing.f3 = spacing.f2 + coords_field(i3.x, i3.y);
    }

    // With free and neumann the value is not set if the point is
    // evolved and it is on the boundary.
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = coords_field(ic.x, ic.y) / 2;
        spacing += t;
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      } else {
        t = coords_field(ic.x, ic.y);
        facs = calc_interp_to_stencil(spacing);
        spacing += t;
      }
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        t = facs.f0 * f(i0.x, i0.y, iz) + facs.f1 * f(i1.x, i1.y, iz)
            + facs.f2 * f(i2.x, i2.y, iz) + facs.f3 * f(i3.x, i3.y, iz);

        f(ic.x, ic.y, iz) = t;
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
