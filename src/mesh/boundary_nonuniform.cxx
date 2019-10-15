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
    BoutReal fac0;
    BoutReal x0;
    BoutReal fac1;
    BoutReal x1;

    const Field2D& dy =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      x0 = 0;
      BoutReal st = 0;
      t = dy[i1];
      x1 = st + t / 2;
      st += t;
    } else {
      x0 = 0; // dy(bndry->x, bndry->y) / 2;
      x1 = x0 + dy[i1];
    }
    if (stagger == -1) {
      i1 = {bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    }
    for (int i = istart; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, fac0, fac1);
        x0 += t;
        x1 += t;
      } else {
        t = dy[ic];
        if (stagger == -1 && i != -1) {
          x0 += t;
          x1 += t;
        }
        calc_interp_to_stencil(x0, x1, fac0, fac1);
        if (stagger == 1) {
          x0 += t;
          x1 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = fac0 * val + fac1 * f[i1];

        f[ic] = t;
      }
    }
  }
}

BoundaryOp* BoundaryDirichletNonUniform_O2::clone(BoundaryRegion* region,
                                                  const std::list<std::string>& args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichletNonUniform_O2(region, newgen);
}

void BoundaryDirichletNonUniform_O2::calc_interp_to_stencil(BoutReal x0, BoutReal x1,
                                                            BoutReal& fac0,
                                                            BoutReal& fac1) const {
  // Stencil Code
  fac0 = -x1 / (x0 - x1);
  fac1 = x0 / (x0 - x1);
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
    BoutReal fac0;
    BoutReal x0;
    BoutReal fac1;
    BoutReal x1;

    const Field2D& dy =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      x0 = 0;
      BoutReal st = 0;
      t = dy[i1];
      x1 = st + t / 2;
      st += t;
    } else {
      x0 = 0; // dy(bndry->x, bndry->y) / 2;
      x1 = x0 + dy[i1];
    }
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, fac0, fac1);
        x0 += t;
        x1 += t;
      } else {
        t = dy[ic];
        if (stagger == -1) {
          x0 += t;
          x1 += t;
        }
        calc_interp_to_stencil(x0, x1, fac0, fac1);
        if (stagger == 1) {
          x0 += t;
          x1 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = fac0 * val + fac1 * f[i1];

        f[ic] = t;
      }
    }
  }
}

BoundaryOp* BoundaryNeumannNonUniform_O2::clone(BoundaryRegion* region,
                                                const std::list<std::string>& args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumannNonUniform_O2(region, newgen);
}

void BoundaryNeumannNonUniform_O2::calc_interp_to_stencil(BoutReal x0, BoutReal x1,
                                                          BoutReal& fac0,
                                                          BoutReal& fac1) const {
  // Stencil Code
  fac0 = -x1;
  fac1 = 1;
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
    BoutReal fac0;
    BoutReal x0;
    BoutReal fac1;
    BoutReal x1;

    const Field2D& dy =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i0{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i1{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    if (stagger == 0) {
      BoutReal st = 0;
      t = dy[i0];
      x0 = st + t / 2;
      st += t;
      t = dy[i1];
      x1 = st + t / 2;
      st += t;
    } else {
      x0 = 0;
      x1 = x0 + dy[i1];
    }

    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, fac0, fac1);
        x0 += t;
        x1 += t;
      } else {
        t = dy[ic];
        if (stagger == -1) {
          x0 += t;
          x1 += t;
        }
        calc_interp_to_stencil(x0, x1, fac0, fac1);
        if (stagger == 1) {
          x0 += t;
          x1 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        t = fac0 * f[i0] + fac1 * f[i1];

        f[ic] = t;
      }
    }
  }
}

BoundaryOp* BoundaryFreeNonUniform_O2::clone(BoundaryRegion* region,
                                             const std::list<std::string>& args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryFreeNonUniform_O2(region, newgen);
}

void BoundaryFreeNonUniform_O2::calc_interp_to_stencil(BoutReal x0, BoutReal x1,
                                                       BoutReal& fac0,
                                                       BoutReal& fac1) const {
  // Stencil Code
  fac0 = -x1 / (x0 - x1);
  fac1 = x0 / (x0 - x1);
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
    BoutReal fac0;
    BoutReal x0;
    BoutReal fac1;
    BoutReal x1;
    BoutReal fac2;
    BoutReal x2;

    const Field2D& dy =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i2{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      x0 = 0;
      BoutReal st = 0;
      t = dy[i1];
      x1 = st + t / 2;
      st += t;
      t = dy[i2];
      x2 = st + t / 2;
      st += t;
    } else {
      x0 = 0; // dy(bndry->x, bndry->y) / 2;
      x1 = x0 + dy[i1];
      x2 = x1 + dy[i2];
    }
    if (stagger == -1) {
      i1 = {bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
      i2 = {bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    }
    for (int i = istart; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        x2 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, x2, fac0, fac1, fac2);
        x0 += t;
        x1 += t;
        x2 += t;
      } else {
        t = dy[ic];
        if (stagger == -1 && i != -1) {
          x0 += t;
          x1 += t;
          x2 += t;
        }
        calc_interp_to_stencil(x0, x1, x2, fac0, fac1, fac2);
        if (stagger == 1) {
          x0 += t;
          x1 += t;
          x2 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = fac0 * val + fac1 * f[i1] + fac2 * f[i2];

        f[ic] = t;
      }
    }
  }
}

BoundaryOp* BoundaryDirichletNonUniform_O3::clone(BoundaryRegion* region,
                                                  const std::list<std::string>& args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichletNonUniform_O3(region, newgen);
}

void BoundaryDirichletNonUniform_O3::calc_interp_to_stencil(BoutReal x0, BoutReal x1,
                                                            BoutReal x2, BoutReal& fac0,
                                                            BoutReal& fac1,
                                                            BoutReal& fac2) const {
  // Stencil Code
  fac0 = x1 * x2 / ((x0 - x1) * (x0 - x2));
  fac1 = -x0 * x2 / ((x0 - x1) * (x1 - x2));
  fac2 = x0 * x1 / ((x0 - x2) * (x1 - x2));
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
    BoutReal fac0;
    BoutReal x0;
    BoutReal fac1;
    BoutReal x1;
    BoutReal fac2;
    BoutReal x2;

    const Field2D& dy =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i2{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      x0 = 0;
      BoutReal st = 0;
      t = dy[i1];
      x1 = st + t / 2;
      st += t;
      t = dy[i2];
      x2 = st + t / 2;
      st += t;
    } else {
      x0 = 0; // dy(bndry->x, bndry->y) / 2;
      x1 = x0 + dy[i1];
      x2 = x1 + dy[i2];
    }
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        x2 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, x2, fac0, fac1, fac2);
        x0 += t;
        x1 += t;
        x2 += t;
      } else {
        t = dy[ic];
        if (stagger == -1) {
          x0 += t;
          x1 += t;
          x2 += t;
        }
        calc_interp_to_stencil(x0, x1, x2, fac0, fac1, fac2);
        if (stagger == 1) {
          x0 += t;
          x1 += t;
          x2 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = fac0 * val + fac1 * f[i1] + fac2 * f[i2];

        f[ic] = t;
      }
    }
  }
}

BoundaryOp* BoundaryNeumannNonUniform_O3::clone(BoundaryRegion* region,
                                                const std::list<std::string>& args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumannNonUniform_O3(region, newgen);
}

void BoundaryNeumannNonUniform_O3::calc_interp_to_stencil(BoutReal x0, BoutReal x1,
                                                          BoutReal x2, BoutReal& fac0,
                                                          BoutReal& fac1,
                                                          BoutReal& fac2) const {
  // Stencil Code
  fac0 = x1 * x2 / (2 * x0 - x1 - x2);
  fac1 = -x2 * (2 * x0 - x2) / ((x1 - x2) * (2 * x0 - x1 - x2));
  fac2 = x1 * (2 * x0 - x1) / ((x1 - x2) * (2 * x0 - x1 - x2));
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
    BoutReal fac0;
    BoutReal x0;
    BoutReal fac1;
    BoutReal x1;
    BoutReal fac2;
    BoutReal x2;

    const Field2D& dy =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i0{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i1{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    Indices i2{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    if (stagger == 0) {
      BoutReal st = 0;
      t = dy[i0];
      x0 = st + t / 2;
      st += t;
      t = dy[i1];
      x1 = st + t / 2;
      st += t;
      t = dy[i2];
      x2 = st + t / 2;
      st += t;
    } else {
      x0 = 0;
      x1 = x0 + dy[i1];
      x2 = x1 + dy[i2];
    }

    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        x2 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, x2, fac0, fac1, fac2);
        x0 += t;
        x1 += t;
        x2 += t;
      } else {
        t = dy[ic];
        if (stagger == -1) {
          x0 += t;
          x1 += t;
          x2 += t;
        }
        calc_interp_to_stencil(x0, x1, x2, fac0, fac1, fac2);
        if (stagger == 1) {
          x0 += t;
          x1 += t;
          x2 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        t = fac0 * f[i0] + fac1 * f[i1] + fac2 * f[i2];

        f[ic] = t;
      }
    }
  }
}

BoundaryOp* BoundaryFreeNonUniform_O3::clone(BoundaryRegion* region,
                                             const std::list<std::string>& args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryFreeNonUniform_O3(region, newgen);
}

void BoundaryFreeNonUniform_O3::calc_interp_to_stencil(BoutReal x0, BoutReal x1,
                                                       BoutReal x2, BoutReal& fac0,
                                                       BoutReal& fac1,
                                                       BoutReal& fac2) const {
  // Stencil Code
  fac0 = x1 * x2 / ((x0 - x1) * (x0 - x2));
  fac1 = -x0 * x2 / ((x0 - x1) * (x1 - x2));
  fac2 = x0 * x1 / ((x0 - x2) * (x1 - x2));
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
    BoutReal fac0;
    BoutReal x0;
    BoutReal fac1;
    BoutReal x1;
    BoutReal fac2;
    BoutReal x2;
    BoutReal fac3;
    BoutReal x3;

    const Field2D& dy =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i2{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    Indices i3{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      x0 = 0;
      BoutReal st = 0;
      t = dy[i1];
      x1 = st + t / 2;
      st += t;
      t = dy[i2];
      x2 = st + t / 2;
      st += t;
      t = dy[i3];
      x3 = st + t / 2;
      st += t;
    } else {
      x0 = 0; // dy(bndry->x, bndry->y) / 2;
      x1 = x0 + dy[i1];
      x2 = x1 + dy[i2];
      x3 = x2 + dy[i3];
    }
    if (stagger == -1) {
      i1 = {bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
      i2 = {bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
      i3 = {bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by, 0};
    }
    for (int i = istart; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        x2 += t;
        x3 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, x2, x3, fac0, fac1, fac2, fac3);
        x0 += t;
        x1 += t;
        x2 += t;
        x3 += t;
      } else {
        t = dy[ic];
        if (stagger == -1 && i != -1) {
          x0 += t;
          x1 += t;
          x2 += t;
          x3 += t;
        }
        calc_interp_to_stencil(x0, x1, x2, x3, fac0, fac1, fac2, fac3);
        if (stagger == 1) {
          x0 += t;
          x1 += t;
          x2 += t;
          x3 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        i3.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = fac0 * val + fac1 * f[i1] + fac2 * f[i2] + fac3 * f[i3];

        f[ic] = t;
      }
    }
  }
}

BoundaryOp* BoundaryDirichletNonUniform_O4::clone(BoundaryRegion* region,
                                                  const std::list<std::string>& args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichletNonUniform_O4(region, newgen);
}

void BoundaryDirichletNonUniform_O4::calc_interp_to_stencil(
    BoutReal x0, BoutReal x1, BoutReal x2, BoutReal x3, BoutReal& fac0, BoutReal& fac1,
    BoutReal& fac2, BoutReal& fac3) const {
  // Stencil Code
  fac0 = -x1 * x2 * x3 / ((x0 - x1) * (x0 - x2) * (x0 - x3));
  fac1 = x0 * x2 * x3 / ((x0 - x1) * (x1 - x2) * (x1 - x3));
  fac2 = -x0 * x1 * x3 / ((x0 - x2) * (x1 - x2) * (x2 - x3));
  fac3 = x0 * x1 * x2 / ((x0 - x3) * (x1 - x3) * (x2 - x3));
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
    BoutReal fac0;
    BoutReal x0;
    BoutReal fac1;
    BoutReal x1;
    BoutReal fac2;
    BoutReal x2;
    BoutReal fac3;
    BoutReal x3;

    const Field2D& dy =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i2{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    Indices i3{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    BoutReal t;
    if (stagger == 0) {
      x0 = 0;
      BoutReal st = 0;
      t = dy[i1];
      x1 = st + t / 2;
      st += t;
      t = dy[i2];
      x2 = st + t / 2;
      st += t;
      t = dy[i3];
      x3 = st + t / 2;
      st += t;
    } else {
      x0 = 0; // dy(bndry->x, bndry->y) / 2;
      x1 = x0 + dy[i1];
      x2 = x1 + dy[i2];
      x3 = x2 + dy[i3];
    }
    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        x2 += t;
        x3 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, x2, x3, fac0, fac1, fac2, fac3);
        x0 += t;
        x1 += t;
        x2 += t;
        x3 += t;
      } else {
        t = dy[ic];
        if (stagger == -1) {
          x0 += t;
          x1 += t;
          x2 += t;
          x3 += t;
        }
        calc_interp_to_stencil(x0, x1, x2, x3, fac0, fac1, fac2, fac3);
        if (stagger == 1) {
          x0 += t;
          x1 += t;
          x2 += t;
          x3 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        i3.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        t = fac0 * val + fac1 * f[i1] + fac2 * f[i2] + fac3 * f[i3];

        f[ic] = t;
      }
    }
  }
}

BoundaryOp* BoundaryNeumannNonUniform_O4::clone(BoundaryRegion* region,
                                                const std::list<std::string>& args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryNeumannNonUniform_O4(region, newgen);
}

void BoundaryNeumannNonUniform_O4::calc_interp_to_stencil(BoutReal x0, BoutReal x1,
                                                          BoutReal x2, BoutReal x3,
                                                          BoutReal& fac0, BoutReal& fac1,
                                                          BoutReal& fac2,
                                                          BoutReal& fac3) const {
  // Stencil Code
  fac0 = -x1 * x2 * x3
         / (3 * pow(x0, 2) - 2 * x0 * x1 - 2 * x0 * x2 - 2 * x0 * x3 + x1 * x2 + x1 * x3
            + x2 * x3);
  fac1 = x2 * x3 * (3 * pow(x0, 2) - 2 * x0 * x2 - 2 * x0 * x3 + x2 * x3)
         / ((x1 - x2) * (x1 - x3)
            * (3 * pow(x0, 2) - 2 * x0 * x1 - 2 * x0 * x2 - 2 * x0 * x3 + x1 * x2
               + x1 * x3 + x2 * x3));
  fac2 = -x1 * x3 * (3 * pow(x0, 2) - 2 * x0 * x1 - 2 * x0 * x3 + x1 * x3)
         / ((x1 - x2) * (x2 - x3)
            * (3 * pow(x0, 2) - 2 * x0 * x1 - 2 * x0 * x2 - 2 * x0 * x3 + x1 * x2
               + x1 * x3 + x2 * x3));
  fac3 = x1 * x2 * (3 * pow(x0, 2) - 2 * x0 * x1 - 2 * x0 * x2 + x1 * x2)
         / ((x1 - x3) * (x2 - x3)
            * (3 * pow(x0, 2) - 2 * x0 * x1 - 2 * x0 * x2 - 2 * x0 * x3 + x1 * x2
               + x1 * x3 + x2 * x3));
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
    BoutReal fac0;
    BoutReal x0;
    BoutReal fac1;
    BoutReal x1;
    BoutReal fac2;
    BoutReal x2;
    BoutReal fac3;
    BoutReal x3;

    const Field2D& dy =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;
    Indices i0{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i1{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    Indices i2{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    Indices i3{bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by, 0};
    if (stagger == 0) {
      BoutReal st = 0;
      t = dy[i0];
      x0 = st + t / 2;
      st += t;
      t = dy[i1];
      x1 = st + t / 2;
      st += t;
      t = dy[i2];
      x2 = st + t / 2;
      st += t;
      t = dy[i3];
      x3 = st + t / 2;
      st += t;
    } else {
      x0 = 0;
      x1 = x0 + dy[i1];
      x2 = x1 + dy[i2];
      x3 = x2 + dy[i3];
    }

    for (int i = 0; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stagger == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        x2 += t;
        x3 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stagger, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, x2, x3, fac0, fac1, fac2, fac3);
        x0 += t;
        x1 += t;
        x2 += t;
        x3 += t;
      } else {
        t = dy[ic];
        if (stagger == -1) {
          x0 += t;
          x1 += t;
          x2 += t;
          x3 += t;
        }
        calc_interp_to_stencil(x0, x1, x2, x3, fac0, fac1, fac2, fac3);
        if (stagger == 1) {
          x0 += t;
          x1 += t;
          x2 += t;
          x3 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = ic.z;
        i2.z = ic.z;
        i3.z = ic.z;
        t = fac0 * f[i0] + fac1 * f[i1] + fac2 * f[i2] + fac3 * f[i3];

        f[ic] = t;
      }
    }
  }
}

BoundaryOp* BoundaryFreeNonUniform_O4::clone(BoundaryRegion* region,
                                             const std::list<std::string>& args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryFreeNonUniform_O4(region, newgen);
}

void BoundaryFreeNonUniform_O4::calc_interp_to_stencil(BoutReal x0, BoutReal x1,
                                                       BoutReal x2, BoutReal x3,
                                                       BoutReal& fac0, BoutReal& fac1,
                                                       BoutReal& fac2,
                                                       BoutReal& fac3) const {
  // Stencil Code
  fac0 = -x1 * x2 * x3 / ((x0 - x1) * (x0 - x2) * (x0 - x3));
  fac1 = x0 * x2 * x3 / ((x0 - x1) * (x1 - x2) * (x1 - x3));
  fac2 = -x0 * x1 * x3 / ((x0 - x2) * (x1 - x2) * (x2 - x3));
  fac3 = x0 * x1 * x2 / ((x0 - x3) * (x1 - x3) * (x2 - x3));
}
