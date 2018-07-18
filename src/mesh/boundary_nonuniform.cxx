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
void BoundaryDirichletNonUniform_O4::apply(Field3D &f, BoutReal t) {
  bndry->first();

  // todo: staggering

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg)
    fg = f.getBndryGenerator(bndry->location);
  fg = nullptr;

  BoutReal val = 0.0;

  Mesh *mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  BoutReal vals[mesh->LocalNz];
  int bx = bndry->bx;
  int stag = 0;
  if (loc == CELL_XLOW) {
    if (bx == 0) {
      bx = -1;
    } else if (bx < 0) {
      stag = -1;
    } else {
      stag = 1;
    }
  }
  int by = bndry->by;
  if (loc == CELL_YLOW) {
    if (by == 0) {
      by = -1;
    } else if (by < 0) {
      stag = -1;
    } else {
      stag = 1;
    }
  }
  int istart = (stag == -1) ? -1 : 0;

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

    BoutReal a, b, c, d;
    BoutReal x0, x1, x2, x3;

    const Field2D &dy =
        bndry->by != 0 ? mesh->coordinates()->dy : mesh->coordinates()->dx;

    Indices i1{bndry->x - 1 * bndry->bx, bndry->y - 1 * bndry->by, 0};
    Indices i2{bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
    Indices i3{bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
    BoutReal t;
    if (stag == 0) {
      x0 = 0; // dy(bndry->x, bndry->y) / 2;
      t = dy[i1];
      x1 = t / 2;
      x2 = t;
      t = dy[i2];
      x3 = x2 + t;
      x2 += t / 2;
      t = dy[i3];
      x3 += t / 2;
    } else {
      x0 = 0; // dy(bndry->x, bndry->y) / 2;
      x1 = dy[i1];
      x2 = x1 + dy[i2];
      x3 = x2 + dy[i3];
    }
    if (stag == -1) {
      i1 = {bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, 0};
      i2 = {bndry->x - 3 * bndry->bx, bndry->y - 3 * bndry->by, 0};
      i3 = {bndry->x - 4 * bndry->bx, bndry->y - 4 * bndry->by, 0};
    }

    for (int i = istart; i < bndry->width; i++) {
      Indices ic{bndry->x + i * bndry->bx, bndry->y + i * bndry->by, 0};
      if (stag == 0) {
        t = dy[ic] / 2;
        x0 += t;
        x1 += t;
        x2 += t;
        x3 += t;
        // printf("%+2d: %d %d %g %g %g %g\n", stag, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, x2, x3, a, b, c, d);
        x0 += t;
        x1 += t;
        x2 += t;
        x3 += t;
      } else {
        t = dy[ic];
        if (stag == -1 && i != -1) {
          x0 += t;
          x1 += t;
          x2 += t;
          x3 += t;
        }
        // printf("%+2d: %d %d %g %g %g %g\n", stag, ic.x, ic.y, x0, x1, x2, x3);
        calc_interp_to_stencil(x0, x1, x2, x3, a, b, c, d);
        if (stag == 1) {
          x0 += t;
          x1 += t;
          x2 += t;
          x3 += t;
        }
      }
      for (ic.z = 0; ic.z < mesh->LocalNz; ic.z++) {
        i1.z = i2.z = i3.z = ic.z;
        val = (fg) ? vals[ic.z] : 0.0;
        if (x0 == 0) {
          t = val;
        } else {
          t = a * val + b * f[i1] + c * f[i2] + d * f[i3];
        }
        f[ic] = t; // a * val + b * f[i1] + c * f[i2] + d * f[i3];
      }
    }
  }
}

BoundaryOp *BoundaryDirichletNonUniform_O4::clone(BoundaryRegion *region,
                                                  const list<string> &args) {
  // verifyNumPoints(region, 3);

  std::shared_ptr<FieldGenerator> newgen;
  if (!args.empty()) {
    // First argument should be an expression
    newgen = FieldFactory::get()->parse(args.front());
  }
  return new BoundaryDirichletNonUniform_O4(region, newgen);
}

void BoundaryDirichletNonUniform_O4::calc_interp_to_stencil(BoutReal a, BoutReal b,
                                                            BoutReal c, BoutReal d,
                                                            BoutReal &e, BoutReal &f,
                                                            BoutReal &g,
                                                            BoutReal &h) const {
  auto ab = a * b;
  auto cd = c * d;
  auto bma = b - a;
  auto cma = c - a;
  auto cmb = c - b;
  auto dma = d - a;
  auto bmd = b - d; // minus intentionally hidden here
  auto dmc = d - c;
  auto s1 = bma * cma * dma;
  auto s2 = bma * cmb * bmd;
  auto s3 = cma * cmb * dmc;
  auto s4 = dma * bmd * dmc;
  e = cd * b / s1;
  f = a * cd / s2;
  g = ab * d / s3;
  h = ab * c / s4;
}

void BoundaryNeumannNonUniform_O4::apply(Field3D &f, BoutReal t) { ... }

void BoundaryNeumannNonUniform_O4::calc_interp_to_stencil(BoutReal a, BoutReal b,
                                                          BoutReal c, BoutReal d,
                                                          BoutReal &e, BoutReal &f,
                                                          BoutReal &g,
                                                          BoutReal &h) const {
  // N4: matrix(
  // 	    [0,1,a,a**2/2],
  // 	    [1,b,b**2/2,b**3/6],
  // 	    [1,c,c**2/2,c**3/6],
  // 	    [1,d,d**2/2,d**3/6]
  // 	    );
  // iN4:invert(N4);
  // r:matrix(
  // 	  [1],
  // 	  [0],
  // 	  [0],
  // 	  [0]
  // 	  );
  // fortran(ratsimp(iN4.r));
  // replace a**2 with a*a
  // and move fraction to top
  BoutReal t = 1 / ((c + b - 2 * a) * d + (b - 2 * a) * c - 2 * a * b + 3 * a * a);
  e = -(b * c * d) * t;
  f = ((c + b) * d + b * c) * t;
  g = -(2 * d + 2 * c + 2 * b) * t;
  h = 6 * t;
}
