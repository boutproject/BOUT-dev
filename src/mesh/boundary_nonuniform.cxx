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

#if !BOUT_USE_METRIC_3D
using IndMetric = Ind2D;
#define IND(var, x, y, z) IndMetric var(x* localNy + y, localNy, 1)
#define IND3D(var, x, y, z) Ind3D var((x * localNy + y) * localNz + z, localNy, localNz)
#else
using IndMetric = Ind3D;
#define IND(var, x, y, z) IndMetric var((x * localNy + y) * localNz + z, localNy, localNz)
#define IND3D(var, x, y, z)
#endif

void BoundaryDirichletNonUniform_O2::apply(Field3D& f, MAYBE_UNUSED(BoutReal t)) {
  bndry->first();
  Mesh* mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  // Decide which generator to use
  std::shared_ptr<FieldGenerator> fg = gen;
  if (!fg) {
    fg = f.getBndryGenerator(bndry->location);
  }

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

    const Coordinates::FieldMetric& coords_field =
        bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;

#if BOUT_USE_METRIC_3D
    for (int iz{0}; iz < mesh->LocalNz; iz++) {
#else
    const int iz = 0;
#endif
      const int localNy = mesh->LocalNy;
      const int localNz = mesh->LocalNz;
      const int offset = (bndry->bx * localNy + bndry->by)
#if BOUT_USE_METRIC_3D
                         * localNz
#endif
          ;
      const IND(temp, bndry->x, bndry->y, 0);
      IND3D(temp3d, bndry->x, bndry->y, 0);
      IndMetric i1{temp - 1 * offset};
#if !BOUT_USE_METRIC_3D
      Ind3D i13d{temp3d - 1 * offset * localNz};
#endif

      if (stagger == 0) {
        spacing.f0 = 0;
        BoutReal total_offset = 0;
        BoutReal offset = coords_field[i1 + iz];
        spacing.f1 = total_offset + offset / 2;
      } else {
        spacing.f0 = 0;
        spacing.f1 = spacing.f0 + coords_field[i1 + iz];
      }
      if (stagger == -1) {
#if BOUT_USE_METRIC_3D
        i1 = temp - 2 * offset;
#else
      i13d = temp3d - 2 * offset;
#endif
      }

      // with dirichlet, we specify the value on the boundary, even if
      // the value is part of the evolving system.
      for (int i = ((stagger == -1) ? -1 : 0); i < bndry->width; i++) {
        IndMetric icm{temp + IndMetric(i * offset)};
#if BOUT_USE_METRIC_3D
        icm + iz;
        IndMetric ic = icm;
#else
      Ind3D ic{temp3d + Ind3D(i * offset * localNz)};
#endif
        vec2 facs;
        if (stagger == 0) {
          BoutReal to_add = coords_field[icm] / 2;
          spacing += to_add;
          facs = calc_interp_to_stencil(spacing);
          spacing += to_add;
        } else {
          if (stagger == -1 && i != -1) {
            spacing += coords_field[icm];
          }
          facs = calc_interp_to_stencil(spacing);
          if (stagger == 1) {
            spacing += coords_field[icm];
          }
        }
#if !BOUT_USE_METRIC_3D
        for (int iz{0}; iz < mesh->LocalNz; iz++) {
#endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
          const BoutReal val = (fg) ? vals[iz] : 0.0;
          f[ic] = facs.f0 * val + facs.f1 * f[MAKE3D(i1) + iz];
        }
      }
    }
  }

  BoundaryOp* BoundaryDirichletNonUniform_O2::clone(BoundaryRegion * region,
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

  void BoundaryNeumannNonUniform_O2::apply(Field3D & f, MAYBE_UNUSED(BoutReal t)) {
    bndry->first();
    Mesh* mesh = f.getMesh();
    CELL_LOC loc = f.getLocation();

    // Decide which generator to use
    std::shared_ptr<FieldGenerator> fg = gen;
    if (!fg) {
      fg = f.getBndryGenerator(bndry->location);
    }

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

      const Coordinates::FieldMetric& coords_field =
          bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;

      const int localNy = mesh->LocalNy;
      const int localNz = mesh->LocalNz;
      const int offset = (bndry->bx * localNy + bndry->by)
#if BOUT_USE_METRIC_3D
                         * localNz
#endif
          ;
      const IND(temp, bndry->x, bndry->y, 0);
      IND3D(temp3d, bndry->x, bndry->y, 0);
      IndMetric i1{temp - 1 * offset};
#if !BOUT_USE_METRIC_3D
      Ind3D i13d{temp3d - 1 * offset * localNz};
#endif

#if BOUT_USE_METRIC_3D
      for (int iz{0}; iz < mesh->LocalNz; iz++) {
#else
    const int iz = 0;
#endif
        if (stagger == 0) {
          spacing.f0 = 0;
          BoutReal total_offset = 0;
          BoutReal offset = coords_field[i1 + iz];
          spacing.f1 = total_offset + offset / 2;
        } else { // stagger != 0
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
            spacing.f1 = spacing.f0 + coords_field[i1 + iz];
          }
        } // stagger != 0
        // With neumann (and free) the value is not set if the point is
        // evolved and it is on the boundary.
        for (int i = 0; i < bndry->width; i++) {
          IndMetric icm{temp + IndMetric(i * offset)};
#if BOUT_USE_METRIC_3D
          icm += iz;
          IndMetric ic = icm;
#else
      Ind3D ic{temp3d + Ind3D(i * offset * localNz)};
#endif
          vec2 facs;
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
#if !BOUT_USE_METRIC_3D
          for (int iz{0}; iz < mesh->LocalNz; iz++) {
#endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
            const BoutReal val = (fg) ? vals[iz] : 0.0;
            f[ic] = facs.f0 * val + facs.f1 * f[MAKE3D(i1) + iz];
          }
        }
      }
    }

    BoundaryOp* BoundaryNeumannNonUniform_O2::clone(BoundaryRegion * region,
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

    void BoundaryFreeNonUniform_O2::apply(Field3D & f, MAYBE_UNUSED(BoutReal t)) {
      bndry->first();
      Mesh* mesh = f.getMesh();
      CELL_LOC loc = f.getLocation();

      int x_boundary_offset = bndry->bx;
      int y_boundary_offset = bndry->by;
      int stagger = 0;
      update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

      for (; !bndry->isDone(); bndry->next1d()) {

        vec2 spacing;

        const Coordinates::FieldMetric& coords_field =
            bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;

        const int localNy = mesh->LocalNy;
        const int localNz = mesh->LocalNz;
        const int offset = (bndry->bx * localNy + bndry->by)
#if BOUT_USE_METRIC_3D
                           * localNz
#endif
            ;
        const IND(temp, bndry->x, bndry->y, 0);
        IND3D(temp3d, bndry->x, bndry->y, 0);
        const IndMetric i0{temp - 1 * offset};
#if !BOUT_USE_METRIC_3D
        const Ind3D i03d{temp3d - 1 * offset * localNz};
#endif
        const IndMetric i1{temp - 2 * offset};
#if !BOUT_USE_METRIC_3D
        const Ind3D i13d{temp3d - 2 * offset * localNz};
#endif

#if BOUT_USE_METRIC_3D
        for (int iz{0}; iz < mesh->LocalNz; iz++) {
#else
    const int iz = 0;
#endif
          if (stagger == 0) {
            BoutReal total_offset = 0;
            BoutReal offset = coords_field[i0 + iz];
            spacing.f0 = total_offset + offset / 2;
            total_offset += offset;
            offset = coords_field[i1 + iz];
            spacing.f1 = total_offset + offset / 2;
          } else {
            spacing.f0 = coords_field[i0 + iz];
            spacing.f1 = spacing.f0 + coords_field[i1 + iz];
          }

          // With free (and neumann) the value is not set if the point is
          // evolved and it is on the boundary.
          for (int i = 0; i < bndry->width; i++) {
            IndMetric icm{temp + IndMetric(i * offset)};
#if BOUT_USE_METRIC_3D
            icm += iz;
            IndMetric ic = icm;
#else
      Ind3D ic{temp3d + Ind3D(i * offset * localNz)};
#endif
            vec2 facs;
            if (stagger == 0) {
              BoutReal to_add = coords_field[icm] / 2;
              spacing += to_add;
              facs = calc_interp_to_stencil(spacing);
              spacing += to_add;
            } else {
              facs = calc_interp_to_stencil(spacing);
              spacing += coords_field[icm];
            }
#if !BOUT_USE_METRIC_3D
            for (int iz{0}; iz < mesh->LocalNz; iz++) {
#endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
              const BoutReal val = f[MAKE3D(i0) + iz];
              f[ic] = facs.f0 * val + facs.f1 * f[MAKE3D(i1) + iz];
            }
          }
        }
      }

      BoundaryOp* BoundaryFreeNonUniform_O2::clone(BoundaryRegion * region,
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

      void BoundaryDirichletNonUniform_O3::apply(Field3D & f, MAYBE_UNUSED(BoutReal t)) {
        bndry->first();
        Mesh* mesh = f.getMesh();
        CELL_LOC loc = f.getLocation();

        // Decide which generator to use
        std::shared_ptr<FieldGenerator> fg = gen;
        if (!fg) {
          fg = f.getBndryGenerator(bndry->location);
        }

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
              vals[zk] = fg->generate(bout::generator::Context().set(
                  "x", xnorm, "y", ynorm, "z", zfac * zk, "t", t));
            }
          }

          vec3 spacing;

          const Coordinates::FieldMetric& coords_field =
              bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;

#if BOUT_USE_METRIC_3D
          for (int iz{0}; iz < mesh->LocalNz; iz++) {
#else
    const int iz = 0;
#endif
            const int localNy = mesh->LocalNy;
            const int localNz = mesh->LocalNz;
            const int offset = (bndry->bx * localNy + bndry->by)
#if BOUT_USE_METRIC_3D
                               * localNz
#endif
                ;
            const IND(temp, bndry->x, bndry->y, 0);
            IND3D(temp3d, bndry->x, bndry->y, 0);
            IndMetric i1{temp - 1 * offset};
#if !BOUT_USE_METRIC_3D
            Ind3D i13d{temp3d - 1 * offset * localNz};
#endif
            IndMetric i2{temp - 2 * offset};
#if !BOUT_USE_METRIC_3D
            Ind3D i23d{temp3d - 2 * offset * localNz};
#endif

            if (stagger == 0) {
              spacing.f0 = 0;
              BoutReal total_offset = 0;
              BoutReal offset = coords_field[i1 + iz];
              spacing.f1 = total_offset + offset / 2;
              total_offset += offset;
              offset = coords_field[i2 + iz];
              spacing.f2 = total_offset + offset / 2;
            } else {
              spacing.f0 = 0;
              spacing.f1 = spacing.f0 + coords_field[i1 + iz];
              spacing.f2 = spacing.f1 + coords_field[i2 + iz];
            }
            if (stagger == -1) {
#if BOUT_USE_METRIC_3D
              i1 = temp - 2 * offset;
#else
      i13d = temp3d - 2 * offset;
#endif
#if BOUT_USE_METRIC_3D
              i2 = temp - 3 * offset;
#else
      i23d = temp3d - 3 * offset;
#endif
            }

            // with dirichlet, we specify the value on the boundary, even if
            // the value is part of the evolving system.
            for (int i = ((stagger == -1) ? -1 : 0); i < bndry->width; i++) {
              IndMetric icm{temp + IndMetric(i * offset)};
#if BOUT_USE_METRIC_3D
              icm + iz;
              IndMetric ic = icm;
#else
      Ind3D ic{temp3d + Ind3D(i * offset * localNz)};
#endif
              vec3 facs;
              if (stagger == 0) {
                BoutReal to_add = coords_field[icm] / 2;
                spacing += to_add;
                facs = calc_interp_to_stencil(spacing);
                spacing += to_add;
              } else {
                if (stagger == -1 && i != -1) {
                  spacing += coords_field[icm];
                }
                facs = calc_interp_to_stencil(spacing);
                if (stagger == 1) {
                  spacing += coords_field[icm];
                }
              }
#if !BOUT_USE_METRIC_3D
              for (int iz{0}; iz < mesh->LocalNz; iz++) {
#endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
                const BoutReal val = (fg) ? vals[iz] : 0.0;
                f[ic] = facs.f0 * val + facs.f1 * f[MAKE3D(i1) + iz]
                        + facs.f2 * f[MAKE3D(i2) + iz];
              }
            }
          }
        }

        BoundaryOp* BoundaryDirichletNonUniform_O3::clone(
            BoundaryRegion * region, const std::list<std::string>& args) {

          std::shared_ptr<FieldGenerator> newgen;
          if (!args.empty()) {
            // First argument should be an expression
            newgen = FieldFactory::get()->parse(args.front());
          }
          return new BoundaryDirichletNonUniform_O3(region, newgen);
        }

        vec3 BoundaryDirichletNonUniform_O3::calc_interp_to_stencil(const vec3& spacing)
            const {
          vec3 facs;
          // Stencil Code
          facs.f0 = spacing.f1 * spacing.f2
                    / ((spacing.f0 - spacing.f1) * (spacing.f0 - spacing.f2));
          facs.f1 = -spacing.f0 * spacing.f2
                    / ((spacing.f0 - spacing.f1) * (spacing.f1 - spacing.f2));
          facs.f2 = spacing.f0 * spacing.f1
                    / ((spacing.f0 - spacing.f2) * (spacing.f1 - spacing.f2));

          return facs;
        }

        void BoundaryNeumannNonUniform_O3::apply(Field3D & f, MAYBE_UNUSED(BoutReal t)) {
          bndry->first();
          Mesh* mesh = f.getMesh();
          CELL_LOC loc = f.getLocation();

          // Decide which generator to use
          std::shared_ptr<FieldGenerator> fg = gen;
          if (!fg) {
            fg = f.getBndryGenerator(bndry->location);
          }

          std::vector<BoutReal> vals;
          vals.reserve(mesh->LocalNz);

          int x_boundary_offset = bndry->bx;
          int y_boundary_offset = bndry->by;
          int stagger = 0;
          update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

          for (; !bndry->isDone(); bndry->next1d()) {
            if (fg) {
              // Calculate the X and Y normalised values half-way between the guard cell
              // and grid cell
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
                vals[zk] = fg->generate(bout::generator::Context().set(
                    "x", xnorm, "y", ynorm, "z", zfac * zk, "t", t));
              }
            }

            vec3 spacing;

            const Coordinates::FieldMetric& coords_field =
                bndry->by != 0 ? mesh->getCoordinates()->dy : mesh->getCoordinates()->dx;

            const int localNy = mesh->LocalNy;
            const int localNz = mesh->LocalNz;
            const int offset = (bndry->bx * localNy + bndry->by)
#if BOUT_USE_METRIC_3D
                               * localNz
#endif
                ;
            const IND(temp, bndry->x, bndry->y, 0);
            IND3D(temp3d, bndry->x, bndry->y, 0);
            IndMetric i1{temp - 1 * offset};
#if !BOUT_USE_METRIC_3D
            Ind3D i13d{temp3d - 1 * offset * localNz};
#endif
            IndMetric i2{temp - 2 * offset};
#if !BOUT_USE_METRIC_3D
            Ind3D i23d{temp3d - 2 * offset * localNz};
#endif

#if BOUT_USE_METRIC_3D
            for (int iz{0}; iz < mesh->LocalNz; iz++) {
#else
    const int iz = 0;
#endif
              if (stagger == 0) {
                spacing.f0 = 0;
                BoutReal total_offset = 0;
                BoutReal offset = coords_field[i1 + iz];
                spacing.f1 = total_offset + offset / 2;
                total_offset += offset;
                offset = coords_field[i2 + iz];
                spacing.f2 = total_offset + offset / 2;
              } else { // stagger != 0
                spacing.f0 = 0;
                // Check if we are staggered and also boundary in low
                //  direction
                // In the case of Neumann we have in this case two values
                //  defined at the same point
                if (stagger == -1
                    && ((bndry->bx && x_boundary_offset == -1)
                        || (bndry->by && y_boundary_offset == -1))) {
                  spacing.f1 = spacing.f0;
                  spacing.f2 = spacing.f1 + coords_field[i1 + iz];
                } else {
                  spacing.f1 = spacing.f0 + coords_field[i1 + iz];
                  spacing.f2 = spacing.f1 + coords_field[i2 + iz];
                }
              } // stagger != 0
              // With neumann (and free) the value is not set if the point is
              // evolved and it is on the boundary.
              for (int i = 0; i < bndry->width; i++) {
                IndMetric icm{temp + IndMetric(i * offset)};
#if BOUT_USE_METRIC_3D
                icm += iz;
                IndMetric ic = icm;
#else
      Ind3D ic{temp3d + Ind3D(i * offset * localNz)};
#endif
                vec3 facs;
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
#if !BOUT_USE_METRIC_3D
                for (int iz{0}; iz < mesh->LocalNz; iz++) {
#endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
                  const BoutReal val = (fg) ? vals[iz] : 0.0;
                  f[ic] = facs.f0 * val + facs.f1 * f[MAKE3D(i1) + iz]
                          + facs.f2 * f[MAKE3D(i2) + iz];
                }
              }
            }
          }

          BoundaryOp* BoundaryNeumannNonUniform_O3::clone(
              BoundaryRegion * region, const std::list<std::string>& args) {

            std::shared_ptr<FieldGenerator> newgen;
            if (!args.empty()) {
              // First argument should be an expression
              newgen = FieldFactory::get()->parse(args.front());
            }
            return new BoundaryNeumannNonUniform_O3(region, newgen);
          }

          vec3 BoundaryNeumannNonUniform_O3::calc_interp_to_stencil(const vec3& spacing)
              const {
            vec3 facs;
            // Stencil Code
            facs.f0 =
                spacing.f1 * spacing.f2 / (2 * spacing.f0 - spacing.f1 - spacing.f2);
            facs.f1 = -spacing.f2 * (2 * spacing.f0 - spacing.f2)
                      / ((spacing.f1 - spacing.f2)
                         * (2 * spacing.f0 - spacing.f1 - spacing.f2));
            facs.f2 = spacing.f1 * (2 * spacing.f0 - spacing.f1)
                      / ((spacing.f1 - spacing.f2)
                         * (2 * spacing.f0 - spacing.f1 - spacing.f2));

            return facs;
          }

          void BoundaryFreeNonUniform_O3::apply(Field3D & f, MAYBE_UNUSED(BoutReal t)) {
            bndry->first();
            Mesh* mesh = f.getMesh();
            CELL_LOC loc = f.getLocation();

            int x_boundary_offset = bndry->bx;
            int y_boundary_offset = bndry->by;
            int stagger = 0;
            update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

            for (; !bndry->isDone(); bndry->next1d()) {

              vec3 spacing;

              const Coordinates::FieldMetric& coords_field =
                  bndry->by != 0 ? mesh->getCoordinates()->dy
                                 : mesh->getCoordinates()->dx;

              const int localNy = mesh->LocalNy;
              const int localNz = mesh->LocalNz;
              const int offset = (bndry->bx * localNy + bndry->by)
#if BOUT_USE_METRIC_3D
                                 * localNz
#endif
                  ;
              const IND(temp, bndry->x, bndry->y, 0);
              IND3D(temp3d, bndry->x, bndry->y, 0);
              const IndMetric i0{temp - 1 * offset};
#if !BOUT_USE_METRIC_3D
              const Ind3D i03d{temp3d - 1 * offset * localNz};
#endif
              const IndMetric i1{temp - 2 * offset};
#if !BOUT_USE_METRIC_3D
              const Ind3D i13d{temp3d - 2 * offset * localNz};
#endif
              const IndMetric i2{temp - 3 * offset};
#if !BOUT_USE_METRIC_3D
              const Ind3D i23d{temp3d - 3 * offset * localNz};
#endif

#if BOUT_USE_METRIC_3D
              for (int iz{0}; iz < mesh->LocalNz; iz++) {
#else
    const int iz = 0;
#endif
                if (stagger == 0) {
                  BoutReal total_offset = 0;
                  BoutReal offset = coords_field[i0 + iz];
                  spacing.f0 = total_offset + offset / 2;
                  total_offset += offset;
                  offset = coords_field[i1 + iz];
                  spacing.f1 = total_offset + offset / 2;
                  total_offset += offset;
                  offset = coords_field[i2 + iz];
                  spacing.f2 = total_offset + offset / 2;
                } else {
                  spacing.f0 = coords_field[i0 + iz];
                  spacing.f1 = spacing.f0 + coords_field[i1 + iz];
                  spacing.f2 = spacing.f1 + coords_field[i2 + iz];
                }

                // With free (and neumann) the value is not set if the point is
                // evolved and it is on the boundary.
                for (int i = 0; i < bndry->width; i++) {
                  IndMetric icm{temp + IndMetric(i * offset)};
#if BOUT_USE_METRIC_3D
                  icm += iz;
                  IndMetric ic = icm;
#else
      Ind3D ic{temp3d + Ind3D(i * offset * localNz)};
#endif
                  vec3 facs;
                  if (stagger == 0) {
                    BoutReal to_add = coords_field[icm] / 2;
                    spacing += to_add;
                    facs = calc_interp_to_stencil(spacing);
                    spacing += to_add;
                  } else {
                    facs = calc_interp_to_stencil(spacing);
                    spacing += coords_field[icm];
                  }
#if !BOUT_USE_METRIC_3D
                  for (int iz{0}; iz < mesh->LocalNz; iz++) {
#endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
                    const BoutReal val = f[MAKE3D(i0) + iz];
                    f[ic] = facs.f0 * val + facs.f1 * f[MAKE3D(i1) + iz]
                            + facs.f2 * f[MAKE3D(i2) + iz];
                  }
                }
              }
            }

            BoundaryOp* BoundaryFreeNonUniform_O3::clone(
                BoundaryRegion * region, const std::list<std::string>& args) {

              std::shared_ptr<FieldGenerator> newgen;
              if (!args.empty()) {
                // First argument should be an expression
                newgen = FieldFactory::get()->parse(args.front());
              }
              return new BoundaryFreeNonUniform_O3(region, newgen);
            }

            vec3 BoundaryFreeNonUniform_O3::calc_interp_to_stencil(const vec3& spacing)
                const {
              vec3 facs;
              // Stencil Code
              facs.f0 = spacing.f1 * spacing.f2
                        / ((spacing.f0 - spacing.f1) * (spacing.f0 - spacing.f2));
              facs.f1 = -spacing.f0 * spacing.f2
                        / ((spacing.f0 - spacing.f1) * (spacing.f1 - spacing.f2));
              facs.f2 = spacing.f0 * spacing.f1
                        / ((spacing.f0 - spacing.f2) * (spacing.f1 - spacing.f2));

              return facs;
            }

            void BoundaryDirichletNonUniform_O4::apply(Field3D & f,
                                                       MAYBE_UNUSED(BoutReal t)) {
              bndry->first();
              Mesh* mesh = f.getMesh();
              CELL_LOC loc = f.getLocation();

              // Decide which generator to use
              std::shared_ptr<FieldGenerator> fg = gen;
              if (!fg) {
                fg = f.getBndryGenerator(bndry->location);
              }

              std::vector<BoutReal> vals;
              vals.reserve(mesh->LocalNz);

              int x_boundary_offset = bndry->bx;
              int y_boundary_offset = bndry->by;
              int stagger = 0;
              update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger, loc);

              for (; !bndry->isDone(); bndry->next1d()) {
                if (fg) {
                  // Calculate the X and Y normalised values half-way between the guard
                  // cell and grid cell
                  const BoutReal xnorm =
                      0.5
                      * (mesh->GlobalX(bndry->x) // In the guard cell
                         + mesh->GlobalX(bndry->x - x_boundary_offset)); // the grid cell
                  const BoutReal ynorm =
                      TWOPI * 0.5
                      * (mesh->GlobalY(bndry->y) // In the guard cell
                         + mesh->GlobalY(bndry->y - y_boundary_offset)); // the grid cell
                  const BoutReal zfac = TWOPI / mesh->LocalNz;
                  for (int zk = 0; zk < mesh->LocalNz; zk++) {
                    vals[zk] = fg->generate(bout::generator::Context().set(
                        "x", xnorm, "y", ynorm, "z", zfac * zk, "t", t));
                  }
                }

                vec4 spacing;

                const Coordinates::FieldMetric& coords_field =
                    bndry->by != 0 ? mesh->getCoordinates()->dy
                                   : mesh->getCoordinates()->dx;

#if BOUT_USE_METRIC_3D
                for (int iz{0}; iz < mesh->LocalNz; iz++) {
#else
    const int iz = 0;
#endif
                  const int localNy = mesh->LocalNy;
                  const int localNz = mesh->LocalNz;
                  const int offset = (bndry->bx * localNy + bndry->by)
#if BOUT_USE_METRIC_3D
                                     * localNz
#endif
                      ;
                  const IND(temp, bndry->x, bndry->y, 0);
                  IND3D(temp3d, bndry->x, bndry->y, 0);
                  IndMetric i1{temp - 1 * offset};
#if !BOUT_USE_METRIC_3D
                  Ind3D i13d{temp3d - 1 * offset * localNz};
#endif
                  IndMetric i2{temp - 2 * offset};
#if !BOUT_USE_METRIC_3D
                  Ind3D i23d{temp3d - 2 * offset * localNz};
#endif
                  IndMetric i3{temp - 3 * offset};
#if !BOUT_USE_METRIC_3D
                  Ind3D i33d{temp3d - 3 * offset * localNz};
#endif

                  if (stagger == 0) {
                    spacing.f0 = 0;
                    BoutReal total_offset = 0;
                    BoutReal offset = coords_field[i1 + iz];
                    spacing.f1 = total_offset + offset / 2;
                    total_offset += offset;
                    offset = coords_field[i2 + iz];
                    spacing.f2 = total_offset + offset / 2;
                    total_offset += offset;
                    offset = coords_field[i3 + iz];
                    spacing.f3 = total_offset + offset / 2;
                  } else {
                    spacing.f0 = 0;
                    spacing.f1 = spacing.f0 + coords_field[i1 + iz];
                    spacing.f2 = spacing.f1 + coords_field[i2 + iz];
                    spacing.f3 = spacing.f2 + coords_field[i3 + iz];
                  }
                  if (stagger == -1) {
#if BOUT_USE_METRIC_3D
                    i1 = temp - 2 * offset;
#else
      i13d = temp3d - 2 * offset;
#endif
#if BOUT_USE_METRIC_3D
                    i2 = temp - 3 * offset;
#else
      i23d = temp3d - 3 * offset;
#endif
#if BOUT_USE_METRIC_3D
                    i3 = temp - 4 * offset;
#else
      i33d = temp3d - 4 * offset;
#endif
                  }

                  // with dirichlet, we specify the value on the boundary, even if
                  // the value is part of the evolving system.
                  for (int i = ((stagger == -1) ? -1 : 0); i < bndry->width; i++) {
                    IndMetric icm{temp + IndMetric(i * offset)};
#if BOUT_USE_METRIC_3D
                    icm + iz;
                    IndMetric ic = icm;
#else
      Ind3D ic{temp3d + Ind3D(i * offset * localNz)};
#endif
                    vec4 facs;
                    if (stagger == 0) {
                      BoutReal to_add = coords_field[icm] / 2;
                      spacing += to_add;
                      facs = calc_interp_to_stencil(spacing);
                      spacing += to_add;
                    } else {
                      if (stagger == -1 && i != -1) {
                        spacing += coords_field[icm];
                      }
                      facs = calc_interp_to_stencil(spacing);
                      if (stagger == 1) {
                        spacing += coords_field[icm];
                      }
                    }
#if !BOUT_USE_METRIC_3D
                    for (int iz{0}; iz < mesh->LocalNz; iz++) {
#endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
                      const BoutReal val = (fg) ? vals[iz] : 0.0;
                      f[ic] = facs.f0 * val + facs.f1 * f[MAKE3D(i1) + iz]
                              + facs.f2 * f[MAKE3D(i2) + iz]
                              + facs.f3 * f[MAKE3D(i3) + iz];
                    }
                  }
                }
              }

              BoundaryOp* BoundaryDirichletNonUniform_O4::clone(
                  BoundaryRegion * region, const std::list<std::string>& args) {

                std::shared_ptr<FieldGenerator> newgen;
                if (!args.empty()) {
                  // First argument should be an expression
                  newgen = FieldFactory::get()->parse(args.front());
                }
                return new BoundaryDirichletNonUniform_O4(region, newgen);
              }

              vec4 BoundaryDirichletNonUniform_O4::calc_interp_to_stencil(
                  const vec4& spacing) const {
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

              void BoundaryNeumannNonUniform_O4::apply(Field3D & f,
                                                       MAYBE_UNUSED(BoutReal t)) {
                bndry->first();
                Mesh* mesh = f.getMesh();
                CELL_LOC loc = f.getLocation();

                // Decide which generator to use
                std::shared_ptr<FieldGenerator> fg = gen;
                if (!fg) {
                  fg = f.getBndryGenerator(bndry->location);
                }

                std::vector<BoutReal> vals;
                vals.reserve(mesh->LocalNz);

                int x_boundary_offset = bndry->bx;
                int y_boundary_offset = bndry->by;
                int stagger = 0;
                update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger,
                                       loc);

                for (; !bndry->isDone(); bndry->next1d()) {
                  if (fg) {
                    // Calculate the X and Y normalised values half-way between the guard
                    // cell and grid cell
                    const BoutReal xnorm =
                        0.5
                        * (mesh->GlobalX(bndry->x) // In the guard cell
                           + mesh->GlobalX(bndry->x
                                           - x_boundary_offset)); // the grid cell
                    const BoutReal ynorm =
                        TWOPI * 0.5
                        * (mesh->GlobalY(bndry->y) // In the guard cell
                           + mesh->GlobalY(bndry->y
                                           - y_boundary_offset)); // the grid cell
                    const BoutReal zfac = TWOPI / mesh->LocalNz;
                    for (int zk = 0; zk < mesh->LocalNz; zk++) {
                      vals[zk] = fg->generate(bout::generator::Context().set(
                          "x", xnorm, "y", ynorm, "z", zfac * zk, "t", t));
                    }
                  }

                  vec4 spacing;

                  const Coordinates::FieldMetric& coords_field =
                      bndry->by != 0 ? mesh->getCoordinates()->dy
                                     : mesh->getCoordinates()->dx;

                  const int localNy = mesh->LocalNy;
                  const int localNz = mesh->LocalNz;
                  const int offset = (bndry->bx * localNy + bndry->by)
#if BOUT_USE_METRIC_3D
                                     * localNz
#endif
                      ;
                  const IND(temp, bndry->x, bndry->y, 0);
                  IND3D(temp3d, bndry->x, bndry->y, 0);
                  IndMetric i1{temp - 1 * offset};
#if !BOUT_USE_METRIC_3D
                  Ind3D i13d{temp3d - 1 * offset * localNz};
#endif
                  IndMetric i2{temp - 2 * offset};
#if !BOUT_USE_METRIC_3D
                  Ind3D i23d{temp3d - 2 * offset * localNz};
#endif
                  IndMetric i3{temp - 3 * offset};
#if !BOUT_USE_METRIC_3D
                  Ind3D i33d{temp3d - 3 * offset * localNz};
#endif

#if BOUT_USE_METRIC_3D
                  for (int iz{0}; iz < mesh->LocalNz; iz++) {
#else
    const int iz = 0;
#endif
                    if (stagger == 0) {
                      spacing.f0 = 0;
                      BoutReal total_offset = 0;
                      BoutReal offset = coords_field[i1 + iz];
                      spacing.f1 = total_offset + offset / 2;
                      total_offset += offset;
                      offset = coords_field[i2 + iz];
                      spacing.f2 = total_offset + offset / 2;
                      total_offset += offset;
                      offset = coords_field[i3 + iz];
                      spacing.f3 = total_offset + offset / 2;
                    } else { // stagger != 0
                      spacing.f0 = 0;
                      // Check if we are staggered and also boundary in low
                      //  direction
                      // In the case of Neumann we have in this case two values
                      //  defined at the same point
                      if (stagger == -1
                          && ((bndry->bx && x_boundary_offset == -1)
                              || (bndry->by && y_boundary_offset == -1))) {
                        spacing.f1 = spacing.f0;
                        spacing.f2 = spacing.f1 + coords_field[i1 + iz];
                        spacing.f3 = spacing.f2 + coords_field[i2 + iz];
                      } else {
                        spacing.f1 = spacing.f0 + coords_field[i1 + iz];
                        spacing.f2 = spacing.f1 + coords_field[i2 + iz];
                        spacing.f3 = spacing.f2 + coords_field[i3 + iz];
                      }
                    } // stagger != 0
                    // With neumann (and free) the value is not set if the point is
                    // evolved and it is on the boundary.
                    for (int i = 0; i < bndry->width; i++) {
                      IndMetric icm{temp + IndMetric(i * offset)};
#if BOUT_USE_METRIC_3D
                      icm += iz;
                      IndMetric ic = icm;
#else
      Ind3D ic{temp3d + Ind3D(i * offset * localNz)};
#endif
                      vec4 facs;
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
#if !BOUT_USE_METRIC_3D
                      for (int iz{0}; iz < mesh->LocalNz; iz++) {
#endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
                        const BoutReal val = (fg) ? vals[iz] : 0.0;
                        f[ic] = facs.f0 * val + facs.f1 * f[MAKE3D(i1) + iz]
                                + facs.f2 * f[MAKE3D(i2) + iz]
                                + facs.f3 * f[MAKE3D(i3) + iz];
                      }
                    }
                  }
                }

                BoundaryOp* BoundaryNeumannNonUniform_O4::clone(
                    BoundaryRegion * region, const std::list<std::string>& args) {

                  std::shared_ptr<FieldGenerator> newgen;
                  if (!args.empty()) {
                    // First argument should be an expression
                    newgen = FieldFactory::get()->parse(args.front());
                  }
                  return new BoundaryNeumannNonUniform_O4(region, newgen);
                }

                vec4 BoundaryNeumannNonUniform_O4::calc_interp_to_stencil(
                    const vec4& spacing) const {
                  vec4 facs;
                  // Stencil Code
                  facs.f0 = -spacing.f1 * spacing.f2 * spacing.f3
                            / (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
                               - 2 * spacing.f0 * spacing.f2 - 2 * spacing.f0 * spacing.f3
                               + spacing.f1 * spacing.f2 + spacing.f1 * spacing.f3
                               + spacing.f2 * spacing.f3);
                  facs.f1 = spacing.f2 * spacing.f3
                            * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f2
                               - 2 * spacing.f0 * spacing.f3 + spacing.f2 * spacing.f3)
                            / ((spacing.f1 - spacing.f2) * (spacing.f1 - spacing.f3)
                               * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
                                  - 2 * spacing.f0 * spacing.f2
                                  - 2 * spacing.f0 * spacing.f3 + spacing.f1 * spacing.f2
                                  + spacing.f1 * spacing.f3 + spacing.f2 * spacing.f3));
                  facs.f2 = -spacing.f1 * spacing.f3
                            * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
                               - 2 * spacing.f0 * spacing.f3 + spacing.f1 * spacing.f3)
                            / ((spacing.f1 - spacing.f2) * (spacing.f2 - spacing.f3)
                               * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
                                  - 2 * spacing.f0 * spacing.f2
                                  - 2 * spacing.f0 * spacing.f3 + spacing.f1 * spacing.f2
                                  + spacing.f1 * spacing.f3 + spacing.f2 * spacing.f3));
                  facs.f3 = spacing.f1 * spacing.f2
                            * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
                               - 2 * spacing.f0 * spacing.f2 + spacing.f1 * spacing.f2)
                            / ((spacing.f1 - spacing.f3) * (spacing.f2 - spacing.f3)
                               * (3 * pow(spacing.f0, 2) - 2 * spacing.f0 * spacing.f1
                                  - 2 * spacing.f0 * spacing.f2
                                  - 2 * spacing.f0 * spacing.f3 + spacing.f1 * spacing.f2
                                  + spacing.f1 * spacing.f3 + spacing.f2 * spacing.f3));

                  return facs;
                }

                void BoundaryFreeNonUniform_O4::apply(Field3D & f,
                                                      MAYBE_UNUSED(BoutReal t)) {
                  bndry->first();
                  Mesh* mesh = f.getMesh();
                  CELL_LOC loc = f.getLocation();

                  int x_boundary_offset = bndry->bx;
                  int y_boundary_offset = bndry->by;
                  int stagger = 0;
                  update_stagger_offsets(x_boundary_offset, y_boundary_offset, stagger,
                                         loc);

                  for (; !bndry->isDone(); bndry->next1d()) {

                    vec4 spacing;

                    const Coordinates::FieldMetric& coords_field =
                        bndry->by != 0 ? mesh->getCoordinates()->dy
                                       : mesh->getCoordinates()->dx;

                    const int localNy = mesh->LocalNy;
                    const int localNz = mesh->LocalNz;
                    const int offset = (bndry->bx * localNy + bndry->by)
#if BOUT_USE_METRIC_3D
                                       * localNz
#endif
                        ;
                    const IND(temp, bndry->x, bndry->y, 0);
                    IND3D(temp3d, bndry->x, bndry->y, 0);
                    const IndMetric i0{temp - 1 * offset};
#if !BOUT_USE_METRIC_3D
                    const Ind3D i03d{temp3d - 1 * offset * localNz};
#endif
                    const IndMetric i1{temp - 2 * offset};
#if !BOUT_USE_METRIC_3D
                    const Ind3D i13d{temp3d - 2 * offset * localNz};
#endif
                    const IndMetric i2{temp - 3 * offset};
#if !BOUT_USE_METRIC_3D
                    const Ind3D i23d{temp3d - 3 * offset * localNz};
#endif
                    const IndMetric i3{temp - 4 * offset};
#if !BOUT_USE_METRIC_3D
                    const Ind3D i33d{temp3d - 4 * offset * localNz};
#endif

#if BOUT_USE_METRIC_3D
                    for (int iz{0}; iz < mesh->LocalNz; iz++) {
#else
    const int iz = 0;
#endif
                      if (stagger == 0) {
                        BoutReal total_offset = 0;
                        BoutReal offset = coords_field[i0 + iz];
                        spacing.f0 = total_offset + offset / 2;
                        total_offset += offset;
                        offset = coords_field[i1 + iz];
                        spacing.f1 = total_offset + offset / 2;
                        total_offset += offset;
                        offset = coords_field[i2 + iz];
                        spacing.f2 = total_offset + offset / 2;
                        total_offset += offset;
                        offset = coords_field[i3 + iz];
                        spacing.f3 = total_offset + offset / 2;
                      } else {
                        spacing.f0 = coords_field[i0 + iz];
                        spacing.f1 = spacing.f0 + coords_field[i1 + iz];
                        spacing.f2 = spacing.f1 + coords_field[i2 + iz];
                        spacing.f3 = spacing.f2 + coords_field[i3 + iz];
                      }

                      // With free (and neumann) the value is not set if the point is
                      // evolved and it is on the boundary.
                      for (int i = 0; i < bndry->width; i++) {
                        IndMetric icm{temp + IndMetric(i * offset)};
#if BOUT_USE_METRIC_3D
                        icm += iz;
                        IndMetric ic = icm;
#else
      Ind3D ic{temp3d + Ind3D(i * offset * localNz)};
#endif
                        vec4 facs;
                        if (stagger == 0) {
                          BoutReal to_add = coords_field[icm] / 2;
                          spacing += to_add;
                          facs = calc_interp_to_stencil(spacing);
                          spacing += to_add;
                        } else {
                          facs = calc_interp_to_stencil(spacing);
                          spacing += coords_field[icm];
                        }
#if !BOUT_USE_METRIC_3D
                        for (int iz{0}; iz < mesh->LocalNz; iz++) {
#endif
#if BOUT_USE_METRIC_3D
#define MAKE3D(x) x
#else
#define MAKE3D(x) x##3d
#endif
                          const BoutReal val = f[MAKE3D(i0) + iz];
                          f[ic] = facs.f0 * val + facs.f1 * f[MAKE3D(i1) + iz]
                                  + facs.f2 * f[MAKE3D(i2) + iz]
                                  + facs.f3 * f[MAKE3D(i3) + iz];
                        }
                      }
                    }
                  }

                  BoundaryOp* BoundaryFreeNonUniform_O4::clone(
                      BoundaryRegion * region, const std::list<std::string>& args) {

                    std::shared_ptr<FieldGenerator> newgen;
                    if (!args.empty()) {
                      // First argument should be an expression
                      newgen = FieldFactory::get()->parse(args.front());
                    }
                    return new BoundaryFreeNonUniform_O4(region, newgen);
                  }

                  vec4 BoundaryFreeNonUniform_O4::calc_interp_to_stencil(
                      const vec4& spacing) const {
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
