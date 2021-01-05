/*
 * Demonstrates how to use staggered grids with boundary conditions
 */

#include <bout/physicsmodel.hxx>
#include <bout.hxx>
#include <derivs.hxx>

class TestStaggered : public PhysicsModel {
  Field3D n, v;
  CELL_LOC maybe_ylow{CELL_CENTRE};

protected:
  int init(bool UNUSED(restart)) override {

    if (mesh->StaggerGrids) {
      maybe_ylow = CELL_YLOW;
    }

    v.setLocation(maybe_ylow); // Staggered relative to n

    SOLVE_FOR(n, v);

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    mesh->communicate(n, v);

    ddt(n) = -n * Grad_par(v, CELL_CENTRE) - Vpar_Grad_par(v, n, CELL_CENTRE);

    ddt(v) = -Grad_par(n, maybe_ylow);

    // Have to manually apply the lower Y boundary region, using a width of 3
    for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++) {
      for (int y = 2; y >= 0; y--) {
        for (int z = 0; z < mesh->LocalNz; z++) {
          ddt(v)(rlow.ind, y, z) = ddt(v)(rlow.ind, y + 1, z);
        }
      }
    }

    return 0;
  }
};

BOUTMAIN(TestStaggered)
