#include <bout/physicsmodel.hxx>

class TestCommunications : public PhysicsModel {
public:
  int init(bool) {
    SOLVE_FOR(f);
    f = -1.;

    // fill non-guard cells:

    // interior cells
    BOUT_FOR(i, f.getRegion("RGN_NOBNDRY")) {
      f[i] = mesh->GlobalNzNoBoundaries*(
               mesh->GlobalNyNoBoundaries*mesh->getGlobalXIndexNoBoundaries(i.x())
               + mesh->getGlobalYIndexNoBoundaries(i.y()))
             + i.z();
    }

    // lower x-boundary cells
    int startind =
      mesh->GlobalNxNoBoundaries*mesh->GlobalNyNoBoundaries*mesh->GlobalNzNoBoundaries;
    if (mesh->firstX()) {
      for (int x = 0; x < mesh->xstart; x++) {
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          for (int z = mesh->zstart; z <= mesh->zend; z++) {
            f(x, y, z) =
              startind
              + mesh->GlobalNzNoBoundaries*(mesh->GlobalNyNoBoundaries*x +
                                            mesh->getGlobalYIndexNoBoundaries(y))
              + z;
          }
        }
      }
    }
    startind += mesh->xstart*mesh->GlobalNyNoBoundaries*mesh->GlobalNzNoBoundaries;

    // upper x-boundary cells
    if (mesh->lastX()) {
      for (int x = 0; x < mesh->xstart; x++) {
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          for (int z = mesh->zstart; z <= mesh->zend; z++) {
            f(mesh->xend + 1 + x, y, z) =
              startind
              + mesh->GlobalNzNoBoundaries*(mesh->GlobalNyNoBoundaries*x +
                                            mesh->getGlobalYIndexNoBoundaries(y))
              + z;
          }
        }
      }
    }
    startind += mesh->xstart*mesh->GlobalNyNoBoundaries*mesh->GlobalNzNoBoundaries;

    // lower y-boundary cells
    for (auto it = mesh->iterateBndryLowerY(); !it.isDone(); it++) {
      int x = it.ind;
      for (int y = 0; y < mesh->ystart; y++) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          f(x, y, z) =
            startind
            + mesh->GlobalNzNoBoundaries*(mesh->getGlobalXIndex(x) + y)
            + z;
        }
      }
    }
    startind += mesh->GlobalNx*mesh->ystart*mesh->GlobalNzNoBoundaries;

    // upper y-boundary cells
    for (auto it = mesh->iterateBndryUpperY(); !it.isDone(); it++) {
      int x = it.ind;
      for (int y = 0; y < mesh->ystart; y++) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          f(x, mesh->yend + 1 + y, z) =
            startind
            + mesh->GlobalNzNoBoundaries*(mesh->getGlobalXIndex(x) + y)
            + z;
        }
      }
    }
    startind += mesh->GlobalNx*mesh->ystart*mesh->GlobalNzNoBoundaries;

    // there are no z-boundaries

    return 0;
  }

  int rhs(BoutReal) {
    // communicate f to fill guard cells
    mesh->communicate(f);
    ddt(f)=0;
    return 0;
  }
private:
  Field3D f;
};

BOUTMAIN(TestCommunications);
