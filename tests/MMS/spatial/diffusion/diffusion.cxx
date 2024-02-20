#include <bout/bout.hxx>
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/physicsmodel.hxx>
#include <math.h>

class Diffusion : public PhysicsModel {
  Field3D N;

  BoutReal Dx, Dy, Dz;
  BoutReal Lx, Ly, Lz;

protected:
  int init(bool UNUSED(restarting)) override {
    // Get the options
    Options* meshoptions = Options::getRoot()->getSection("mesh");

    Coordinates* coords = mesh->getCoordinates();

    meshoptions->get("Lx", Lx, 1.0);
    meshoptions->get("Ly", Ly, 1.0);

    /*this assumes equidistant grid*/
    coords->setDx(Lx / (mesh->GlobalNx - 2 * mesh->xstart));

    coords->setDy(Ly / (mesh->GlobalNy - 2 * mesh->ystart));

    output.write("SIZES: {:d}, {:d}, {:e}\n", mesh->GlobalNy,
                 (mesh->GlobalNy - 2 * mesh->ystart), coords->dy()(0, 0, 0));

    SAVE_ONCE2(Lx, Ly)

    Options* cytooptions = Options::getRoot()->getSection("cyto");
    OPTION(cytooptions, Dx, 1.0);
    OPTION(cytooptions, Dy, -1.0);
    OPTION(cytooptions, Dz, -1.0);

    SAVE_ONCE3(Dx, Dy, Dz)

    // set mesh
    auto contravariant_metric_tensor = ContravariantMetricTensor(1.1, 1.1, 1.1, 0.0, 0.0, 0.0);
    coords->setContravariantMetricTensor(contravariant_metric_tensor);

    auto covariant_metric_tensor = CovariantMetricTensor(1.1, 1.1, 1.1, 0.0, 0.0, 0.0);
    coords->setCovariantMetricTensor(covariant_metric_tensor);

    // Tell BOUT++ to solve N
    SOLVE_FOR(N)

    return 0;
  }

  int rhs(BoutReal UNUSED(t)) override {
    mesh->communicate(N); // Communicate guard cells

    ddt(N) = 0.0;

    if (Dx > 0.0) {
      ddt(N) += Dx * D2DX2(N);
    }

    if (Dy > 0.0) {
      ddt(N) += Dy * D2DY2(N);
    }

    if (Dz > 0.0) {
      ddt(N) += Dz * D2DZ2(N);
    }

    return 0;
  }
};

BOUTMAIN(Diffusion)
