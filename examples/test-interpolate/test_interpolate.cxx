/*
 * Interpolation regression test
 *
 * Test the Interpolation class
 *
 */

// for std::bind
#include <functional>
#include <random>
#include <string>

#include "bout.hxx"
#include "bout/physicsmodel.hxx"
#include "field_factory.hxx"
#include "interpolation_factory.hxx"
#include "bout/constants.hxx"

class TestInterpolate : public PhysicsModel {

  std::default_random_engine generator;
  std::uniform_real_distribution<BoutReal> distribution{0.0, 1.0};

public:

  TestInterpolate() {}

  FieldGenerator* getGeneratorFromOptions(const std::string varname, std::string &func) {
    Options *options = Options::getRoot()->getSection(varname);
    options->get("solution", func, "0.0");

    if(!func.empty()) {
      return FieldFactory::get()->parse(func);
    }
  }

  int init(bool restarting) {
    FieldFactory f(mesh);

    mesh->coordinates();

    std::string a_func;
    FieldGenerator* a_gen = getGeneratorFromOptions("a", a_func);
    Field3D a = f.create3D(a_func);
    mesh->communicate(a);
    Field3D a_solution = 0;
    Field3D a_interp = 0;

    std::string b_func;
    FieldGenerator* b_gen = getGeneratorFromOptions("b", b_func);
    Field3D b = f.create3D(b_func);
    mesh->communicate(b);
    Field3D b_solution = 0;
    Field3D b_interp = 0;

    std::string c_func;
    FieldGenerator* c_gen = getGeneratorFromOptions("c", c_func);
    Field3D c = f.create3D(c_func);
    mesh->communicate(c);
    Field3D c_solution = 0;
    Field3D c_interp = 0;

    Field3D deltax = 0;
    Field3D deltaz = 0;

    Field3D meshdx = 0;
    Field3D meshdy = 0;
    Field3D meshdz = 0;

    auto dice = std::bind(distribution, generator);

    for (const auto index : deltax) {
      BoutReal dx = index.x + dice();
      BoutReal dz = index.z + dice();
      deltax[index] = dx;
      deltaz[index] = dz;
      BoutReal x = mesh->GlobalX(dx);
      BoutReal y = TWOPI*mesh->GlobalY(index.y);
      BoutReal z = TWOPI*((BoutReal) dz) / ((BoutReal) (mesh->ngz-1));
      a_solution[index] = a_gen->generate(x, y, z, 0.0);
      b_solution[index] = b_gen->generate(x, y, z, 0.0);
      c_solution[index] = c_gen->generate(x, y, z, 0.0);
      meshdx[index] = x;
      meshdy[index] = y;
      meshdz[index] = z;
    }

    mesh->communicate(a_solution);
    mesh->communicate(b_solution);
    mesh->communicate(c_solution);

    Interpolation *interp = InterpolationFactory::getInstance()->create();

    a_interp = interp->interpolate(a, deltax, deltaz);
    b_interp = interp->interpolate(b, deltax, deltaz);
    c_interp = interp->interpolate(c, deltax, deltaz);

    SAVE_ONCE3(a, a_interp, a_solution);
    SAVE_ONCE3(b, b_interp, b_solution);
    SAVE_ONCE3(c, c_interp, c_solution);
    SAVE_ONCE2(deltax, deltaz);
    SAVE_ONCE3(meshdx, meshdy, meshdz);

    // Write data to file
    dump.write();
    dump.close();

    // Need to wait for all processes to finish writing
    MPI_Barrier(BoutComm::get());

    // Send an error code so quits
    return 1;
  }

  int rhs(BoutReal time);

};

BOUTMAIN(TestInterpolate);

int TestInterpolate::rhs(BoutReal time) {
  // Doesn't do anything
  return 1;
}
