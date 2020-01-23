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
#include "bout/constants.hxx"
#include "field_factory.hxx"
#include "bout/sys/generator_context.hxx"
#include "interpolation_xz.hxx"

/// Get a FieldGenerator from the options for a variable
std::shared_ptr<FieldGenerator> getGeneratorFromOptions(const std::string& varname,
                                                        std::string& func) {
  Options *options = Options::getRoot()->getSection(varname);
  options->get("solution", func, "0.0");

  if (func.empty()) {
    throw BoutException("Couldn't read 'solution' option");
  }
  return FieldFactory::get()->parse(func);
}

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  // Random number generator
  std::default_random_engine generator;
  // Uniform distribution of BoutReals from 0 to 1
  std::uniform_real_distribution<BoutReal> distribution{0.0, 1.0};

  FieldFactory f(mesh);

  // Set up generators and solutions for three different analtyic functions
  std::string a_func;
  auto a_gen = getGeneratorFromOptions("a", a_func);
  Field3D a = f.create3D(a_func);
  Field3D a_solution = 0.0;
  Field3D a_interp = 0.0;

  std::string b_func;
  auto b_gen = getGeneratorFromOptions("b", b_func);
  Field3D b = f.create3D(b_func);
  Field3D b_solution = 0.0;
  Field3D b_interp = 0.0;

  std::string c_func;
  auto c_gen = getGeneratorFromOptions("c", c_func);
  Field3D c = f.create3D(c_func);
  Field3D c_solution = 0.0;
  Field3D c_interp = 0.0;

  // x and z displacements
  Field3D deltax = 0.0;
  Field3D deltaz = 0.0;

  // Bind the random number generator and distribution into a single function
  auto dice = std::bind(distribution, generator);

  for (const auto &index : deltax) {
    // Get some random displacements
    BoutReal dx = index.x() + dice();
    BoutReal dz = index.z() + dice();
    // For the last point, put the displacement inwards
    // Otherwise we try to interpolate in the guard cells, which doesn't work so well
    if (index.x() >= mesh->xend) {
      dx = index.x() - dice();
    }
    deltax[index] = dx;
    deltaz[index] = dz;
    // Get the global indices
    bout::generator::Context pos{index, CELL_CENTRE, deltax.getMesh(), 0.0};
    pos.set("x", mesh->GlobalX(dx),
            "z", TWOPI * static_cast<BoutReal>(dz) / static_cast<BoutReal>(mesh->LocalNz));
    // Generate the analytic solution at the displacements
    a_solution[index] = a_gen->generate(pos);
    b_solution[index] = b_gen->generate(pos);
    c_solution[index] = c_gen->generate(pos);
  }

  // Create the interpolation object from the input options
  auto interp = XZInterpolationFactory::getInstance().create();

  // Interpolate the analytic functions at the displacements
  a_interp = interp->interpolate(a, deltax, deltaz);
  b_interp = interp->interpolate(b, deltax, deltaz);
  c_interp = interp->interpolate(c, deltax, deltaz);

  SAVE_ONCE3(a, a_interp, a_solution);
  SAVE_ONCE3(b, b_interp, b_solution);
  SAVE_ONCE3(c, c_interp, c_solution);

  dump.write();

  BoutFinalise();

  return 0;
}
