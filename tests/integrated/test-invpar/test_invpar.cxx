/*
 * Test parallel inversion routines
 *
 */

#include <bout.hxx>
#include <derivs.hxx>
#include <invert_parderiv.hxx>
#include <field_factory.hxx>
#include <utils.hxx>

using bout::globals::mesh;

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f(mesh);

  // Get options
  Options &options = Options::root();
  std::string acoef, bcoef, ccoef, dcoef, ecoef, func;
  options.get("acoef", acoef, "1.0");
  options.get("bcoef", bcoef, "-1.0");
  options.get("ccoef", ccoef, "0.0");
  options.get("dcoef", dcoef, "0.0");
  options.get("ecoef", ecoef, "0.0");
  options.get("input", func, "sin(2*y)*(1. + 0.2*exp(cos(z)))");
  auto location = CELL_LOCFromString(options["test_location"].withDefault("CELL_CENTRE"));
  BoutReal tol = options["tol"].withDefault(1e-10);

  auto inv = InvertParFactory::getInstance().create(nullptr, location, mesh);

  Field2D A = f.create2D(acoef, nullptr, nullptr, location);
  Field2D B = f.create2D(bcoef, nullptr, nullptr, location);
  Field2D C = f.create2D(ccoef, nullptr, nullptr, location);
  Field2D D = f.create2D(dcoef, nullptr, nullptr, location);
  Field2D E = f.create2D(ecoef, nullptr, nullptr, location);

  inv->setCoefA(A);
  inv->setCoefB(B);
  inv->setCoefC(C);
  inv->setCoefD(D);
  inv->setCoefE(E);

  Field3D input = f.create3D(func, nullptr, nullptr, location);
  Field3D result = inv->solve(input);
  mesh->communicate(result);

  Field3D deriv = A*result + B*Grad2_par2(result) + C*D2DYDZ(result)
	  + D*D2DZ2(result) + E*DDY(result);

  // Check the result
  int passed = 1;
  for (int y = mesh->ystart; y < mesh->yend; y++) {
    for (int z = 0; z < mesh->LocalNz; z++) {
      output.write("result: [{:d},{:d}] : {:e}, {:e}, {:e}\n", y, z,
                   input(mesh->xstart, y, z), result(mesh->xstart, y, z),
                   deriv(mesh->xstart, y, z));
      if (std::abs(input(mesh->xstart, y, z) - deriv(mesh->xstart, y, z)) > tol) {
        passed = 0;
      }
    }
  }

  int allpassed;
  MPI_Allreduce(&passed, &allpassed, 1, MPI_INT, MPI_MIN, BoutComm::get());

  output << "******* Parallel inversion test case: ";
  if (allpassed) {
    output << "PASSED" << endl;
  } else {
    output << "FAILED" << endl;
  }

  // Saving state to file
  SAVE_ONCE(allpassed);

  // Write data to file
  bout::globals::dump.write();
  bout::globals::dump.close();

  MPI_Barrier(BoutComm::get());

  BoutFinalise();
  return 0;
}
