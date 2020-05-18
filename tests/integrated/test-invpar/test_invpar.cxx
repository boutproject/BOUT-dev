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

  auto inv = InvertPar::Create();
  FieldFactory f(mesh);

  // Get options
  Options *options = Options::getRoot();
  std::string acoef, bcoef, ccoef, dcoef, ecoef, func;
  options->get("acoef", acoef, "1.0");
  options->get("bcoef", bcoef, "-1.0");
  options->get("ccoef", ccoef, "0.0");
  options->get("dcoef", dcoef, "0.0");
  options->get("ecoef", ecoef, "0.0");
  options->get("input", func, "sin(2*y)");
  BoutReal tol;
  OPTION(options, tol, 1e-10);

  Field2D A = f.create2D(acoef);
  Field2D B = f.create2D(bcoef);
  Field2D C = f.create2D(ccoef);
  Field2D D = f.create2D(dcoef);
  Field2D E = f.create2D(ecoef);

  inv->setCoefA(A);
  inv->setCoefB(B);
  inv->setCoefC(C);
  inv->setCoefD(D);
  inv->setCoefE(E);

  Field3D input = f.create3D(func);
  Field3D result = inv->solve(input);
  mesh->communicate(result);

  Field3D deriv = A*result + B*Grad2_par2(result) + C*D2DYDZ(result)
	  + D*D2DZ2(result) + E*DDY(result);

  // Check the result
  int passed = 1;
  for (int y = 2; y < mesh->LocalNy - 2; y++) {
    for (int z = 0; z < mesh->LocalNz; z++) {
      output.write("result: [{:d},{:d}] : {:e}, {:e}, {:e}\n", y, z, input(2, y, z),
                   result(2, y, z), deriv(2, y, z));
      if (std::abs(input(2, y, z) - deriv(2, y, z)) > tol)
        passed = 0;
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
