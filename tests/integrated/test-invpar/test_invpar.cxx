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

int test(const std::string& acoef, const std::string& bcoef, const std::string& ccoef,
         const std::string& dcoef, const std::string& ecoef, const std::string& func,
         BoutReal tol, FieldFactory& f, CELL_LOC location) {
  output_error.write(
      "acoef={:s} bcoef={:s} ccoef={:s} dcoef={:s} ecoef{:s} with {:s} at {:s} ...",
      acoef, bcoef, ccoef, dcoef, ecoef, func, toString(location));
  auto inv = InvertPar::create(nullptr, location, mesh);

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
  bool success{true};
  int local_ystart = mesh->ystart;
  if (location == CELL_YLOW) {
    // Point at mesh->ystart in 'result' is set by the Neumann boundary condition, so may
    // not agree with 'deriv'
    local_ystart = mesh->ystart + 1;
  }
  for (int y = local_ystart; y < mesh->yend; y++) {
    for (int z = 0; z < mesh->LocalNz; z++) {
      output.write("result: [{:d},{:d}] : {:e}, {:e}, {:e}\n", y, z,
                   input(mesh->xstart, y, z), result(mesh->xstart, y, z),
                   deriv(mesh->xstart, y, z));
      if (std::abs(input(mesh->xstart, y, z) - deriv(mesh->xstart, y, z)) > tol) {
        success = false;
      }
    }
  }
  if (success) {
    output_error.write(" success\n");
  } else {
    output_error.write(" failure\n");
  }
  return success;
}

int main(int argc, char** argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f(mesh);

  // Get options
  Options& options = Options::root();

  BoutReal tol = options["tol"].withDefault(1e-10);

  bool success{true};
  int i = 0;
  while (true) {
    if (not options.isSet(fmt::format("{}_{}", "acoef", i))) {
      break;
    }
    const std::string acoef = options[fmt::format("{}_{}", "acoef", i)].as<std::string>();
    const std::string bcoef =
        options[fmt::format("{}_{}", "bcoef", i)].withDefault("-1.0");
    const std::string ccoef =
        options[fmt::format("{}_{}", "ccoef", i)].withDefault("0.0");
    const std::string dcoef =
        options[fmt::format("{}_{}", "dcoef", i)].withDefault("0.0");
    const std::string ecoef =
        options[fmt::format("{}_{}", "ecoef", i)].withDefault("0.0");
    const std::string func = options[fmt::format("{}_{}", "input_field", i)].withDefault(
        "sin(2*y)*(1. + 0.2*exp(cos(z)))");

    for (auto location : {CELL_CENTRE, CELL_YLOW, CELL_XLOW, CELL_ZLOW}) {
      success &= test(acoef, bcoef, ccoef, dcoef, ecoef, func, tol, f, location);
    }
    i += 1;
  }
  bool allsuccess{true};
  MPI_Allreduce(&success, &allsuccess, 1, MPI_C_BOOL, MPI_LAND, BoutComm::get());

  output << "******* Parallel inversion test case: ";
  if (allsuccess) {
    output << "PASSED" << endl;
  } else {
    output << "FAILED" << endl;
  }

  MPI_Barrier(BoutComm::get());

  BoutFinalise();
  return !allsuccess;
}
