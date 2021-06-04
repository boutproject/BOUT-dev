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

int test(std::string acoef, std::string bcoef, std::string ccoef, std::string dcoef,
         std::string ecoef, std::string func, BoutReal tol, FieldFactory& f,
         CELL_LOC location) {
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
  int exit = 0;
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
        exit = 1;
      }
    }
  }
  return exit;
}

std::string get(Options& options, const char* name, int i, const std::string& def) {
  char* fname = new char[64];
  snprintf(fname, 63, "%s_%d", name, i);
  std::string val;
  options.get(fname, val, def);
  return val;
}

int main(int argc, char** argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f(mesh);

  // Get options
  Options& options = Options::root();
  std::string acoef, bcoef, ccoef, dcoef, ecoef, func;

  BoutReal tol = options["tol"].withDefault(1e-10);

  int exit = 0;
  int i = 0;
  while (true) {
    acoef = get(options, "acoef", i, "not_set");
    if (acoef == "not_set") {
      break;
    }
    bcoef = get(options, "bcoef", i, "-1.0");
    ccoef = get(options, "ccoef", i, "0.0");
    dcoef = get(options, "dcoef", i, "0.0");
    ecoef = get(options, "ecoef", i, "0.0");
    func = get(options, "input_field", i, "sin(2*y)*(1. + 0.2*exp(cos(z)))");

    for (auto location : {CELL_CENTRE, CELL_YLOW, CELL_XLOW, CELL_ZLOW}) {
      exit += test(acoef, bcoef, ccoef, dcoef, ecoef, func, tol, f, location);
    }
    i += 1;
  }
  int allexit;
  MPI_Allreduce(&exit, &allexit, 1, MPI_INT, MPI_MAX, BoutComm::get());

  output << "******* Parallel inversion test case: ";
  if (allexit == 0) {
    output << "PASSED" << endl;
  } else {
    output << "FAILED" << endl;
    allexit = 1;
  }

  MPI_Barrier(BoutComm::get());

  BoutFinalise();
  return allexit;
}
