/*
 * Laplacian inversion
 *
 */

#include <bout.hxx>
#include <invert_laplace.hxx>
#include <field_factory.hxx>

int main(int argc, char** argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f{bout::globals::mesh};

  Field3D input = f.create3D("(1-gauss(x-0.5,0.2))*gauss(y-pi)*gauss(z-pi)");
  Field2D a = f.create2D("gauss(x) * sin(y)");
  Field2D c = f.create2D("sin(x) * gauss(x-0.5) * gauss(y-pi)");
  Field2D d = f.create2D("y - pi/2");
  Options::root()["input"] = input;
  Options::root()["a"] = a;
  Options::root()["c"] = c;
  Options::root()["d"] = d;

  auto lap = std::unique_ptr<Laplacian>{Laplacian::create()};

  lap->setCoefA(0.0);
  lap->setCoefC(1.0);
  lap->setCoefD(1.0);
  Options::root()["flag0"] = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Options::root()["flag3"] = lap->solve(input);

  lap->setCoefA(a);
  lap->setInnerBoundaryFlags(0);
  Options::root()["flag0a"] = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Options::root()["flag3a"] = lap->solve(input);

  lap->setCoefC(c);
  lap->setInnerBoundaryFlags(0);
  Options::root()["flag0ac"] = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Options::root()["flag3ac"] = lap->solve(input);

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(0);
  Options::root()["flag0ad"] = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Options::root()["flag3ad"] = lap->solve(input);

  /// Test new interface and INVERT_IN/OUT_SET flags

  Field2D set_to = f.create2D("cos(2*y)*(x - 0.5)");
  Options::root()["set_to"] = set_to;

  lap->setCoefA(0.0);
  lap->setCoefC(1.0);
  lap->setCoefD(1.0);

  lap->setInnerBoundaryFlags(INVERT_SET);
  Options::root()["flagis"] = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Options::root()["flagos"] = lap->solve(input, set_to);

  lap->setCoefA(a);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  lap->setOuterBoundaryFlags(0);
  Options::root()["flagisa"] = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Options::root()["flagosa"] = lap->solve(input, set_to);

  lap->setCoefC(c);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  Options::root()["flagisac"] = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Options::root()["flagosac"] = lap->solve(input, set_to);

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  Options::root()["flagisad"] = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Options::root()["flagosad"] = lap->solve(input, set_to);

  bout::writeDefaultOutputFile();

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  bout::checkForUnusedOptions();
  BoutFinalise();
  return 0;
}
