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
  Options dump;

  Field3D input = f.create3D("(1-gauss(x-0.5,0.2))*gauss(y-pi)*gauss(z-pi)");
  Field2D a = f.create2D("gauss(x) * sin(y)");
  Field2D c = f.create2D("sin(x) * gauss(x-0.5) * gauss(y-pi)");
  Field2D d = f.create2D("y - pi/2");
  dump["input"] = input;
  dump["a"] = a;
  dump["c"] = c;
  dump["d"] = d;

  auto lap = std::unique_ptr<Laplacian>{Laplacian::create()};

  lap->setCoefA(0.0);
  lap->setCoefC(1.0);
  lap->setCoefD(1.0);
  dump["flag0"] = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  dump["flag3"] = lap->solve(input);

  lap->setCoefA(a);
  lap->setInnerBoundaryFlags(0);
  dump["flag0a"] = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  dump["flag3a"] = lap->solve(input);

  lap->setCoefC(c);
  lap->setInnerBoundaryFlags(0);
  dump["flag0ac"] = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  dump["flag3ac"] = lap->solve(input);

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(0);
  dump["flag0ad"] = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  dump["flag3ad"] = lap->solve(input);

  /// Test new interface and INVERT_IN/OUT_SET flags

  Field2D set_to = f.create2D("cos(2*y)*(x - 0.5)");
  dump["set_to"] = set_to;

  lap->setCoefA(0.0);
  lap->setCoefC(1.0);
  lap->setCoefD(1.0);

  lap->setInnerBoundaryFlags(INVERT_SET);
  dump["flagis"] = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  dump["flagos"] = lap->solve(input, set_to);

  lap->setCoefA(a);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  lap->setOuterBoundaryFlags(0);
  dump["flagisa"] = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  dump["flagosa"] = lap->solve(input, set_to);

  lap->setCoefC(c);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  dump["flagisac"] = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  dump["flagosac"] = lap->solve(input, set_to);

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  dump["flagisad"] = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  dump["flagosad"] = lap->solve(input, set_to);

  bout::writeDefaultOutputFile(dump);

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  bout::checkForUnusedOptions();
  BoutFinalise();
  return 0;
}
