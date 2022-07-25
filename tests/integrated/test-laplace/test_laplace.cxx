/*
 * Laplacian inversion
 *
 */

#include <bout.hxx>
#include <invert_laplace.hxx>
#include <field_factory.hxx>

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f{bout::globals::mesh};

  Field3D input = f.create3D("(1-gauss(x-0.5,0.2))*gauss(y-pi)*gauss(z-pi)");
  Field2D a = f.create2D("gauss(x) * sin(y)");
  Field2D c = f.create2D("sin(x) * gauss(x-0.5) * gauss(y-pi)");
  Field2D d = f.create2D("y - pi/2");
  SAVE_ONCE4(input, a, c, d);

  auto lap = std::unique_ptr<Laplacian>{Laplacian::create()};

  lap->setCoefA(0.0);
  lap->setCoefC(1.0);
  lap->setCoefD(1.0);
  Field3D flag0 = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Field3D flag3 = lap->solve(input);
  SAVE_ONCE2(flag0, flag3);

  lap->setCoefA(a);
  lap->setInnerBoundaryFlags(0);
  Field3D flag0a = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Field3D flag3a = lap->solve(input);
  SAVE_ONCE2(flag0a, flag3a);

  lap->setCoefC(c);
  lap->setInnerBoundaryFlags(0);
  Field3D flag0ac = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Field3D flag3ac = lap->solve(input);
  SAVE_ONCE2(flag0ac, flag3ac);

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(0);
  Field3D flag0ad = lap->solve(input);
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Field3D flag3ad = lap->solve(input);
  SAVE_ONCE2(flag0ad, flag3ad);

  /// Test new interface and INVERT_IN/OUT_SET flags

  Field2D set_to = f.create2D("cos(2*y)*(x - 0.5)");
  SAVE_ONCE(set_to);

  lap->setCoefA(0.0);
  lap->setCoefC(1.0);
  lap->setCoefD(1.0);

  lap->setInnerBoundaryFlags(INVERT_SET);
  Field3D flagis = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Field3D flagos = lap->solve(input, set_to);
  SAVE_ONCE2(flagis, flagos);

  lap->setCoefA(a);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  lap->setOuterBoundaryFlags(0);
  Field3D flagisa = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Field3D flagosa = lap->solve(input, set_to);
  SAVE_ONCE2(flagisa, flagosa);

  lap->setCoefC(c);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  Field3D flagisac = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Field3D flagosac = lap->solve(input, set_to);
  SAVE_ONCE2(flagisac, flagosac);

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  Field3D flagisad = lap->solve(input, set_to);
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Field3D flagosad = lap->solve(input, set_to);
  SAVE_ONCE2(flagisad, flagosad);

  // Write and close the output file
  bout::globals::dump.write();
  bout::globals::dump.close();

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  bout::checkForUnusedOptions();
  BoutFinalise();
  return 0;
}
