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

  FieldFactory f(mesh);

  Field3D input = f.create3D("(1-gauss(x-0.5,0.2))*gauss(y-pi)*gauss(z-pi)");
  Field2D a = f.create2D("gauss(x,0.24) * sin(y)");
  Field2D c = f.create2D("sin(x) * gauss(x-0.5) * gauss(y-pi)");
  Field2D d = f.create2D("y - pi/2");
  SAVE_ONCE4(input, a, c, d);

  Field3D flag0 = invert_laplace(input, 0);
  Field3D flag3 = invert_laplace(input, 3);
  SAVE_ONCE2(flag0, flag3);

  Field3D flag0a = invert_laplace(input, 0, &a);
  Field3D flag3a = invert_laplace(input, 3, &a);
  SAVE_ONCE2(flag0a, flag3a);

  Field3D res0a  = Delp2(flag0a);
  SAVE_ONCE(res0a);

  Field3D flag0ac = invert_laplace(input, 0, &a, &c);
  Field3D flag3ac = invert_laplace(input, 3, &a, &c);
  SAVE_ONCE2(flag0ac, flag3ac);

  Field3D flag0ad = invert_laplace(input, 0, &a, nullptr, &d);
  Field3D flag3ad = invert_laplace(input, 3, &a, nullptr, &d);
  SAVE_ONCE2(flag0ad, flag3ad);

  /// Test new interface and INVERT_IN/OUT_SET flags

  Field2D set_to = f.create2D("cos(2*y)*(x - 0.5)");
  SAVE_ONCE(set_to);

  Laplacian *lap = Laplacian::create();
  lap->setFlags(4096);
  Field3D flagis = lap->solve(input, set_to);
  lap->setFlags(8192);
  Field3D flagos = lap->solve(input, set_to);
  SAVE_ONCE2(flagis, flagos);

  lap->setCoefA(a);
  lap->setFlags(4096);
  Field3D flagisa = lap->solve(input, set_to);
  lap->setFlags(8192);
  Field3D flagosa = lap->solve(input, set_to);
  SAVE_ONCE2(flagisa, flagosa);

  lap->setCoefC(c);
  lap->setFlags(4096);
  Field3D flagisac = lap->solve(input, set_to);
  lap->setFlags(8192);
  Field3D flagosac = lap->solve(input, set_to);
  SAVE_ONCE2(flagisac, flagosac);

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setFlags(4096);
  Field3D flagisad = lap->solve(input, set_to);
  lap->setFlags(8192);
  Field3D flagosad = lap->solve(input, set_to);
  SAVE_ONCE2(flagisad, flagosad);

  // Delete Laplacian when done
  delete lap;

  // Write and close the output file

  dump.write();
  dump.close();

  output << "\nFinished running test. Triggering error to quit\n\n";

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  BoutFinalise();
  return 0;
}
