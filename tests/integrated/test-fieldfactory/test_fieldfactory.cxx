/*
 * FieldFactory regression test
 * 
 * Test the FieldFactory class
 *
 */

#include <bout.hxx>
#include <field_factory.hxx>

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f(bout::globals::mesh);

  Field2D a = f.create2D("2.");
  Field2D b = f.create2D("1 - x");
  Field3D c = f.create3D("sin(3*z)");
  Field3D d = f.create3D("gauss(x-0.5,0.2)*gauss(y)*sin(z)");
  SAVE_ONCE4(a, b, c, d);

  // Write data to file
  bout::globals::dump.write();
  bout::globals::dump.close();

  // Need to wait for all processes to finish writing
  MPI_Barrier(BoutComm::get());

  /// Finished, tidy up and free memory
  BoutFinalise();

  return 0;
}
