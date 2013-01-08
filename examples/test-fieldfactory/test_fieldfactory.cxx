/*
 * FieldFactory regression test
 * 
 * Test the FieldFactory class
 *
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include <field_factory.hxx>

int physics_init(bool restarting) {
  FieldFactory f(mesh);
  
  Field2D a = f.create2D("2.");
  Field2D b = f.create2D("1 - x");
  Field3D c = f.create3D("sin(3*z)");
  Field3D d = f.create3D("gauss(x-0.5,0.2)*gauss(y)*sin(z)");
  SAVE_ONCE4(a, b, c, d);

  
  // Write data to file
  dump.write();
  dump.close();
  
  // Need to wait for all processes to finish writing
  MPI_Barrier(BoutComm::get());

  // Send an error code so quits
  return 1;
}

int physics_run(BoutReal t) {
  // Doesn't do anything
  return 1;
}
