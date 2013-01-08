/*
 * Smoothing operator tests
 * 
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include <field_factory.hxx>
#include <gyro_average.hxx>

int physics_init(bool restarting) {
  FieldFactory f(mesh);
  
  Field3D input3d = f.create3D("gauss(x-0.5,0.2)*gauss(y-pi)*sin(3*y - z)");
  SAVE_ONCE(input3d);
  
  // Gyro-average
  Field3D pade1 = gyroPade1(input3d, 0.5);
  Field3D pade2 = gyroPade2(input3d, 0.5);
  SAVE_ONCE2(pade1, pade2);
  
  // Write data
  dump.write();
  dump.close();
  
  output << "\nFinished running test. Triggering error to quit\n\n";
  
  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data
  
  return 1;
}

int physics_run(BoutReal t) {
  // Doesn't do anything
  return 1;
}
