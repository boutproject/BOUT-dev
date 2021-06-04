/*
 * Gyro operator tests
 * 
 */

#include <bout.hxx>
#include <field_factory.hxx>
#include <gyro_average.hxx>

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f(bout::globals::mesh);

  Field3D input3d = f.create3D("gauss(x-0.5,0.2)*gauss(y-pi)*sin(3*y - z)");
  
  // Gyro-average
  Options::root()["pade1"] = gyroPade1(input3d, 0.5);
  Options::root()["pade2"] = gyroPade2(input3d, 0.5);

  Options::root()["input3d"] = input3d;
  bout::writeDefaultOutputFile();

  output.write("\nFinished running test.\n\n");
  
  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data
  
  bout::checkForUnusedOptions();
  BoutFinalise();
  return 0;
}
