/*
 * Smoothing operator tests
 * 
 */

#include <bout.hxx>
#include <smoothing.hxx>
#include <field_factory.hxx>

int main(int argc, char** argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f(bout::globals::mesh);

  Field2D input2d = f.create2D("1 + sin(2*y)");
  Field3D input3d = f.create3D("gauss(x-0.5,0.2)*gauss(y-pi)*sin(3*y - z)");

  input3d.calcParallelSlices();

  Options::root()["input2d"] = input2d;
  Options::root()["input3d"] = input3d;

  // Average in 3D
  Options::root()["yavg2d"] = averageY(input2d);
  Options::root()["yavg3d"] = averageY(input3d);

  Options::root()["sm3d"] = smooth_y(input3d);

  bout::writeDefaultOutputFile();

  bout::checkForUnusedOptions();

  BoutFinalise();
  return 0;
}
