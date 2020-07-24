/*
 * Smoothing operator tests
 * 
 */

#include <bout.hxx>
#include <smoothing.hxx>
#include <field_factory.hxx>

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f(bout::globals::mesh);

  Field2D input2d = f.create2D("1 + sin(2*y)");
  Field3D input3d = f.create3D("gauss(x-0.5,0.2)*gauss(y-pi)*sin(3*y - z)");
  
  input3d.calcParallelSlices();

  SAVE_ONCE2(input2d, input3d);
  
  // Average in 3D
  Field2D yavg2d = averageY(input2d);
  Field3D yavg3d = averageY(input3d);  
  SAVE_ONCE2(yavg2d, yavg3d);
  
  Field3D sm3d = smooth_y(input3d);
  SAVE_ONCE(sm3d);
  
  // Output data
  bout::globals::dump.write();

  BoutFinalise();
  return 0;
}
