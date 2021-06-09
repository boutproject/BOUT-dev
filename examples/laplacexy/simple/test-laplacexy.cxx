
#include <bout.hxx>
#include <bout/invert/laplacexy.hxx>
#include <field_factory.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  /// Create a LaplaceXY object
  LaplaceXY laplacexy(bout::globals::mesh);

  /// Generate rhs function
  Field2D rhs = FieldFactory::get()->create2D("laplacexy:rhs", Options::getRoot(),
                                              bout::globals::mesh);

  /// Solution
  Field2D x = 0.0;
  
  x = laplacexy.solve(rhs, x);

  Options dump;
  dump["rhs"] = rhs;
  dump["x"] = x;
  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
  return 0;
}

