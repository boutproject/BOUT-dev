
#include <bout/bout.hxx>
#include <bout/field_factory.hxx>
#include <bout/invert/laplacexy.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  /// Create a LaplaceXY object
  auto laplacexy = LaplaceXY::create(bout::globals::mesh);

  /// Generate rhs function
  Field2D rhs = FieldFactory::get()->create2D("laplacexy:rhs", Options::getRoot(),
                                              bout::globals::mesh);

  /// Solution
  Field2D result = laplacexy->solve(rhs, 0.0);

  Options dump;
  dump["rhs"] = rhs;
  dump["result"] = result;
  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
  return 0;
}
