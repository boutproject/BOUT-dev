#include <bout.hxx>
#include <bout/invert/laplacexy2_hypre.hxx>
#include <field_factory.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  { 
     /// Create a LaplaceXY object
     LaplaceXY2Hypre laplacexy(bout::globals::mesh);
     
     /// Generate rhs function
     Field2D rhs = FieldFactory::get()->create2D("laplacexy:rhs", Options::getRoot(), bout::globals::mesh);
     
     /// Solution
     Field2D x = 0.0;
     
     x = laplacexy.solve(rhs, x);
     
     SAVE_ONCE2(rhs, x);
     bout::globals::dump.write();  // Save output file
  }
  BoutFinalise();
  return 0;
}

