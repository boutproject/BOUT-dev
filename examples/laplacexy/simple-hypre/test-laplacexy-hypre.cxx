#include <bout.hxx>
#include <bout/invert/laplacexy2_hypre.hxx>
#include <field_factory.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  { 
     /// Create a LaplaceXY object
     LaplaceXY2Hypre laplacexy(mesh);
     
     /// Generate rhs function
     Field2D rhs = FieldFactory::get()->create2D("laplacexy:rhs", Options::getRoot(), mesh);
     
     /// Solution
     Field2D x = 0.0;
     
     x = laplacexy.solve(rhs, x);
     
     SAVE_ONCE2(rhs, x);
     dump.write();  // Save output file
  }
  BoutFinalise();
  return 0;
}

