
#include <bout.hxx>
#include <bout/invert/laplacexy.hxx>
#include <field_factory.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  { // braces make sure laplacexy's destructor is called before BoutFinalise(). Otherwise PETSc 3.7 raises an error.
    /// Create a LaplaceXY object
    LaplaceXY laplacexy(mesh);
    
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

