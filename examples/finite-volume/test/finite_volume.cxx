
#include <bout.hxx>

#include <bout/fv_ops.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  Field3D f = 0.0;
  
  f.getMesh()->communicate(f);
  
  Field3D g = FV::D4DY4_Index(f);
  
  BoutFinalise();
  return 0;
};
