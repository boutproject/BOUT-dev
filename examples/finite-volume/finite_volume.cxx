
#include <bout.hxx>

#include <bout/fv_ops.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  Field3D f = 0.0;
  
  mesh->communicate(f);
  
  Field3D g = FV::Grad_par(f);
  
  BoutFinalise();
  return 0;
};
