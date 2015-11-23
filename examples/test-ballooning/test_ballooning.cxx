#include <bout.hxx>
#include <initialprofiles.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  Field3D var0, var1;
  
  initial_profile("var0", var0);
  initial_profile("var1", var1);

  mesh->communicate(var0, var1);
  
  SAVE_ONCE2(var0, var1);
  dump.write();
  
  BoutFinalise();
  return 0;
}
