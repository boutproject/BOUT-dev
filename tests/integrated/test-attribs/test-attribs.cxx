#include "bout.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field3D a = 0.0;

  SAVE_ONCE(a);
  bout::globals::dump.setAttribute("a", "meta", "something useful");
  bout::globals::dump.setAttribute("g12", "value", 42);
  bout::globals::dump.write();

  BoutFinalise();
  return 0;
}
