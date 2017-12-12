#include "bout.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field3D a = 0.0;

  SAVE_ONCE(a);
  dump.setAttribute("a", "meta", "something useful");
  dump.setAttribute("g12", "value", 42);
  dump.write();

  BoutFinalise();
  return 0;
}
