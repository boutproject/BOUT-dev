#include "bout.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field3D a = 1.23;

  SAVE_ONCE(a);
  dump.write();

  BoutFinalise();
  return 0;
}
