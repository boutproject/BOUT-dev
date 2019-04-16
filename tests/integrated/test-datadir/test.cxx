#include <bout.hxx>

int main(int argc, char **argv) {
  auto ret = BoutInitialise(argc, argv);
  BoutFinalise(false);
  return ret;
}
