/*
 * Check for backtrace after throw
 */

#include <bout/physicsmodel.hxx>
#include <bout.hxx>

void f1() {
  BoutReal a = 1;
  BoutReal b = 0;
  BoutReal c = a / b;
  output.write("c is {:f}\n", c);
  throw BoutException("Tomatoes are red?\n");
}
void f2(int UNUSED(a)) { f1(); }

int f3() {
  f2(1);
  return 0;
}

class Backtrace : public PhysicsModel {
protected:
  int init(bool UNUSED(restarting)) override {
    f3();

    return 1;
  }
  int rhs(BoutReal UNUSED(time)) override { return 1; }
};

BOUTMAIN(Backtrace)
