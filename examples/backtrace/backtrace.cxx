/*
 * Check for backtrace after throw 
 */

#include <bout.hxx>
#include <bout/physicsmodel.hxx>

class Backtrace : public PhysicsModel {
protected:
  int init(bool UNUSED(restarting)) override;
  int rhs(BoutReal UNUSED(time)) override;
};



void f1(){
  BoutReal a=1;
  BoutReal b=0;
  BoutReal c = a/b;
  output.write("c is {:f}\n",c);
  throw BoutException("Tomatoes are red?\n");
}
void f2(int UNUSED(a)) {
  f1();
}

int f3(){
  f2(1);
  return 0;
}

int Backtrace::init(bool UNUSED(restarting)) {
  f3();
  
  return 1;
}

int Backtrace::rhs(BoutReal UNUSED(time)) {
  
  
  return 1;
}


BOUTMAIN(Backtrace)
