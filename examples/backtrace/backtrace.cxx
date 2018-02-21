/*
 * Check for backtrace after throw 
 */

#include <bout.hxx>
#include <bout/boutmain.hxx>


void f1(){
  BoutReal a=1;
  BoutReal b=0;
  BoutReal c = a/b;
  output.write("c is %f\n",c);
  throw BoutException("Tomatoes are red?\n");
}
void f2(int a){
  f1();
}

int f3(){
  f2(1);
  return 0;
}

int physics_init(bool restarting) {
  f3();
  
  return 1;
}

int physics_run(BoutReal time) {
  
  
  return 1;
}
