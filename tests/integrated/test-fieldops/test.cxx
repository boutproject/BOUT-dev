#include <boutmain.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>
#include <bout/mesh.hxx>
#include <assert.h>

#include <interpolation.hxx>

Field3D n;

template <typename T>
void compare(T a, T b, std::string what){
  for (auto i: a){
    if (abs(a[i]-b[i])>1e-13){
      throw BoutException("Error detected in %s at element (%d,%d,%d)\n\tExpected %.30g but got %.30g",what.c_str(),i.x,i.y,i.z,a[i],b[i]);
    }
  }
}


#include "fieldops.cxx"

int physics_init(bool restart) {
  SOLVE_FOR(n); 

#include "test_fieldops.cxx"
  
  return 0;
}


int physics_run(BoutReal time) {
  ddt(n)=0;
  return 0;
}
