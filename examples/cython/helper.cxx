#include <field3d.hxx>

void c_set_f3d_all(Field3D * f3d, double * data){
  int j=0;
  for (auto i : f3d){
    (*f3d)[i]=data[j++];
  }
}
