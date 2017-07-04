#include <field3d.hxx>
#include <globals.hxx>

void c_set_f3d_all(Field3D * f3d, double * data){
  int j=0;
  f3d->allocate();
  auto k=data[j];
  for (auto i : (*f3d)){
    (*f3d)[i]=data[j++];
  }
  f3d->operator()(0,0,0)=k;
  printf("%d %d %d  ",f3d->getNx(),f3d->getNy(),f3d->getNz());
  printf("%d written\n",j);
}
void c_get_f3d_all(Field3D * f3d, double * data){
  int j=0;
  for (auto i : (*f3d)){
    data[j++]=(*f3d)[i];
  }
  printf("%d read\n",j);
}

Field3D * fadd( Field3D*a,Field3D*b){
  Field3D * r=new Field3D(*a);
  *r += *b;
  return r;
}
Field3D * fmul( Field3D*a,Field3D*b){
  Field3D * r=new Field3D(*a);
  *r *= *b;
  return r;
}

int getNx( Field3D * a){
  return a->getNx();
}
int getNy( Field3D * a){
  //printf("%d\n",a->getNy());
  return a->getNy();
}
int getNz( Field3D * a){
  return a->getNz();
}

Mesh * c_get_global_mesh(){
  return mesh;
}
