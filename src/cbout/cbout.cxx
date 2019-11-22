#include <cbout/bout.h>

#include <field3d.hxx>

extern "C" Field3D* Field3D_new() {
  return new Field3D();
}

extern "C" Field3D* Field3D_zerofrom(Field3D* field) {
  return new Field3D(zeroFrom(*field));
}

extern "C" Field3D* Field3D_emptyfrom(Field3D* field) {
  return new Field3D(emptyFrom(*field));
}

extern "C" void Field3D_delete(Field3D* field) {
  delete field;
}

