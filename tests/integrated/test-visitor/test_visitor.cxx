
#include <bout.hxx>
#include <bout/field_visitor.hxx>

class MyVisitor : public FieldVisitor {
  void accept(Field2D &UNUSED(f)) override { output << "Field2D" << endl; }

  void accept(Field3D &UNUSED(f)) override { output << "Field3D" << endl; }
  
  void accept(FieldPerp &UNUSED(f)) override { output << "FieldPerp" << endl; }
  
  void accept(Vector2D &UNUSED(f)) override { output << "Vector2D\n"; }
  
  void accept(Vector3D &UNUSED(f)) override { output << "Vector3D" << endl; }
};

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  MyVisitor v;

  Field3D f {0.0};

  f.accept(v);
  
  
  BoutFinalise();
  return 0;
}
