#include <bout.hxx>
#include <bout/fieldgroup.hxx>

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  Field3D a;
  Field2D b;

  FieldGroup g(a);
  
  g.add(b);

  for(auto &i : g) {
    output << "FieldData\n";
  }

  for(auto &i : g.field3d()) {
    output << "Field3D\n";
  }
  
  BoutFinalise();
  return 0;
}
