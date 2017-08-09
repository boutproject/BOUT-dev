#include <bout.hxx>
#include <bout/fieldgroup.hxx>
#include <bout/assert.hxx>

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  Field3D a;
  Field2D b;

  /// Create a group with one Field3D
  FieldGroup g(a);

  // Add a field2D
  g.add(b);

  // Should have two FieldData objects
  int count = 0;
  for(auto &i : g) {
    output << "FieldData\n";
    count++;
  }
  ASSERT0(count == 2);

  // Should have one Field3D
  count = 0;
  for(auto &i : g.field3d()) {
    output << "Field3D\n";
    count++;
  }
  ASSERT0(count == 1);
  
  BoutFinalise();
  return 0;
}
