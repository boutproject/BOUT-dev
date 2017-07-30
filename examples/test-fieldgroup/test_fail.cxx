/// This is a test of FieldGroup which should fail to compile
///
///

#include <bout/fieldgroup.hxx>

int main(int argc, char **argv) {

  // Construct a FieldGroup with an integer
  // Should fail to compile
  // (hopefully with a useful error message)
  
  int i;
  FieldGroup g(i); 

  return 0;
}
