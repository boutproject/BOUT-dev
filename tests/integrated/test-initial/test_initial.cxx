/*
 * Initial profiles regression test
 *
 * Check that initial profiles are valid, and do
 * not depend on number of processors
 *
 */

#include "bout/physicsmodel.hxx"
#include "bout/initialprofiles.hxx"

void create_and_dump(Field3D& field, const char* name) {
  initial_profile(name, field);
  dump.add(field, name, false);
}
void create_and_dump(Field2D& field, const char* name) {
  initial_profile(name, field);
  dump.add(field, name, false);
}

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  // Save the coordinate arrays to make the python bit easier
  Field3D var_x;
  create_and_dump(var_x, "var_x");

  Field3D var_y;
  create_and_dump(var_y, "var_y");

  Field3D var_z;
  create_and_dump(var_z, "var_z");

  // Include the functions to be tested
  // ./runtest generates this file by reading the list of variables in
  // data/BOUT.inp, excluding var_{x,y,z}
#include "test_functions.cxx"

  dump.write();

  BoutFinalise();

  return 0;
}
