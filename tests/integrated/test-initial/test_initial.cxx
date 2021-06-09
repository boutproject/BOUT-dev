/*
 * Initial profiles regression test
 *
 * Check that initial profiles are valid, and do
 * not depend on number of processors
 *
 */

#include "initialprofiles.hxx"
#include "bout/physicsmodel.hxx"

#include <algorithm>
#include <vector>

void create_and_dump(Field3D& field, const char* name) {
  initial_profile(name, field);
  bout::globals::dump.add(field, name, false);
}

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  const auto& sections = Options::root().subsections();

  Options dump;

  for (const auto& section : sections) {
    if (!section.second->isSet("function")) {
      continue;
    }
    Field3D field;
    initial_profile(section.first, field);
    dump[section.first] = field;
  }

  bout::writeDefaultOutputFile(dump);

  BoutFinalise();

  return 0;
}
