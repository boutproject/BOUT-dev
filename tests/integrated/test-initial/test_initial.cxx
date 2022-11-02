/*
 * Initial profiles regression test
 *
 * Check that initial profiles are valid, and do
 * not depend on number of processors
 *
 */

#include "initialprofiles.hxx"
#include "bout/physicsmodel.hxx"

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
