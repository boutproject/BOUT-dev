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

  // We need a vector of Fields because:
  //   1) we don't know at compile time how many we need
  //   2) using a local variable inside the loop causes problems with
  //      dump when the variable goes out of scope
  // We also need to reserve the size to avoid allocations
  // invalidating the pointers the output file has stored. Sections is
  // too large as it includes sections we don't want, but that's ok
  std::vector<Field3D> fields(sections.size());

  for (const auto& section : sections) {
    if (!section.second->isSet("function")) {
      continue;
    }
    fields.emplace_back();
    auto& field = fields.back();
    create_and_dump(field, section.first.c_str());
    bout::globals::dump.write();
  }

  BoutFinalise();

  return 0;
}
