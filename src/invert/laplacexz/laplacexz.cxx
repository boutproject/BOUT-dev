#include "impls/cyclic/laplacexz-cyclic.hxx"
#include "impls/petsc/laplacexz-petsc.hxx"

#include <boutexception.hxx>
#include <bout/invert/laplacexz.hxx>

#include <strings.h>

LaplaceXZ* LaplaceXZ::create(Mesh *m, Options *options, const CELL_LOC loc) {
  if (m == nullptr) {
    // use global mesh
    m = bout::globals::mesh;
  }

  if (options == nullptr) {
    options = &(Options::root()["laplacexz"]);
  }

  std::string type = (*options)["type"].withDefault("cyclic");

  if (strcasecmp(type.c_str(), "cyclic") == 0) {
    return new LaplaceXZcyclic(m, options, loc);
  } else if(strcasecmp(type.c_str(), "petsc") == 0) {
    return new LaplaceXZpetsc(m, options, loc);
  } else {
    throw BoutException("Unknown LaplaceXZ solver type '%s'", type.c_str());
  }
  return nullptr;
}
