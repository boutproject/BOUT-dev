#include "impls/cyclic/laplacexz-cyclic.hxx"
#include "impls/petsc/laplacexz-petsc.hxx"

#include <bout/boutexception.hxx>
#include <bout/invert/laplacexz.hxx>

#include <strings.h>

LaplaceXZ* LaplaceXZ::create(Mesh *m, Options *options) {
  if(options == NULL)
    options = Options::getRoot()->getSection("laplacexz");

  string type;
  options->get("type", type, "cyclic");

  if(strcasecmp(type.c_str(), "cyclic") == 0) {
    return new LaplaceXZcyclic(m, options);
  }else if(strcasecmp(type.c_str(), "petsc") == 0) {
    return new LaplaceXZpetsc(m, options);
  }else {
    throw BoutException("Unknown LaplaceXZ solver type '%s'", type.c_str());
  }
  return 0;
}
