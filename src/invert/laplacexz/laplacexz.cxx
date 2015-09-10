#include "impls/cyclic/cyclic.hxx"

#include <boutexception.hxx>
#include <bout/invert/laplacexz.hxx>

#include <strings.h>

LaplaceXZ* LaplaceXZ::create(Mesh *m, Options *options) {
  if(options == NULL)
    options = Options::getRoot()->getSection("laplacexz");

  string type;
  options->get("type", type, "cyclic");

  if(strcasecmp(type.c_str(), "cyclic") == 0) {
    return new LaplaceXYcyclic(m, options);
  }else {
    throw BoutException("Unknown LaplaceXZ solver type '%s'", type.c_str());
  }
  return 0;
}
