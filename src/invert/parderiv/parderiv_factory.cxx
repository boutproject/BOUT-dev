
#include <boutcomm.hxx>
#include <invert_parderiv.hxx>
#include <boutexception.hxx>

#include <strings.h>

#include <invert_parderiv.hxx>

#include "impls/cyclic/cyclic.hxx"

ParDerivFactory *ParDerivFactory::instance = nullptr;

/// Default Options section to look for configuration
static const char* default_section = "parderiv";

ParDerivFactory* ParDerivFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new ParDerivFactory();
  }
  return instance;
}

InvertPar* ParDerivFactory::createInvertPar(CELL_LOC location, Mesh *mesh_in) {
  // Get the default options section
  Options *opt = Options::getRoot()->getSection(default_section);

  return createInvertPar(opt, location, mesh_in);
}

InvertPar* ParDerivFactory::createInvertPar(const char* type, Options *opt,
                                            CELL_LOC location, Mesh *mesh_in) {
  int NPES;
  MPI_Comm_size(BoutComm::get(), &NPES);

  if (opt == nullptr)
    opt = Options::getRoot()->getSection(default_section);

  if (!strcasecmp(type, PARDERIVCYCLIC)) {
    return new InvertParCR(opt, location, mesh_in);
  }
  
  throw BoutException("No such ParDeriv solver exists in this build, type: %s", type);
}

InvertPar* ParDerivFactory::createInvertPar(Options *opts, CELL_LOC location,
                                            Mesh *mesh_in) {
  std::string type;
  opts->get("type", type, "cyclic");
  
  return createInvertPar(type.c_str(), opts, location, mesh_in);
}
