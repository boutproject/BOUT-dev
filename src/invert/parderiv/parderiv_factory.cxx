
#include <boutcomm.hxx>
#include <invert_parderiv.hxx>
#include <boutexception.hxx>

#include <strings.h>

#include "parderiv_factory.hxx"

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

InvertPar* ParDerivFactory::createInvertPar(Mesh *mesh_in) {
  // Get the default options section
  Options *opt = Options::getRoot()->getSection(default_section);

  return createInvertPar(opt, mesh_in);
}

InvertPar* ParDerivFactory::createInvertPar(const char* type, Options *opt, Mesh *mesh_in) {
  int NPES;
  MPI_Comm_size(BoutComm::get(), &NPES);

  if (opt == nullptr)
    opt = Options::getRoot()->getSection(default_section);

  if (!strcasecmp(type, PARDERIVCYCLIC)) {
    return new InvertParCR(opt, mesh_in);
  }
  
  throw BoutException("No such ParDeriv solver exists in this build, type: %s", type);
}

InvertPar* ParDerivFactory::createInvertPar(Options *opts, Mesh *mesh_in) {
  std::string type;
  opts->get("type", type, "cyclic");
  
  return createInvertPar(type.c_str(), opts, mesh_in);
}
