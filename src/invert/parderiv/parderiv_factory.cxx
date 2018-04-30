
#include <boutcomm.hxx>
#include <invert_parderiv.hxx>
#include <boutexception.hxx>

#include <strings.h>

#include "parderiv_factory.hxx"

#include "impls/serial/serial.hxx"
#include "impls/cyclic/cyclic.hxx"

ParDerivFactory* ParDerivFactory::instance = NULL;

/// Default Options section to look for configuration
static const char* default_section = "parderiv";

ParDerivFactory* ParDerivFactory::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new ParDerivFactory();
  }
  return instance;
}

InvertPar* ParDerivFactory::createInvertPar() {
  // Get the default options section
  Options *opt = Options::getRoot()->getSection(default_section);

  return createInvertPar( opt );
}

InvertPar* ParDerivFactory::createInvertPar(const char* type, Options *opt) {
  int NPES;
  MPI_Comm_size(BoutComm::get(), &NPES);
  
  if(opt == NULL)
    opt = Options::getRoot()->getSection(default_section);

  if(!strcasecmp(type, PARDERIVSERIAL)) {
    return new InvertParSerial(opt);
  }else if(!strcasecmp(type, PARDERIVCYCLIC)) {
    return new InvertParCR(opt);
  }
  
  throw BoutException("No such ParDeriv solver exists in this build, type: %s", type);
}

InvertPar* ParDerivFactory::createInvertPar(Options *opts) {
  string type;
  opts->get("type", type, "cyclic");
  
  return createInvertPar(type.c_str(), opts);
}
