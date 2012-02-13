
#include <globals.hxx>
#include <invert_parderiv.hxx>
#include <boutexception.hxx>

#include <strings.h>

#include "parderiv_factory.hxx"

#include "impls/serial/serial.hxx"
#include "impls/cyclic/cyclic.hxx"

ParDerivFactory* ParDerivFactory::instance = NULL;

ParDerivFactory* ParDerivFactory::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new ParDerivFactory();
  }
  return instance;
}

InvertPar* ParDerivFactory::createInvertPar() {
  int NPES;
  MPI_Comm_size(BoutComm::get(), &NPES);
  
  if(NPES > 1) {
    // Need a parallel solver
    return new InvertParCR();
  }
  // Only one processor, so serial solver will do
  return new InvertParSerial();
}

InvertPar* ParDerivFactory::createInvertPar(const char* type) {
  int NPES;
  MPI_Comm_size(BoutComm::get(), &NPES);
  
  if(!strcasecmp(type, PARDERIVSERIAL)) {
    return new InvertParSerial;
  }else if(!strcasecmp(type, PARDERIVCYCLIC)) {
    return new InvertParCR;
  }
  
  throw BoutException("No such ParDeriv solver exists in this build, type: %s", type);
}

InvertPar* ParDerivFactory::createInvertPar(Options *opts) {
  string type;
  opts->get("type", type, "cyclic");
  
  return createInvertPar(type.c_str());
}
