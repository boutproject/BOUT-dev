
#include <globals.hxx>
#include <boutexception.hxx>
#include <strings.h>

#include "laplacefactory.hxx"

#include "impls/serial_tri/serial_tri.hxx"
#include "impls/serial_band/serial_band.hxx"
#include "impls/pdd/pdd.hxx"
#include "impls/spt/spt.hxx"
#include "impls/petsc/petsc_laplace.hxx"
#include "impls/mumps/mumps_laplace.hxx"

#define LAPLACE_SPT  "spt"
#define LAPLACE_PDD  "pdd"
#define LAPLACE_TRI  "tri"
#define LAPLACE_BAND "band"
#define LAPLACE_PETSC "petsc"
#define LAPLACE_MUMPS "mumps"

LaplaceFactory* LaplaceFactory::instance = NULL;

LaplaceFactory* LaplaceFactory::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new LaplaceFactory();
  }
  return instance;
}

Laplacian* LaplaceFactory::createLaplacian(Options *options) {
  if(options == NULL)
    options = Options::getRoot()->getSection("laplace");
  
  string type;
  
  if(mesh->firstX() && mesh->lastX()) {
    // Can use serial algorithm
  
    options->get("type", type, LAPLACE_TRI);
    
    if(strcasecmp(type.c_str(), LAPLACE_TRI) == 0) {
      return new LaplaceSerialTri(options);
    }else if(strcasecmp(type.c_str(), LAPLACE_BAND) == 0) {
      return new LaplaceSerialBand(options);
    }else if(strcasecmp(type.c_str(), LAPLACE_SPT) == 0) {
      return new LaplaceSPT(options);
    }else if(strcasecmp(type.c_str(), LAPLACE_PETSC) == 0) {
      return new LaplacePetsc(options);
    }else if(strcasecmp(type.c_str(), LAPLACE_MUMPS) == 0) {
      return new LaplaceMumps(options);
    }else {
      throw BoutException("Unknown serial Laplacian solver type '%s'", type.c_str());
    }
  }
  
  options->get("type", type, LAPLACE_SPT);

  if(strcasecmp(type.c_str(), LAPLACE_PDD) == 0) {
    return new LaplacePDD(options);
  }else if(strcasecmp(type.c_str(), LAPLACE_SPT) == 0) {
    return new LaplaceSPT(options);
  }else if(strcasecmp(type.c_str(), LAPLACE_PETSC) == 0) {
    return new LaplacePetsc(options);
  }else if(strcasecmp(type.c_str(), LAPLACE_MUMPS) == 0) {
    return new LaplaceMumps(options);
  }else {
    throw BoutException("Unknown parallel Laplacian solver type '%s'", type.c_str());
  }
}

