
#include <bout/options.hxx>
#include <bout/boutexception.hxx>
#include <strings.h>

#include "laplace3d_factory.hxx"

// Include implementations here
#include "impls/petsc/petsc_laplace3d.hxx"

Laplace3DFactory* Laplace3DFactory::instance = 0;

Laplace3DFactory* Laplace3DFactory::getInstance() {
  if(instance == 0) {
    // Create the singleton object
    instance = new Laplace3DFactory();
  }
  return instance;
}

Laplace3D* Laplace3DFactory::createLaplace3D(Options *options) {
  if(!options)
    options = Options::getRoot()->getSection("laplace3d");
  
  string type;
  options->get("type", type, "petsc");

  // Add tests for available solvers here. See src/invert/laplace/laplacefactory.cxx
  if(strcasecmp(type.c_str(), "petsc") == 0) {
    return new Laplace3DPetsc(options);
  }
  
  throw BoutException("Unknown Laplace3D solver type '%s'", type.c_str());
}
