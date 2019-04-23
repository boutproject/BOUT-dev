
#include <options.hxx>
#include <boutexception.hxx>
#include <strings.h>

#include "laplace3d_factory.hxx"

// Include implementations here
#include "impls/petsc/petsc_laplace3d.hxx"

Laplace3DFactory* Laplace3DFactory::instance = nullptr;

Laplace3DFactory* Laplace3DFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new Laplace3DFactory();
  }
  return instance;
}

Laplace3D* Laplace3DFactory::createLaplace3D(Options* options) {
  if (!options)
    options = Options::getRoot()->getSection("laplace3d");

  std::string type = (*options)["type"].withDefault("petsc");

  // Add tests for available solvers here. See src/invert/laplace/laplacefactory.cxx
  if (type == "petsc") {
    return new Laplace3DPetsc(options);
  }

  throw BoutException("Unknown Laplace3D solver type '%s'", type.c_str());
}
