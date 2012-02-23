
#include "meshfactory.hxx"

#include <options.hxx>
#include "impls/bout/boutmesh.hxx"

MeshFactory *MeshFactory::instance = NULL;

MeshFactory* MeshFactory::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new MeshFactory();
  }
  return instance;
}

Mesh* MeshFactory::createMesh(Options *options) {
  if(options == NULL)
    options = Options::getRoot()->getSection("mesh");
  
  return new BoutMesh();
}
