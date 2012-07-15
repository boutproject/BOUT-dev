
#include "meshfactory.hxx"

#include <options.hxx>
#include <strings.h>
#include "impls/bout/boutmesh.hxx"
#include "impls/quilt/quiltmesh.hxx"

MeshFactory *MeshFactory::instance = NULL;

#define MESH_BOUT  "bout"
#define MESH_QUILT "quilt"

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
  
  // Get the type of mesh
  string type;
  options->get("type", type, MESH_BOUT);
  
  if(!strcasecmp(type.c_str(), MESH_QUILT)) {
    return new QuiltMesh();
  }
  
  return new BoutMesh();
}
