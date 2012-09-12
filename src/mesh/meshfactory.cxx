
#include "meshfactory.hxx"

#include <options.hxx>
#include <strings.h>
#include <output.hxx>
#include <dataformat.hxx>

#include "impls/bout/boutmesh.hxx"
#include "impls/quilt/quiltmesh.hxx"

MeshFactory *MeshFactory::instance = NULL;

#define MESH_BOUT  "bout"
#define MESH_QUILT "quilt"

static char DEFAULT_GRID[] = "data/bout.grd.nc";

MeshFactory* MeshFactory::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new MeshFactory();
  }
  return instance;
}

Mesh* MeshFactory::createMesh(GridDataSource *source, Options *options) {
  if(options == NULL)
    options = Options::getRoot()->getSection("mesh");
  
  if(source == NULL) {
    output.write("\nGetting grid data source from options\n");
    
    string grid_name;
    if(options->isSet("file")) {
      // Specified mesh file
      options->get("file", grid_name, DEFAULT_GRID);
    }else {
      // Get the global option
      Options::getRoot()->get("grid", grid_name, DEFAULT_GRID);
    }
    
    /// Create a grid file, using specified format if given
    string grid_ext;
    options->get("format", grid_ext, "");
    
    // Create a grid file
    source = (GridDataSource*) new GridFile(data_format( (grid_ext.empty()) ? grid_name.c_str() : grid_ext.c_str() ), 
                                            grid_name.c_str());
  }
  
  // Get the type of mesh
  string type;
  options->get("type", type, MESH_BOUT);
  
  if(!strcasecmp(type.c_str(), MESH_QUILT)) {
    return new QuiltMesh(source, options);
  }
  
  return new BoutMesh(source, options);
}
