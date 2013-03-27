
#include "meshfactory.hxx"

#include <options.hxx>
#include <strings.h>
#include <output.hxx>
#include <dataformat.hxx>
#include <boutexception.hxx>

#include "impls/bout/boutmesh.hxx"
//#include "impls/quilt/quiltmesh.hxx"

MeshFactory *MeshFactory::instance = NULL;

#define MESH_BOUT  "bout"
//#define MESH_QUILT "quilt"

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
    string grid_name;
    if(options->isSet("file")) {
      // Specified mesh file
      options->get("file", grid_name, "");
      output << "\nGetting grid data from file " << grid_name << endl; 

      /// Create a grid file, using specified format if given
      string grid_ext;
      options->get("format", grid_ext, "");
      
      /// Create a grid file
      source = (GridDataSource*) new GridFile(data_format( (grid_ext.empty()) ? grid_name.c_str() : grid_ext.c_str() ), 
                                              grid_name.c_str());
    }else if(Options::getRoot()->isSet("grid")){
      // Get the global option
      Options::getRoot()->get("grid", grid_name, "");
      output << "\nGetting grid data from file " << grid_name << endl; 
      string grid_ext;
      Options::getRoot()->get("format", grid_ext, "");
      
      source = (GridDataSource*) new GridFile(data_format( (grid_ext.empty()) ? grid_name.c_str() : grid_ext.c_str() ), 
                                              grid_name.c_str());
    }else {
      output << "\nGetting grid data from options\n"; 
      source = (GridDataSource*) new GridFromOptions(options);
    }
  }
  
  if(!source->isValid())
    throw BoutException("Invalid data source. Maybe the wrong grid file name?");

  // Get the type of mesh
  string type;
  options->get("type", type, MESH_BOUT);
  
/*
  if(!strcasecmp(type.c_str(), MESH_QUILT)) {
    return new QuiltMesh(source, options);
  }
*/
  
  return new BoutMesh(source, options);
}
