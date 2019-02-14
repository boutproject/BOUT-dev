
#include "meshfactory.hxx"

#include <options.hxx>
#include <strings.h>
#include <output.hxx>
#include <dataformat.hxx>
#include <boutexception.hxx>

#include "impls/bout/boutmesh.hxx"

MeshFactory *MeshFactory::instance = nullptr;

/// Name of BoutMesh for Options
#define MESH_BOUT  "bout"

MeshFactory* MeshFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new MeshFactory();
  }
  return instance;
}

Mesh* MeshFactory::createMesh(GridDataSource *source, Options *options) {
  if (options == nullptr)
    options = Options::getRoot()->getSection("mesh");

  if (source == nullptr) {
    std::string grid_name;
    if(options->isSet("file")) {
      // Specified mesh file
      options->get("file", grid_name, "");
      output << "\nGetting grid data from file " << grid_name << endl; 

      /// Create a grid file, using specified format if given
      std::string grid_ext;
      options->get("format", grid_ext, "");
      
      /// Create a grid file
      source = static_cast<GridDataSource *>(new GridFile(
          data_format((grid_ext.empty()) ? grid_name.c_str() : grid_ext.c_str()),
          grid_name.c_str()));
    }else if(Options::getRoot()->isSet("grid")){
      // Get the global option
      Options::getRoot()->get("grid", grid_name, "");
      output << "\nGetting grid data from file " << grid_name << endl; 
      std::string grid_ext;
      Options::getRoot()->get("format", grid_ext, "");

      source = static_cast<GridDataSource *>(new GridFile(
          data_format((grid_ext.empty()) ? grid_name.c_str() : grid_ext.c_str()),
          grid_name.c_str()));
    }else {
      output << "\nGetting grid data from options\n";
      source = static_cast<GridDataSource *>(new GridFromOptions(options));
    }
  }

  // Get the type of mesh
  std::string type;
  options->get("type", type, MESH_BOUT);
  
  if(!strcasecmp(type.c_str(), MESH_BOUT)) {
    return new BoutMesh(source, options);
  }

  throw BoutException("Mesh type not implemented: %s",type.c_str());
}
