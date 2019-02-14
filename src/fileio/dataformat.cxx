
#include <bout/mesh.hxx>
#include <globals.hxx>
#include <dataformat.hxx>
#include <utils.hxx>

DataFormat::DataFormat(Mesh* mesh_in)
  : mesh(mesh_in==nullptr ? bout::globals::mesh : mesh_in) {}

bool DataFormat::openr(const std::string &name, int mype) {
  // Split into base name and extension
  size_t pos = name.find_last_of('.');
  std::string base(name.substr(0, pos));
  std::string ext(name.substr(pos+1));
  
  // Insert the processor number between base and extension
  return openr(base + "." + toString(mype) + "." + ext);
}

bool DataFormat::openw(const std::string &name, int mype, bool append) {
  // Split into base name and extension
  size_t pos = name.find_last_of('.');
  std::string base(name.substr(0, pos));
  std::string ext(name.substr(pos+1));
  
  // Insert the processor number between base and extension
  return openw(base + "." + toString(mype) + "." + ext, append);
}

bool DataFormat::setLocalOrigin(int x, int y, int z, int UNUSED(offset_x),
                                int UNUSED(offset_y), int UNUSED(offset_z)) {
  return setGlobalOrigin(x + mesh->OffsetX, y + mesh->OffsetY, z + mesh->OffsetZ);
}
