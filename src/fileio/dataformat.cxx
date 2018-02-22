
#include <bout/globals.hxx>
#include <bout/dataformat.hxx>
#include <bout/utils.hxx>

bool DataFormat::openr(const string &name, int mype) {
  // Split into base name and extension
  size_t pos = name.find_last_of(".");
  string base(name.substr(0, pos));
  string ext(name.substr(pos+1));
  
  // Insert the processor number between base and extension
  return openr(base + "." + toString(mype) + "." + ext);
}

bool DataFormat::openw(const string &name, int mype, bool append) {
  // Split into base name and extension
  size_t pos = name.find_last_of(".");
  string base(name.substr(0, pos));
  string ext(name.substr(pos+1));
  
  // Insert the processor number between base and extension
  return openw(base + "." + toString(mype) + "." + ext, append);
}

bool DataFormat::setLocalOrigin(int x, int y, int z, int offset_x, int offset_y, int offset_z) {
  return setGlobalOrigin(x + mesh->OffsetX, y + mesh->OffsetY, z + mesh->OffsetZ);
}

