
#include <globals.hxx>
#include <dataformat.hxx>
#include <utils.hxx>

bool DataFormat::openr(const string &name, int mype) {
  // Split into base name and extension
  size_t pos = name.find_last_of('.');
  string base(name.substr(0, pos));
  string ext(name.substr(pos+1));
  
  // Insert the processor number between base and extension
  return openr(base + "." + toString(mype) + "." + ext);
}

bool DataFormat::openw(const string &name, int mype, bool append) {
  // Split into base name and extension
  size_t pos = name.find_last_of('.');
  string base(name.substr(0, pos));
  string ext(name.substr(pos+1));
  
  // Insert the processor number between base and extension
  return openw(base + "." + toString(mype) + "." + ext, append);
}

bool DataFormat::setLocalOrigin(int x, int y, int z, int UNUSED(offset_x),
                                int UNUSED(offset_y), int UNUSED(offset_z)) {
  return setGlobalOrigin(x + mesh->OffsetX, y + mesh->OffsetY, z + mesh->OffsetZ);
}
