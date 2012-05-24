
#include <globals.hxx>
#include <dataformat.hxx>


bool DataFormat::setLocalOrigin(int x, int y, int z) {
  return setGlobalOrigin(x + mesh->OffsetX, y + mesh->OffsetY, z + mesh->OffsetZ);
}

