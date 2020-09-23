
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
  // This function should not be called from the DataFormat in GridFromFile, which is
  // created before the Mesh - then DataFormat::mesh would be nullptr.
  ASSERT1(mesh != nullptr);
  return setGlobalOrigin(x + mesh->OffsetX, y + mesh->OffsetY, z + mesh->OffsetZ);
}

void DataFormat::writeFieldAttributes(const std::string& name, const Field& f) {
  setAttribute(name, "cell_location", toString(f.getLocation()));
  setAttribute(name, "direction_y", toString(f.getDirectionY()));
  setAttribute(name, "direction_z", toString(f.getDirectionZ()));
}

void DataFormat::writeFieldAttributes(const std::string& name, const FieldPerp& f) {
  writeFieldAttributes(name, static_cast<const Field&>(f));

  auto& fieldmesh = *f.getMesh();
  int yindex = f.getIndex();
  if (yindex >= fieldmesh.ystart and yindex <= fieldmesh.yend) {
    // write global y-index as attribute
    setAttribute(name, "yindex_global", fieldmesh.getGlobalYIndex(yindex));
  } else {
    // y-index is not valid, set global y-index to -1 to indicate 'not-valid'
    setAttribute(name, "yindex_global", -1);
  }
}

void DataFormat::readFieldAttributes(const std::string& name, Field& f) {
  std::string location_string;
  if (getAttribute(name, "cell_location", location_string)) {
    f.setLocation(CELL_LOCFromString(location_string));
  }

  std::string direction_y_string;
  if (getAttribute(name, "direction_y", direction_y_string)) {
    f.setDirectionY(YDirectionTypeFromString(direction_y_string));
  }

  std::string direction_z_string;
  if (getAttribute(name, "direction_z", direction_z_string)) {
    f.setDirectionZ(ZDirectionTypeFromString(direction_z_string));
  }
}

void DataFormat::readFieldAttributes(const std::string& name, FieldPerp& f) {
  readFieldAttributes(name, static_cast<Field&>(f));

  int yindex_global = 0;
  // Note: don't use DataFormat::mesh variable, because it may be null if the DataFormat
  // is part of a GridFromFile, which is created before the Mesh.
  if (getAttribute(name, "yindex_global", yindex_global)) {
    f.setIndex(f.getMesh()->getLocalYIndex(yindex_global));
  } else {
    // No boundary form here, so default value is on a grid cell
    f.setIndex(f.getMesh()->getLocalYIndexNoBoundaries(0));
  }
}
