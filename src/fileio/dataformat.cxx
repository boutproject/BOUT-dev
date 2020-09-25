
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

void DataFormat::writeFieldAttributes(const std::string& name, const Field& f, bool shiftOutput) {
  // If shiftOutput is true, the data will be written in field-aligned form
  auto direction_y = shiftOutput ? YDirectionType::Aligned : f.getDirectionY();

  setAttribute(name, "cell_location", toString(f.getLocation()));
  setAttribute(name, "direction_y", toString(direction_y));
  setAttribute(name, "direction_z", toString(f.getDirectionZ()));
}

void DataFormat::writeFieldAttributes(const std::string& name, const FieldPerp& f, bool shiftOutput) {
  writeFieldAttributes(name, static_cast<const Field&>(f), shiftOutput);

  auto& fieldmesh = *f.getMesh();
  const int yindex = f.getIndex();
  const int start = fieldmesh.hasBndryLowerY() ? 0 : fieldmesh.ystart;
  const int end = fieldmesh.hasBndryUpperY() ? fieldmesh.LocalNy : fieldmesh.yend + 1;

  // Only use the global y index if it's either an interior (grid)
  // point, or a boundary point. Otherwise, use -1 to indicate a guard
  // cell or an invalid value. The actual FieldPerp value is still
  // written to file
  const int global_yindex =
      (yindex >= start and yindex < end) ? fieldmesh.getGlobalYIndex(yindex) : -1;

  setAttribute(name, "yindex_global", global_yindex);
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

    auto& fieldmesh = *f.getMesh();
    const int start = fieldmesh.hasBndryLowerY() ? 0 : fieldmesh.ystart;
    const int end = fieldmesh.hasBndryUpperY() ? fieldmesh.LocalNy : fieldmesh.yend + 1;

    // Only use the global y index if it's either an interior (grid)
    // point, or a boundary point. Otherwise, use -1 to indicate a
    // guard cell or an invalid value. This may mean that `f` does not
    // get allocated
    const int yindex_local = fieldmesh.getLocalYIndex(yindex_global);
    const int yindex = (yindex_local >= start and yindex_local < end) ? yindex_local : -1;
    f.setIndex(yindex);
  } else {
    // "yindex_global" wasn't present, so this might be an older
    // file. We use the no-boundary form here, such that we get a
    // default value on a grid cell
    f.setIndex(f.getMesh()->getLocalYIndexNoBoundaries(0));
  }
}
