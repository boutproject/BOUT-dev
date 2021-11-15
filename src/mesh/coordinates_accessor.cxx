#include "bout/coordinates_accessor.hxx"

#include "bout/mesh.hxx"

#include <map>

namespace {
  /// Associate each Coordinates object with an Array object
  /// which contains the coordinates data in striped form.
  ///
  /// Note: This association could perhaps be done by putting
  ///       the Array inside Coordinates, but this keeps things decoupled
  std::map<const Coordinates*, Array<BoutReal>> coords_store;
}

CoordinatesAccessor::CoordinatesAccessor(const Coordinates* coords) {
  ASSERT0(coords != nullptr);

  // Size of the mesh in Z. Used to convert 3D -> 2D index
  Mesh* mesh = coords->dx.getMesh();
  mesh_nz = mesh->LocalNz;

  auto search = coords_store.find(coords);
  if (search != coords_store.end()) {
    // Found, so get the pointer to the data
    data = search->second.begin();
    return;
  }

  // Not yet created, so create the array and copy the data

  // Work out the size needed
  int array_size = stripe_size * mesh->LocalNx * mesh->LocalNy;
#if BOUT_USE_METRIC_3D
  array_size *= mesh->LocalNz; // 3D metrics
#endif

  // Create the array and get the underlying data
  data = coords_store.emplace(coords, array_size).first->second.begin();

  // Copy data from Coordinates variable into data array
  // Uses the symbol to look up the corresponding Offset
#define COPY_STRIPE1(symbol)                                           \
  data[stripe_size * ind.ind + static_cast<int>(Offset::symbol)] = coords->symbol[ind];

  // Implement copy for each argument
#define COPY_STRIPE(...)                                \
  { MACRO_FOR_EACH(COPY_STRIPE1, __VA_ARGS__) }

  // Iterate over all points in the field
  // Note this could be 2D or 3D, depending on FieldMetric type
  for (const auto &ind : coords->dx.getRegion("RGN_ALL")) {
    COPY_STRIPE(dx, dy, dz);
    COPY_STRIPE(d1_dx, d1_dy, d1_dz);
    COPY_STRIPE(J);

    data[stripe_size * ind.ind + static_cast<int>(Offset::B)] = coords->Bxy[ind];
    data[stripe_size * ind.ind + static_cast<int>(Offset::Byup)] = coords->Bxy.yup()[ind];
    data[stripe_size * ind.ind + static_cast<int>(Offset::Bydown)] =
        coords->Bxy.ydown()[ind];

    COPY_STRIPE(G1, G3);
    COPY_STRIPE(g11, g12, g13, g22, g23, g33);
    COPY_STRIPE(g_11, g_12, g_13, g_22, g_23, g_33);
  }
}

std::size_t CoordinatesAccessor::clear(const Coordinates* coords) {
  if (coords == nullptr) {
    // clear all
    std::size_t num_removed = coords_store.size();
    coords_store.clear();
    return num_removed;
  }
  // Coordinates specified, so only remove one
  return coords_store.erase(coords);
}
