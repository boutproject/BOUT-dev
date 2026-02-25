#include "bout/output_bout_types.hxx"

#include "bout/mesh.hxx"
#include "bout/region.hxx"

#include <algorithm>

namespace bout::details {

template <class T>
auto region_transpose(const Region<T>& region) -> Region<T> {
  auto indices = region.getIndices();

  std::sort(indices.begin(), indices.end(), [](const T& lhs, const T& rhs) {
    const auto lx = lhs.x();
    const auto ly = lhs.y();
    const auto lz = lhs.z();

    const auto rx = rhs.x();
    const auto ry = rhs.y();
    const auto rz = rhs.z();

    // Z is now outer scale, so put it in largest blocks
    if (lz != rz) {
      return lz < rz;
    }
    if (ly != ry) {
      return ly < ry;
    }
    return lx < rx;
  });

  return Region<T>{indices};
}

template auto region_transpose(const Region<Ind2D>& region) -> Region<Ind2D>;
template auto region_transpose(const Region<Ind3D>& region) -> Region<Ind3D>;
template auto region_transpose(const Region<IndPerp>& region) -> Region<IndPerp>;
} // namespace bout::details
