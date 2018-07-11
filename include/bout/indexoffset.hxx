#ifndef __INDEXOFFSET_H__
#define __INDEXOFFSET_H__

#include "bout/mesh.hxx"
#include "bout/region.hxx"

/// Helper class for offsetting Ind2D/Ind3D/IndPerp
///
/// Provides methods for offsetting by fixed amounts in x, y, z, as
/// well as a generic method for offsetting by any amount in multiple
/// directions
///
/// Also provides helper methods for converting Ind2D/Ind3D/IndPerp to x, y, z
/// indices
///
/// Examples
/// --------
///
///     Field3D field, result;
///     auto index = std::begin(region);
///     IndexOffset<Ind3D> offset(*mesh);
///
///     result = field[offset.yp(*index)] - field[offset.ym(*index)];
///
template <typename T = Ind3D>
struct IndexOffset {

  static_assert(std::is_base_of<Ind2D, T>::value 
             || std::is_base_of<Ind3D, T>::value
             || std::is_base_of<IndPerp, T>::value,
                "IndexOffset must be templated with Ind2D, Ind3D or IndPerp");

  Mesh &mesh;
  const int nx, ny, nz;

  IndexOffset(Mesh &mesh)
      : mesh(mesh), nx(mesh.LocalNx), 
        ny(std::is_base_of<IndPerp,T>::value ? 1 : mesh.LocalNy),
        nz(std::is_base_of<Ind2D, T>::value ? 1 : mesh.LocalNz) {}

  /// Convenience functions for converting \p index to (x, y, z)
  int x(T index) const { return (index.ind / nz) / ny; }
  int y(T index) const { return (index.ind / nz) % ny; }
  int z(T index) const { return (index.ind % nz); }

  /// Positive offset \p index by \p dx in x, default to offset by one
  const inline T xp(T index, int dx = 1) const {
#if CHECK > 3
    if (x(index) + dx < 0 or x(index) + dx >= nx) {
      throw BoutException("Offset in x (%d) would go out of bounds at %d", dx, index.ind);
    }
#endif
    ASSERT2(std::abs(dx) < nx);
    return index + (dx * ny * nz);
  }
  /// Negative offset \p index by \p dx in x, default to offset by one
  const inline T xm(T index, int dx = 1) const { return xp(index, -dx); }
  /// Positive offset \p index by \p dy in y, default to offset by one
  const inline T yp(T index, int dy = 1) const {
#if CHECK > 3
    if (y(index) + dy < 0 or y(index) + dy >= ny) {
      throw BoutException("Offset in y (%d) would go out of bounds at %d", dy, index.ind);
    }
#endif
    ASSERT2(std::abs(dy) < ny);
    return index + (dy * nz);
  }
  /// Negative offset \p index by \p dy in y, default to offset by one
  const inline T ym(T index, int dy = 1) const { return yp(index, -dy); }
  /// Positive offset \p index by \p dz in z, default to offset by one.
  /// Periodic in z, cannot handle negative offsets
  const inline T zp(T index, int dz = 1) const {
    ASSERT2(dz > 0);
    ASSERT2(dz <= nz);
    return (index + dz) % nz < dz ? index - nz + dz : index + dz;
  }
  /// Negative offset \p index by \p dz in z, default to offset by one.
  /// Periodic in z, cannot handle negative offsets
  const inline T zm(T index, int dz = 1) const {
    ASSERT2(dz > 0);
    ASSERT2(dz <= nz);
    return (index) % nz < dz ? index + nz - dz : index - dz;
  }

  /// Helper functions for offsetting \p index by two
  const inline T xpp(T index) const { return xp(index, 2); }
  const inline T xmm(T index) const { return xm(index, 2); }
  const inline T ypp(T index) const { return yp(index, 2); }
  const inline T ymm(T index) const { return ym(index, 2); }
  const inline T zpp(T index) const { return zp(index, 2); }
  const inline T zmm(T index) const { return zm(index, 2); }

  /// Generic offset of \p index in multiple directions simultaneously
  const inline T offset(T index, int dx, int dy, int dz) {
    auto temp = (dz > 0) ? zp(index, dz) : zm(index, -dz);
    return xp(yp(temp, dy), dx);
  }
};

#endif // __INDEXOFFSET_H__
