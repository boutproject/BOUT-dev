#ifndef __INDEXOFFSET_H__
#define __INDEXOFFSET_H__

#include "bout/mesh.hxx"
#include "bout/region.hxx"

/// Helper class for offsetting Ind2D/Ind3D
///
/// Provides methods for offsetting by fixed amounts in x, y, z, as
/// well as a generic method for offsetting by any amount in multiple
/// directions
///
/// Also provides helper methods for converting Ind2D/Ind3D to x, y, z
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

  static_assert(std::is_base_of<Ind2D, T>::value || std::is_base_of<Ind3D, T>::value,
                "IndexOffset must be templated with either Ind2D or Ind3D");

  Mesh &mesh;
  const int nx, ny, nz;

  IndexOffset(Mesh &mesh)
      : mesh(mesh), nx(mesh.LocalNx), ny(mesh.LocalNy),
        nz(std::is_base_of<Ind3D, T>::value ? mesh.LocalNz : 1) {}

  /// Convenience functions for converting to (x, y, z)
  int x(T index) const { return (index.ind / nz) / ny; }
  int y(T index) const { return (index.ind / nz) % ny; }
  int z(T index) const { return (index.ind % nz); }

  /// Positive offset in x, default to offset by one
  const inline T xp(T index, int i = 1) const { return index + (i * ny * nz); }
  const inline T xm(T index, int i = 1) const { return index - (i * ny * nz); }
  /// Negative offset in x, default to offset by one
  /// Positive offset in y, default to offset by one
  const inline T yp(T index, int i = 1) const { return index + (i * nz); }
  const inline T ym(T index, int i = 1) const { return index - (i * nz); }
  /// Negative offset in y, default to offset by one
  /// Positive offset in z, default to offset by one.
  /// Periodic in z, cannot handle negative offsets
  const inline T zp(T index, int i = 1) const {
    ASSERT3(i > 0);
    if (nz == 1) {
      return index;
    }
    return (index + i) % nz < i ? index - nz + i : index + i;
  }
  /// Negative offset in z, default to offset by one.
  /// Periodic in z, cannot handle negative offsets
  const inline T zm(T index, int i = 1) const {
    ASSERT3(i > 0);
    if (nz == 1) {
      return index;
    }
    return (index) % nz < i ? index + nz - i : index - i;
  }

  /// Helper functions for offsets by two
  const inline T xpp(T index) const { return xp(index, 2); }
  const inline T xmm(T index) const { return xm(index, 2); }
  const inline T ypp(T index) const { return yp(index, 2); }
  const inline T ymm(T index) const { return ym(index, 2); }
  const inline T zpp(T index) const { return zp(index, 2); }
  const inline T zmm(T index) const { return zm(index, 2); }

  /// Generic offset in multiple direction simultaneously
  const inline T offset(T index, int dx, int dy, int dz) {
    auto temp = (dz > 0) ? zp(index, dz) : zm(index, -dz);
    return xp(yp(temp, dy), dx);
  }
};

#endif // __INDEXOFFSET_H__
