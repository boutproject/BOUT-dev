/// Test communication of Z guard cells, including the X-Z corner cells
///
/// Every cell of a Field3D is given a value that depends only on its global
/// (x, z) index. X boundary cells are set only for the interior Z points, as the
/// boundary conditions do. After `communicate`, every cell including the Z guard
/// cells and the X-Z corners must hold the value of the periodic image of that
/// cell. Run with different NXPE, NZPE and MZG.

#include "bout/bout.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"

#include <cmath>

namespace {
/// Value of a cell with global X index \p xg and global (no-guard) Z index \p zg
BoutReal expected(const Mesh& mesh, int x, int zg) {
  const int xg = mesh.getGlobalXIndex(x);
  const BoutReal value = (1000. * xg) + zg;
  // X boundary cells get a different (negative) value, as if set by a boundary
  // condition applied on the processor that owns them
  const bool x_boundary = (mesh.firstX() and x < mesh.xstart)
                          or (mesh.lastX() and x > mesh.xend);
  return x_boundary ? -value - 0.5 : value;
}

/// Global Z index (no guards) of local index \p z, including the periodic wrap
int globalZ(const Mesh& mesh, int z) {
  const int nz = mesh.GlobalNzNoBoundaries;
  return (((mesh.getGlobalZIndexNoBoundaries(z)) % nz) + nz) % nz;
}
} // namespace

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  using bout::globals::mesh;

  Field3D field{-1.0, mesh}; // Guards are -1 until communicated
  for (int x = 0; x < mesh->LocalNx; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = mesh->zstart; z <= mesh->zend; ++z) {
        const bool x_guard = (x < mesh->xstart and not mesh->firstX())
                             or (x > mesh->xend and not mesh->lastX());
        if (not x_guard) {
          field(x, y, z) = expected(*mesh, x, globalZ(*mesh, z));
        }
      }
    }
  }

  mesh->communicate(field);

  int failures = 0;
  for (int x = 0; x < mesh->LocalNx; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        const BoutReal want = expected(*mesh, x, globalZ(*mesh, z));
        if (field(x, y, z) != want) {
          if (failures < 20) {
            output.write("Rank {}: ({}, {}, {}) is {} but expected {}\n",
                         BoutComm::rank(), x, y, z, field(x, y, z), want);
          }
          ++failures;
        }
      }
    }
  }

  int total_failures = 0;
  MPI_Allreduce(&failures, &total_failures, 1, MPI_INT, MPI_SUM, BoutComm::get());
  output.write("NXPE={} NZPE={} MZG={}: {} wrong cells\n", mesh->getNXPE(),
               mesh->getNZPE(), mesh->zstart, total_failures);
  if (total_failures == 0) {
    output.write("=> Z communication test passed\n");
  }

  BoutFinalise();
  return total_failures == 0 ? 0 : 1;
}
