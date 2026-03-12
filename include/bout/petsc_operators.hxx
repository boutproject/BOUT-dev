/// Represent operators using PETSc matrices
///

#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/petsclib.hxx"
#include "bout/region.hxx"

#include <petscmat.h>

#include <string>

/// Handles mapping between global cell index used in the mesh file,
/// and the global indexing used in PETSc that depends on the mesh
/// decomposition across MPI ranks.
///
/// Each processor reads a chunk of the global arrays into Field3Ds:
///  - cell_number
///  - forward_cell_number
///  - backward_cell_number
/// Non-negative values in these arrays indicate a cell number.
///
/// A Petsc Mat is constructed that maps between petsc indices and
/// mesh indices.
class PetscMapping {
public:
  PetscMapping(const Field3D& cell_number, const Field3D& forward_cell_number,
               const Field3D& backward_cell_number);

  // Static functions for unit testing
  static Region<Ind3D> create_region(const Field3D& cell_number);
  static Region<Ind3D> create_region_xin(const Field3D& cell_number);
  static Region<Ind3D> create_region_xout(const Field3D& cell_number);

  /// Total number of cells, including X and Y boundaries
  unsigned int size() const {
    return evolving_region.size() + xin_region.size() + xout_region.size()
           + yup_region.size() + ydown_region.size();
  }

private:
  PetscLib lib; // Initialize and finalize PETSc

  Region<Ind3D> evolving_region; ///< Evolving cells
  Region<Ind3D> xin_region;      ///< X boundary cells in inner boundary
  Region<Ind3D> xout_region;     ///< X boundary cells in outer boundary
  Region<Ind3D> yup_region;      ///< Y boundary cells in forward direction
  Region<Ind3D> ydown_region;    ///< Y boundary cells in backward direction

  Mat mat_mesh_to_petsc;
};

class PetscOperators {
public:
  PetscOperators(Mesh* mesh);

private:
  /// Read a Field3D from the mesh or throw BoutException if not found.
  Field3D meshGetField3D(Mesh* mesh, const std::string& name);

  PetscMapping mapping;
};
