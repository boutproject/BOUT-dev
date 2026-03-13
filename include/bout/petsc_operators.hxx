/// Represent operators using PETSc matrices
///

#include "bout/array.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/petsclib.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"

#include <petscmat.h>

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

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

  /// Loops over cells and calls the given function with arguments:
  ///   - PetscInt  PETSc row index
  ///   - Ind3D     Index into Field3D variables
  ///   - int       Global cell number in mesh
  ///
  template <typename Function>
  void map_evolving(Function func) const {
    PetscInt row = this->row_start;
    BOUT_FOR_SERIAL(i, this->evolving_region) {
      func(row, i, ROUND(cell_number[i]));
      ++row;
    }
  }

  /// Note: It doesn't make sense to pass the Ind3D because
  /// some reference Field3D, some the yup/ydown fields.
  template <typename Function>
  void map(Function func) const {
    const std::vector<std::reference_wrapper<const Region<Ind3D>>> regions = {
        evolving_region, xin_region, xout_region, yup_region, ydown_region};
    PetscInt row = row_start;
    for (const auto& region : regions) {
      BOUT_FOR_SERIAL(i, region.get()) {
        func(row, ROUND(cell_number[i]));
        ++row;
      }
    }
  }

  template <typename Function>
  void map_local_field(Function func) const {
    const std::vector<std::reference_wrapper<const Region<Ind3D>>> regions = {
        evolving_region, xin_region, xout_region};
    PetscInt row = 0; // Starting from 0, not row_start
    for (const auto& region : regions) {
      BOUT_FOR_SERIAL(i, region.get()) {
        func(row, i);
        ++row;
      }
    }
  }

  template <typename Function>
  void map_local_yup(Function func) const {
    PetscInt row = evolving_region.size() + xin_region.size() + xout_region.size();
    BOUT_FOR_SERIAL(i, yup_region) {
      func(row, i);
      ++row;
    }
  }

  template <typename Function>
  void map_local_ydown(Function func) const {
    PetscInt row = evolving_region.size() + xin_region.size() + xout_region.size()
                   + yup_region.size();
    BOUT_FOR_SERIAL(i, ydown_region) {
      func(row, i);
      ++row;
    }
  }

  /// Return a matrix that reorders vectors from petsc global index to mesh index
  Mat getPetscToMesh() const { return mat_petsc_to_mesh; }

private:
  PetscLib lib; // Initialize and finalize PETSc

  Field3D cell_number;

  Region<Ind3D> evolving_region; ///< Evolving cells
  Region<Ind3D> xin_region;      ///< X boundary cells in inner boundary
  Region<Ind3D> xout_region;     ///< X boundary cells in outer boundary
  Region<Ind3D> yup_region;      ///< Y boundary cells in forward direction
  Region<Ind3D> ydown_region;    ///< Y boundary cells in backward direction

  PetscInt row_start, row_end; ///< Local row indices
  Mat mat_mesh_to_petsc;
  Mat mat_petsc_to_mesh;
};

/// Shared pointer to const PetscMapping
/// Mapping is const after creation.
using PetscMappingPtr = std::shared_ptr<const PetscMapping>;

/// Represents an operator on Field3Ds
class PetscOperator {
public:
  /// Create using CSR format arrays
  PetscOperator(PetscMappingPtr mapping, Array<int> rows, Array<int> cols,
                Array<BoutReal> weights);

  /// Perform operation
  Field3D operator()(const Field3D& rhs) const;

  /// Operator composition
  PetscOperator operator*(const PetscOperator& rhs) const;

  /// Operator addition
  PetscOperator operator+(const PetscOperator& rhs) const;

  /// Operator subtraction
  PetscOperator operator-(const PetscOperator& rhs) const;

private:
  PetscOperator(PetscMappingPtr mapping, Mat mat)
      : mapping(std::move(mapping)), mat_operator(mat) {}
  PetscMappingPtr mapping;
  Mat mat_operator;
};

class PetscOperators {
public:
  ///
  /// Reads from the mesh:
  /// - cell_number : Field3D
  /// - forward_cell_number : Field3D
  /// - backward_cell_number : Field3D
  /// - total_cells : int, optional
  PetscOperators(Mesh* mesh);

private:
  /// Read a Field3D from the mesh or throw BoutException if not found.
  Field3D meshGetField3D(Mesh* mesh, const std::string& name);

  /// Read a 1D Array<T> from the mesh or throw BoutException
  template <typename T>
  Array<T> meshGetArray(Mesh* mesh, const std::string& name) {
    Array<T> result;
    if (mesh->get(result, name) != 0) {
      throw BoutException("PetscOperators requires Array<int> '{}'", name);
    }
    return result;
  }

  PetscMappingPtr mapping;
};
