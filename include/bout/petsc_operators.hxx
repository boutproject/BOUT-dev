/// Represent finite-difference operators using PETSc matrices.
///
/// This header defines:
/// - `PetscMapping`, which maps between mesh-global numbering,
///   PETSc numbering, and `Field3D` storage;
/// - `PetscOperator`, which wraps a PETSc `Mat` as a linear operator on
///   `Field3D` objects; and
/// - `PetscOperators`, a helper for reading mappings and operators from
///   a `Mesh`.
///

#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/petsc_interface.hxx"
#include "bout/petsclib.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"

#include <petscmat.h>

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

/// Map between mesh-global cell numbering and PETSc numbering.
///
/// The mesh file stores operator stencils and weights using a mesh-global
/// numbering scheme. This numbering includes:
/// - evolving cells in the main domain;
/// - X-boundary cells; and
/// - Y-boundary cells represented by the `yup` and `ydown` fields.
///
/// PETSc uses a different global numbering based on the MPI distribution
/// of rows and vector entries. `PetscMapping` records the correspondence
/// between these numbering systems and provides helper methods to iterate
/// over the local cells in a consistent order.
///
/// The mapping order used by this class is:
/// 1. evolving cells;
/// 2. inner X-boundary cells;
/// 3. outer X-boundary cells;
/// 4. `yup` boundary cells; and
/// 5. `ydown` boundary cells.
///
/// A pair of PETSc matrices is constructed internally to reorder vectors
/// between mesh-global and PETSc-global numbering.
class PetscMapping {
public:
  /// Construct the mapping from mesh numbering fields.
  ///
  /// @param cell_number Mesh-global cell numbers for evolving cells.
  /// @param forward_cell_number Mesh-global cell numbers associated with
  ///   forward Y-boundary (`yup`) values.
  /// @param backward_cell_number Mesh-global cell numbers associated with
  ///   backward Y-boundary (`ydown`) values.
  ///
  /// Non-negative values in the numbering fields identify valid cells.
  /// The constructor builds the local regions and the PETSc permutation
  /// matrices needed to translate between mesh and PETSc numbering.
  PetscMapping(const Field3D& cell_number, const Field3D& forward_cell_number,
               const Field3D& backward_cell_number);

  /// Destroy the PETSc matrices owned by this mapping.
  ~PetscMapping();

  /// Build a region containing all valid cells in a numbering field.
  ///
  /// @param cell_number Field whose entries contain mesh-global cell numbers.
  /// @return A region containing all indices with non-negative cell numbers.
  ///
  /// This helper is public for unit testing.
  static Region<Ind3D> create_region(const Field3D& cell_number);

  /// Build a region containing valid cells on the inner X boundary.
  ///
  /// @param cell_number Field whose entries contain mesh-global cell numbers.
  /// @return A region containing valid cells adjacent to the inner X boundary.
  static Region<Ind3D> create_region_xin(const Field3D& cell_number);

  /// Build a region containing valid cells on the outer X boundary.
  ///
  /// @param cell_number Field whose entries contain mesh-global cell numbers.
  /// @return A region containing valid cells adjacent to the outer X boundary.
  static Region<Ind3D> create_region_xout(const Field3D& cell_number);

  /// Return the total number of locally mapped cells.
  ///
  /// The count includes evolving cells and all boundary cells represented
  /// in the mapping on this processor.
  unsigned int localSize() const {
    return evolving_region.size() + xin_region.size() + xout_region.size()
           + yup_region.size() + ydown_region.size();
  }

  unsigned int globalSize() const {
    // Get global size
    PetscInt globalRows{0};
    PetscInt globalCols{0};
    BOUT_DO_PETSC(MatGetSize(this->mat_mesh_to_petsc, &globalRows, &globalCols));
    ASSERT1(globalRows == globalCols);
    return globalRows;
  }

  /// Iterate over locally owned evolving cells in PETSc row order.
  ///
  /// @tparam Function Callable with signature
  ///   `func(PetscInt petsc_row, Ind3D index, int mesh_global_index)`.
  /// @param func Function called once for each evolving cell owned locally.
  ///
  /// The PETSc row index starts at the locally owned global row `row_start`.
  /// The supplied `Ind3D` index refers to the main `Field3D` storage.
  template <typename Function>
  void map_evolving(Function func) const {
    PetscInt row = this->row_start;
    BOUT_FOR_SERIAL(i, this->evolving_region) {
      func(row, i, ROUND(cell_number[i]));
      ++row;
    }
  }

  /// Iterate over all locally mapped cells in PETSc row order.
  ///
  /// @tparam Function Callable with signature
  ///   `func(PetscInt petsc_row, int mesh_global_index)`.
  /// @param func Function called once for each mapped cell owned locally.
  ///
  /// The iteration order is evolving, inner X boundary, outer X boundary,
  /// `yup`, then `ydown`.
  ///
  /// No `Ind3D` is passed because some entries refer to interior `Field3D`
  /// storage while others refer to boundary storage in `yup` or `ydown`.
  template <typename Function>
  void map(Function func) const {
    const std::vector<std::reference_wrapper<const Region<Ind3D>>> regions = {
        evolving_region, xin_region, xout_region};
    PetscInt row = row_start;
    for (const auto& region : regions) {
      BOUT_FOR_SERIAL(i, region.get()) {
        func(row, ROUND(cell_number[i]));
        ++row;
      }
    }
    BOUT_FOR_SERIAL(i, yup_region) {
      func(row, ROUND(forward_cell_number[i]));
      ++row;
    }
    BOUT_FOR_SERIAL(i, ydown_region) {
      func(row, ROUND(backward_cell_number[i]));
      ++row;
    }
    ASSERT0(row == row_end);
  }

  /// Iterate over locally stored entries packed from the main `Field3D`.
  ///
  /// @tparam Function Callable with signature `func(PetscInt local_index, Ind3D index)`.
  /// @param func Function called once for each evolving or X-boundary cell
  ///   stored in the local PETSc vector layout.
  ///
  /// The index passed to `func` is a local vector index beginning at zero.
  /// Only the main field storage is traversed here; Y-boundary values are
  /// handled separately by `map_local_yup()` and `map_local_ydown()`.
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

  /// Iterate over local vector entries corresponding to `yup` boundary values.
  ///
  /// @tparam Function Callable with signature `func(PetscInt local_index, Ind3D index)`.
  /// @param func Function called once for each locally stored `yup` entry.
  ///
  /// The local index is the offset into the packed local PETSc vector.
  template <typename Function>
  void map_local_yup(Function func) const {
    PetscInt row = evolving_region.size() + xin_region.size() + xout_region.size();
    BOUT_FOR_SERIAL(i, yup_region) {
      func(row, i.yp());
      ++row;
    }
  }

  /// Iterate over local vector entries corresponding to `ydown` boundary values.
  ///
  /// @tparam Function Callable with signature `func(PetscInt local_index, Ind3D index)`.
  /// @param func Function called once for each locally stored `ydown` entry.
  ///
  /// The local index is the offset into the packed local PETSc vector.
  template <typename Function>
  void map_local_ydown(Function func) const {
    PetscInt row = evolving_region.size() + xin_region.size() + xout_region.size()
                   + yup_region.size();
    BOUT_FOR_SERIAL(i, ydown_region) {
      func(row, i.ym());
      ++row;
    }
  }

  /// Return the PETSc matrix that maps PETSc ordering to mesh ordering.
  ///
  /// @return PETSc matrix representing the permutation from PETSc global
  ///   indices to mesh-global indices.
  Mat getPetscToMesh() const { return mat_petsc_to_mesh; }

private:
  PetscLib lib; // Initialize and finalize PETSc

  /// Mesh-global numbering for cells stored on this rank.
  Field3D cell_number;
  Field3D forward_cell_number;
  Field3D backward_cell_number;

  /// Evolving cells in the main domain.
  Region<Ind3D> evolving_region;

  /// X-boundary cells on the inner radial boundary.
  Region<Ind3D> xin_region;

  /// X-boundary cells on the outer radial boundary.
  Region<Ind3D> xout_region;

  /// Y-boundary cells in the forward direction.
  Region<Ind3D> yup_region;

  /// Y-boundary cells in the backward direction.
  Region<Ind3D> ydown_region;

  /// First and one-past-last PETSc global rows owned by this rank.
  PetscInt row_start, row_end;

  /// Permutation matrix mapping mesh-global ordering to PETSc ordering.
  Mat mat_mesh_to_petsc;

  /// Permutation matrix mapping PETSc ordering to mesh-global ordering.
  Mat mat_petsc_to_mesh;
};

/// Shared pointer to an immutable `PetscMapping`.
///
/// The mapping is constructed once and then shared by operators that use it.
using PetscMappingPtr = std::shared_ptr<const PetscMapping>;

/// Linear operator on `Field3D` backed by a PETSc matrix.
///
/// `PetscOperator` wraps a PETSc `Mat` whose rows and columns use the PETSc
/// numbering defined by a corresponding `PetscMapping`. Operators are usually
/// constructed from CSR arrays read from the mesh file and can then be:
/// - applied to a `Field3D`;
/// - composed with other operators; and
/// - added or subtracted.
class PetscOperator {
  using UniqueVec = bout::petsc::UniqueVec; // unique_ptr to Vec
  using UniqueMat = bout::petsc::UniqueMat; // unique_ptr to Mat

public:
  /// Construct an operator from CSR data in mesh-global numbering.
  ///
  /// @param mapping Shared mapping between mesh and PETSc numbering.
  /// @param rows CSR row offsets.
  /// @param cols CSR column indices in mesh-global numbering.
  /// @param weights CSR non-zero values.
  ///
  /// The CSR arrays define an operator in mesh-global numbering. The
  /// constructor converts this representation into a PETSc matrix using
  /// the supplied mapping.
  PetscOperator(PetscMappingPtr mapping, Array<int> rows, Array<int> cols,
                Array<BoutReal> weights);

  // No copying
  PetscOperator(const PetscOperator&) = delete;
  PetscOperator& operator=(const PetscOperator&) = delete;

  // Allow Moving
  PetscOperator(PetscOperator&&) noexcept = default;
  PetscOperator& operator=(PetscOperator&&) noexcept = default;

  /// Destroy the PETSc matrix and working vectors owned by this operator.
  ~PetscOperator() = default;

  /// Create a diagonal operator
  static PetscOperator diagonal(PetscMappingPtr mapping, const Field3D& f);

  /// Apply the operator to a field.
  ///
  /// @param rhs Input field, including any required `yup` and `ydown`
  ///   boundary values.
  /// @return Result of applying the operator to `rhs`.
  ///
  /// The input field is packed into PETSc vector storage using the mapping,
  /// multiplied by the PETSc matrix, then unpacked back into a `Field3D`.
  Field3D operator()(const Field3D& rhs) const;

  /// Compose this operator with another operator.
  ///
  /// @param rhs Right-hand-side operator.
  /// @return The composed operator `(*this) * rhs`.
  ///
  /// Both operators must use compatible mappings.
  PetscOperator operator*(const PetscOperator& rhs) const;

  /// Add two operators.
  ///
  /// @param rhs Right-hand-side operator.
  /// @return The sum of the two operators.
  ///
  /// Both operators must use compatible mappings.
  PetscOperator operator+(const PetscOperator& rhs) const;

  /// Subtract one operator from another.
  ///
  /// @param rhs Right-hand-side operator.
  /// @return The difference of the two operators.
  ///
  /// Both operators must use compatible mappings.
  PetscOperator operator-(const PetscOperator& rhs) const;

  /// Calculate transpose as a new matrix
  PetscOperator transpose() const;

  void view() const { MatView(*this->mat_operator, PETSC_VIEWER_STDOUT_WORLD); }

  /// Multiply by a scalar
  /// This version allocates a new Mat and so is for lvalues.
  PetscOperator operator*(BoutReal scalar) const& {
    UniqueMat mat{new Mat{}};
    BOUT_DO_PETSC(MatDuplicate(*mat_operator, MAT_COPY_VALUES, mat.get()));
    BOUT_DO_PETSC(MatScale(*mat, scalar));
    return PetscOperator(this->mapping, std::move(mat));
  }

  /// Multiply by a scalar.
  /// This version is for rvalue temporaries and modifies in-place.
  PetscOperator operator*(BoutReal scalar) && {
    BOUT_DO_PETSC(MatScale(*this->mat_operator, scalar));
    return std::move(*this);
  }

private:
  /// Construct directly from an existing PETSc matrix.
  ///
  /// @param mapping Shared mapping between mesh and PETSc numbering.
  /// @param mat PETSc matrix implementing the operator.
  ///
  /// This constructor is used internally when creating operators from
  /// PETSc matrix algebra such as composition, addition, or subtraction.
  PetscOperator(PetscMappingPtr mapping, UniqueMat mat)
      : mapping(std::move(mapping)), mat_operator(std::move(mat)),
        rhs_vec(createVec(this->mapping->localSize())),
        result_vec(createVec(this->mapping->localSize())) {}

  /// Shared mapping between mesh and PETSc numbering.
  PetscMappingPtr mapping;

  /// PETSc matrix implementing the operator.
  UniqueMat mat_operator;

  /// Work vector holding the packed input field.
  UniqueVec rhs_vec;

  /// Work vector holding the packed result.
  UniqueVec result_vec;

  /// Copy contents of f, f.yup() and f.ydown() into vec.
  /// Assumes that vec has already been allocated.
  static void copyToVec(PetscMappingPtr mapping, const Field3D& f, Vec vec);

  // Allocate MPI vector
  static UniqueVec createVec(PetscInt local_size);
};

/// Collection of PETSc operators read from a mesh.
///
/// `PetscOperators` owns the `PetscMapping` derived from the mesh numbering
/// fields and provides access to named operators stored in the mesh file
/// as CSR arrays.
class PetscOperators {
public:
  /// Construct PETSc operator support from a mesh.
  ///
  /// Reads the following mesh variables:
  /// - `cell_number`
  /// - `forward_cell_number`
  /// - `backward_cell_number`
  ///
  /// and uses them to construct the shared `PetscMapping`.
  ///
  /// @param mesh Mesh from which numbering fields and operators are read.
  PetscOperators(Mesh* mesh);

  /// Retrieve a named PETSc operator from the mesh.
  ///
  /// @param name Base name of the operator in the mesh file.
  /// @return A `PetscOperator` constructed from the arrays
  ///   `{name}_rows`, `{name}_columns`, and `{name}_weights`.
  ///
  /// Throws `BoutException` if any of the required arrays are missing.
  PetscOperator get(const std::string& name) const;

  /// Create a diagonal operator
  PetscOperator diagonal(const Field3D& f) const {
    return PetscOperator::diagonal(mapping, f);
  }

private:
  /// Read a `Field3D` from the mesh.
  ///
  /// @param mesh Mesh from which to read.
  /// @param name Name of the field.
  /// @return The requested `Field3D`.
  ///
  /// Throws `BoutException` if the field is not present.
  static Field3D meshGetField3D(Mesh* mesh, const std::string& name);

  /// Read a one-dimensional `Array<T>` from the mesh.
  ///
  /// @tparam T Array element type.
  /// @param mesh Mesh from which to read.
  /// @param name Name of the array.
  /// @return The requested array.
  ///
  /// Throws `BoutException` if the array is not present.
  template <typename T>
  Array<T> meshGetArray(Mesh* mesh, const std::string& name) const {
    Array<T> result;
    if (mesh->get(result, name) != 0) {
      throw BoutException("PetscOperators requires Array<int> '{}'", name);
    }
    return result;
  }

  /// Mesh from which operator data are read.
  Mesh* mesh;

  /// Shared mapping used by all operators created from this mesh.
  PetscMappingPtr mapping;
};

#else // BOUT_HAS_PETSC

#warning PETSc not available. No PetscOperators.

#endif // BOUT_HAS_PETSC
