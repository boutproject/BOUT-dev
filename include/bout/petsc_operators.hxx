#pragma once

#ifndef BOUT_PETSC_OPERATORS
#define BOUT_PETSC_OPERATORS

#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/petsc_interface.hxx"
#include "bout/petsclib.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"

#include <petscmat.h>
#include <petscvec.h>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

/// @brief Tag type identifying the cell space C.
///
/// Used as a compile-time template parameter to distinguish vectors and operators
/// that live in the cell space from those in the forward or backward leg spaces.
/// Mismatched space tags produce a compile error rather than a runtime failure.
struct CellSpaceTag {};

/// @brief Tag type identifying the forward leg space L+.
///
/// Operators and vectors tagged with this type map to or from the forward (y+1)
/// leg space. See @ref BackwardLegSpaceTag and @ref CellSpaceTag.
struct ForwardLegSpaceTag {};

/// @brief Tag type identifying the backward leg space L-.
///
/// Operators and vectors tagged with this type map to or from the backward (y-1)
/// leg space. See @ref ForwardLegSpaceTag and @ref CellSpaceTag.
struct BackwardLegSpaceTag {};

template <typename Out, typename In>
class PetscOperator;

/// Used to indicate whether forward or backward boundary 'legs' should be used
enum class BoundaryDirection { Forward, Backward, Both };

class PetscCellMapping;

/// @brief Shared-pointer alias for a const PetscCellMapping.
using PetscCellMappingPtr = std::shared_ptr<const PetscCellMapping>;

/// @brief Bidirectional index mapping between mesh-file stored numbering and PETSc
///        row ownership.
///
/// Every vector space (cell space, forward leg space, backward leg space) is defined
/// by a global stored numbering written to the mesh file by the Python weights module.
/// PETSc distributes rows contiguously across MPI ranks, so the local PETSc row range
/// [rowStart(), rowEnd()) does not in general coincide with the stored indices.
///
/// This class maintains:
/// - the list of stored indices owned by the local rank, in PETSc row order;
/// - a hash map for O(1) stored → PETSc lookup;
/// - a permutation matrix @c mat_stored_to_petsc (PETSc ordering → stored numbering)
///   and its transpose @c mat_petsc_to_stored (stored numbering → PETSc ordering),
///   used to translate operator column indices during matrix assembly.
///
/// Subclasses specialise the constructor to populate these structures for a particular
/// space from mesh metadata.
class PetscIndexMapping {
public:
  PetscIndexMapping() = default;
  virtual ~PetscIndexMapping();

  PetscIndexMapping(const PetscIndexMapping&) = delete;
  PetscIndexMapping& operator=(const PetscIndexMapping&) = delete;
  PetscIndexMapping(PetscIndexMapping&&) = delete;
  PetscIndexMapping& operator=(PetscIndexMapping&&) = delete;

  /// @brief Number of rows owned locally by this MPI rank.
  PetscInt localSize() const {
    return static_cast<PetscInt>(local_stored_indices.size());
  }

  /// @brief Total number of rows in this space across all MPI ranks.
  PetscInt globalSize() const { return global_size; }

  /// @brief Global PETSc index of the first row owned by this MPI rank.
  PetscInt rowStart() const { return row_start; }

  /// @brief Global PETSc index one past the last row owned by this MPI rank.
  ///
  /// The locally owned rows are the half-open interval [rowStart(), rowEnd()).
  PetscInt rowEnd() const { return row_end; }

  /// @brief Stored mesh-file indices owned locally, listed in PETSc row order.
  ///
  /// Entry @c i corresponds to global PETSc row <tt>rowStart() + i</tt>.
  const std::vector<int>& localStoredIndices() const { return local_stored_indices; }

  /// @brief Convert a locally owned stored index to its global PETSc row index.
  ///
  /// @param stored_index A mesh-file stored index that must be owned by this rank.
  /// @returns The corresponding global PETSc row index.
  /// @throws BoutException if @p stored_index is not owned locally.
  PetscInt storedToPetsc(int stored_index) const;

  /// @brief Return the permutation matrix that maps PETSc row ordering to stored
  ///        mesh-file numbering.
  ///
  /// This is the transpose of @c mat_stored_to_petsc and is used in
  /// @ref PetscOperator::assembleFromCSR to translate stored column indices into
  /// PETSc column indices via a post-multiplication.
  Mat getPetscToStored() const { return mat_petsc_to_stored; }

protected:
  /// @brief BOUT++ PETSc library handle; ensures PETSc is initialised.
  PetscLib lib;

  /// @brief Construct the permutation matrices and populate internal lookup structures.
  ///
  /// Assigns contiguous PETSc rows starting from the rank-local offset determined by
  /// PETSc's ownership partitioning. Each stored index in @p stored_indices maps to
  /// one PETSc row, in the order given.
  ///
  /// @param nlocal         Number of rows owned by this rank.
  /// @param nglobal        Total number of rows across all ranks.
  /// @param stored_indices Mesh-file stored indices for the locally owned rows, in
  ///                       the order they should appear in the PETSc vector.
  void buildPermutation(PetscInt nlocal, PetscInt nglobal,
                        const std::vector<int>& stored_indices);

  PetscInt global_size{0}; ///< Total global row count.
  PetscInt row_start{0};   ///< First global PETSc row owned locally.
  PetscInt row_end{0};     ///< One-past-last global PETSc row owned locally.
  std::vector<int> local_stored_indices; ///< Stored indices in local PETSc row order.
  std::unordered_map<int, PetscInt> local_stored_to_petsc; ///< Stored → PETSc lookup.
  Mat mat_stored_to_petsc{nullptr}; ///< Permutation: stored numbering → PETSc ordering.
  Mat mat_petsc_to_stored{nullptr}; ///< Permutation: PETSc ordering → stored numbering.
};

/// @brief Cell-space index mapping that bridges stored cell numbers, PETSc row
///        ownership, and BOUT++ Field3D index arithmetic.
///
/// The cell space C contains five disjoint subsets of cells, appended in order:
///   1. Evolving interior cells (the primary unknowns).
///   2. Inner radial X-boundary cells.
///   3. Outer radial X-boundary cells.
///   4. Forward parallel boundary cells (one virtual cell per interior source cell
///      whose forward map intersects the radial boundary).
///   5. Backward parallel boundary cells (analogous for the backward map).
///
/// The local PETSc vector stores these in the same order. The @c mapLocal* methods
/// iterate over subsets 1–3 (interior + X-boundary) while @c mapLocalYup and
/// @c mapLocalYdown cover subsets 4 and 5 respectively.
///
/// @c mapOwnedInteriorCells iterates only the evolving interior cells (subset 1) and
/// provides global PETSc row indices, making it suitable for assembling matrices.
class PetscCellMapping : public PetscIndexMapping {
public:
  /// @brief Construct the cell-space mapping from mesh metadata fields.
  ///
  /// Populates all five cell subsets and calls @ref buildPermutation with the
  /// concatenated list of stored indices.
  ///
  /// @param cell_number          Field3D of stored cell numbers for interior and
  ///                             X-boundary cells. Negative entries are skipped.
  /// @param forward_cell_number  Field3D of stored cell numbers for forward parallel
  ///                             boundary virtual cells. Negative entries are skipped.
  /// @param backward_cell_number Field3D of stored cell numbers for backward parallel
  ///                             boundary virtual cells. Negative entries are skipped.
  /// @param total_cells          Total number of cells in the global cell space,
  ///                             as computed by the Python weights module.
  PetscCellMapping(const Field3D& cell_number, const Field3D& forward_cell_number,
                   const Field3D& backward_cell_number, int total_cells);

  /// @brief Build the region of evolving interior cells from a cell-number field.
  ///
  /// Iterates RGN_NOBNDRY and retains indices where @p cell_number is non-negative.
  ///
  /// @param cell_number Field3D of stored cell numbers.
  /// @returns Region containing the Ind3D indices of all valid interior cells.
  static Region<Ind3D> create_region(const Field3D& cell_number);

  /// @brief Build the region of inner radial X-boundary cells.
  ///
  /// Iterates RGN_INNER_X and retains indices where @p cell_number is non-negative.
  ///
  /// @param cell_number Field3D of stored cell numbers.
  /// @returns Region containing the Ind3D indices of valid inner-X boundary cells.
  static Region<Ind3D> create_region_xin(const Field3D& cell_number);

  /// @brief Build the region of outer radial X-boundary cells.
  ///
  /// Iterates RGN_OUTER_X and retains indices where @p cell_number is non-negative.
  ///
  /// @param cell_number Field3D of stored cell numbers.
  /// @returns Region containing the Ind3D indices of valid outer-X boundary cells.
  static Region<Ind3D> create_region_xout(const Field3D& cell_number);

  /// @brief Iterate over locally owned evolving interior cells, providing global
  ///        PETSc row indices.
  ///
  /// Calls @p func(petsc_row, field_index, stored_cell_number) for each cell in
  /// the evolving region that is owned by this MPI rank. @p petsc_row is the global
  /// PETSc index and is suitable for use as a row or column argument to MatSetValues.
  ///
  /// @tparam Function Callable with signature
  ///         <tt>(PetscInt row, Ind3D i, int stored) -> void</tt>.
  /// @param func Callback invoked for each owned interior cell.
  template <typename Function>
  void mapOwnedInteriorCells(Function func) const {
    PetscInt row = rowStart();
    BOUT_FOR_SERIAL(i, evolving_region) {
      func(row, i, ROUND(cell_number[i]));
      ++row;
    }
  }

  /// @brief Iterate over all locally stored field cells (interior + X-boundary),
  ///        providing local PETSc vector indices.
  ///
  /// Covers the evolving interior region, inner-X boundary, and outer-X boundary
  /// in that order. The local index starts at zero on each rank and increases
  /// contiguously. Use this to read from or write to a local PETSc Vec array
  /// obtained via VecGetArray.
  ///
  /// @tparam Function Callable with signature
  ///         <tt>(PetscInt local_row, Ind3D i) -> void</tt>.
  /// @param func Callback invoked for each locally stored field cell.
  template <typename Function>
  void mapLocalField(Function func) const {
    PetscInt row = 0;
    const std::vector<std::reference_wrapper<const Region<Ind3D>>> regions = {
        evolving_region, xin_region, xout_region};
    for (const auto& region : regions) {
      BOUT_FOR_SERIAL(i, region.get()) {
        func(row, i);
        ++row;
      }
    }
  }

  /// @brief Iterate over locally stored forward parallel boundary cells, providing
  ///        local PETSc vector indices.
  ///
  /// The local index starts immediately after the field cells covered by
  /// @ref mapLocalField. The Field3D index passed to @p func is the y+1 neighbour
  /// of the source cell (i.yp()), matching how field.yup() is indexed.
  ///
  /// @tparam Function Callable with signature
  ///         <tt>(PetscInt local_row, Ind3D i_yp) -> void</tt>.
  /// @param func Callback invoked for each locally stored forward boundary cell.
  template <typename Function>
  void mapLocalYup(Function func) const {
    PetscInt row = evolving_region.size() + xin_region.size() + xout_region.size();
    BOUT_FOR_SERIAL(i, yup_region) {
      func(row, i.yp());
      ++row;
    }
  }

  /// @brief Iterate over locally stored backward parallel boundary cells, providing
  ///        local PETSc vector indices.
  ///
  /// The local index starts immediately after the forward boundary cells covered by
  /// @ref mapLocalYup. The Field3D index passed to @p func is the y-1 neighbour
  /// of the source cell (i.ym()), matching how field.ydown() is indexed.
  ///
  /// @tparam Function Callable with signature
  ///         <tt>(PetscInt local_row, Ind3D i_ym) -> void</tt>.
  /// @param func Callback invoked for each locally stored backward boundary cell.
  template <typename Function>
  void mapLocalYdown(Function func) const {
    PetscInt row = evolving_region.size() + xin_region.size() + xout_region.size()
                   + yup_region.size();
    BOUT_FOR_SERIAL(i, ydown_region) {
      func(row, i.ym());
      ++row;
    }
  }

  friend PetscOperator<CellSpaceTag, CellSpaceTag>
  makeNeumannOperator(const PetscCellMappingPtr& mapping, BoundaryDirection direction);

private:
  Field3D cell_number; ///< Stored cell numbers for interior/X-boundary cells.
  Field3D
      forward_cell_number; ///< Stored cell numbers for forward boundary virtual cells.
  Field3D
      backward_cell_number; ///< Stored cell numbers for backward boundary virtual cells.
  Region<Ind3D> evolving_region; ///< Interior evolving cells.
  Region<Ind3D> xin_region;      ///< Inner radial X-boundary cells.
  Region<Ind3D> xout_region;     ///< Outer radial X-boundary cells.
  Region<Ind3D> yup_region;      ///< Forward parallel boundary virtual cells.
  Region<Ind3D> ydown_region;    ///< Backward parallel boundary virtual cells.
};

/// @brief Index mapping for forward (L+) or backward (L-) leg spaces.
///
/// Each interior source cell contributes one or two leg rows depending on whether
/// its mapped sub-samples are entirely interior, entirely boundary, or mixed.
/// The leg indices are assigned by the Python weights module and stored in the mesh
/// file as Field3D arrays (e.g. forward_leg_interior_number).
///
/// On construction, duplicate indices are removed and the remaining indices are
/// sorted before being passed to @ref buildPermutation.
class PetscLegMapping : public PetscIndexMapping {
public:
  /// @brief Construct a leg-space mapping.
  ///
  /// @param total_legs        Global number of leg rows in this space.
  /// @param local_leg_indices Stored leg indices owned by this MPI rank. Duplicates
  ///                          are removed internally before building the permutation.
  PetscLegMapping(int total_legs, std::vector<int> local_leg_indices);
};

/// @brief Shared-pointer alias for a const PetscLegMapping.
using PetscLegMappingPtr = std::shared_ptr<const PetscLegMapping>;

/// @brief Type-safe wrapper around a PETSc Vec belonging to a particular vector space.
///
/// The compile-time @p SpaceTag (one of @ref CellSpaceTag, @ref ForwardLegSpaceTag,
/// @ref BackwardLegSpaceTag) prevents accidental arithmetic between vectors from
/// incompatible spaces. All binary operations assert that both operands share the
/// same mapping pointer, catching dimension mismatches at runtime.
///
/// Vectors are non-copyable but movable, matching the ownership semantics of the
/// underlying PETSc Vec.
///
/// @tparam SpaceTag One of CellSpaceTag, ForwardLegSpaceTag, BackwardLegSpaceTag.
template <typename SpaceTag>
class PetscSpaceVector {
  using UniqueVec = bout::petsc::UniqueVec;

public:
  /// @brief Default constructor; produces an empty (unusable) vector.
  PetscSpaceVector() = default;

  /// @brief Construct a zero-initialised vector for the given space.
  ///
  /// Allocates a PETSc MPI vector sized according to @p mapping->localSize() and
  /// @p mapping->globalSize().
  ///
  /// @param mapping Shared ownership of the space mapping.
  explicit PetscSpaceVector(std::shared_ptr<const PetscIndexMapping> mapping)
      : mapping(std::move(mapping)),
        vec(createVec(this->mapping->localSize(), this->mapping->globalSize())) {}

  /// @brief Construct from an existing mapping and a pre-built PETSc Vec.
  ///
  /// Takes ownership of @p vec.
  ///
  /// @param mapping Shared ownership of the space mapping.
  /// @param vec     An already-assembled PETSc Vec. Ownership is transferred.
  PetscSpaceVector(std::shared_ptr<const PetscIndexMapping> mapping, UniqueVec vec)
      : mapping(std::move(mapping)), vec(std::move(vec)) {}

  PetscSpaceVector(const PetscSpaceVector&) = delete;
  PetscSpaceVector& operator=(const PetscSpaceVector&) = delete;
  PetscSpaceVector(PetscSpaceVector&&) noexcept = default;
  PetscSpaceVector& operator=(PetscSpaceVector&&) noexcept = default;

  /// @brief Return the underlying PETSc Vec handle (non-owning).
  Vec raw() const { return *this->vec; }

  /// @brief Return the shared space mapping.
  const std::shared_ptr<const PetscIndexMapping>& getMapping() const { return mapping; }

  /// @brief Return a deep copy of this vector.
  ///
  /// Allocates a new PETSc Vec with the same size and mapping, then copies values.
  PetscSpaceVector copy() const {
    UniqueVec out{new Vec{nullptr}};
    BOUT_DO_PETSC(VecDuplicate(*this->vec, out.get()));
    BOUT_DO_PETSC(VecCopy(*this->vec, *out.get()));
    return PetscSpaceVector(this->mapping, std::move(out));
  }

  /// @brief Return a new vector equal to @c (*this) * @p scalar.
  ///
  /// @param scalar Scalar multiplier.
  PetscSpaceVector operator*(BoutReal scalar) const& {
    auto out = this->copy();
    BOUT_DO_PETSC(VecScale(out.raw(), scalar));
    return out;
  }

  /// Multiply by a scalar, modifying a temporary in-place
  PetscSpaceVector operator*(BoutReal scalar) && {
    BOUT_DO_PETSC(VecScale(this->raw(), scalar));
    return std::move(*this);
  }

  /// @brief Return a new vector equal to @c (*this) + @p rhs.
  ///
  /// @param rhs Must share the same mapping as @c *this.
  PetscSpaceVector operator+(const PetscSpaceVector& rhs) const& {
    ASSERT0(mapping == rhs.mapping);
    auto out = this->copy();
    BOUT_DO_PETSC(VecAXPY(out.raw(), 1.0, rhs.raw()));
    return out;
  }

  PetscSpaceVector operator+(const PetscSpaceVector& rhs) && {
    ASSERT0(mapping == rhs.mapping);
    BOUT_DO_PETSC(VecAXPY(this->raw(), 1.0, rhs.raw()));
    return std::move(*this);
  }

  /// @brief Return a new vector equal to @c (*this) - @p rhs.
  ///
  /// @param rhs Must share the same mapping as @c *this.
  PetscSpaceVector operator-(const PetscSpaceVector& rhs) const& {
    ASSERT0(mapping == rhs.mapping);
    auto out = this->copy();
    BOUT_DO_PETSC(VecAXPY(out.raw(), -1.0, rhs.raw()));
    return out;
  }

  PetscSpaceVector operator-(const PetscSpaceVector& rhs) && {
    ASSERT0(mapping == rhs.mapping);
    BOUT_DO_PETSC(VecAXPY(this->raw(), -1.0, rhs.raw()));
    return std::move(*this);
  }

  /// @brief Return a new vector whose entries are the element-wise product
  ///        @c lhs[i] * rhs[i].
  ///
  /// @param lhs Left operand; must share a mapping with @p rhs.
  /// @param rhs Right operand.
  static PetscSpaceVector pointwiseMultiply(const PetscSpaceVector& lhs,
                                            const PetscSpaceVector& rhs) {
    ASSERT0(lhs.mapping == rhs.mapping);
    auto out = lhs.copy();
    BOUT_DO_PETSC(VecPointwiseMult(out.raw(), lhs.raw(), rhs.raw()));
    return out;
  }

  /// @brief Return a new vector whose entries are @c 1 / in[i].
  ///
  /// @param in Input vector. No entry may be zero.
  static PetscSpaceVector reciprocal(const PetscSpaceVector& in) {
    auto out = in.copy();
    BOUT_DO_PETSC(VecReciprocal(out.raw()));
    return out;
  }

private:
  /// @brief Allocate a new zeroed PETSc MPI vector of the given size.
  ///
  /// @param local_size  Number of locally owned entries.
  /// @param global_size Total number of entries across all ranks.
  static UniqueVec createVec(PetscInt local_size, PetscInt global_size) {
    UniqueVec out{new Vec{nullptr}};
    BOUT_DO_PETSC(VecCreate(BoutComm::get(), out.get()));
    BOUT_DO_PETSC(VecSetType(*out, VECMPI));
    BOUT_DO_PETSC(VecSetSizes(*out, local_size, global_size));
    return out;
  }

  std::shared_ptr<const PetscIndexMapping> mapping; ///< Space mapping (shared ownership).
  UniqueVec vec;                                    ///< Underlying PETSc Vec (owned).
};

// Tagged vector aliases — these provide compile-time checks of compatible operations.
/// @brief PETSc vector in the cell space C.
using PetscCellVector = PetscSpaceVector<CellSpaceTag>;
/// @brief PETSc vector in the forward leg space L+.
using PetscForwardLegVector = PetscSpaceVector<ForwardLegSpaceTag>;
/// @brief PETSc vector in the backward leg space L-.
using PetscBackwardLegVector = PetscSpaceVector<BackwardLegSpaceTag>;

/// @brief Populate a PetscCellVector from a Field3D and its parallel slices.
///
/// The cell vector is filled in the order defined by @ref PetscCellMapping:
/// evolving interior cells (from @p field), inner/outer X-boundary cells
/// (from @p field), forward boundary virtual cells (from @p field.yup()), and
/// backward boundary virtual cells (from @p field.ydown()).
///
/// @param mapping Shared cell-space mapping.
/// @param field   Source field; must have parallel slices split and allocated.
/// @returns A fully assembled PetscCellVector.
PetscCellVector makePetscCellVector(const PetscCellMappingPtr& mapping,
                                    const Field3D& field);

/// @brief Extract a Field3D (with parallel slices) from a PetscCellVector.
///
/// Reverses @ref makePetscCellVector: interior and X-boundary values are written to
/// @c result, forward boundary values to @c result.yup(), and backward boundary
/// values to @c result.ydown().
///
/// @param vec       Source vector; must use a PetscCellMapping.
/// @param prototype Field3D whose mesh, coordinates, and metadata are copied to
///                  initialise the result.
/// @returns A newly allocated Field3D with split parallel slices.
Field3D toField3D(const PetscCellVector& vec, const Field3D& prototype);

/// @brief Type-safe sparse linear operator mapping from one space to another.
///
/// Wraps a PETSc Mat and enforces compatible space pairings at compile time via
/// the @p OutSpaceTag and @p InSpaceTag template parameters. The underlying matrix
/// stores values in PETSc's internal ordering; translation from mesh-file stored
/// column indices to PETSc column indices is performed during construction via a
/// post-multiplication by @ref PetscIndexMapping::getPetscToStored.
///
/// Operators are non-copyable but movable. Binary operations (addition, subtraction,
/// scalar multiplication, composition) return new operators without modifying
/// their operands, unless the rvalue-qualified @c operator*(BoutReal) overload is
/// selected, which scales in place.
///
/// @tparam OutSpaceTag Tag for the output (row) space.
/// @tparam InSpaceTag  Tag for the input (column) space.
template <typename OutSpaceTag, typename InSpaceTag>
class PetscOperator {
  using UniqueMat = bout::petsc::UniqueMat;

public:
  /// @brief Default constructor; produces an empty (unusable) operator.
  PetscOperator() = default;
  ~PetscOperator() = default;

  /// @brief Construct an operator from a CSR-format sparse matrix stored in the mesh.
  ///
  /// Assembles the PETSc matrix from stored (non-PETSc) row and column indices by
  /// calling @ref assembleFromCSR.
  ///
  /// @param out_mapping Shared mapping for the output (row) space.
  /// @param in_mapping  Shared mapping for the input (column) space.
  /// @param rows        CSR row-pointer array (size = number_of_rows + 1).
  /// @param cols        CSR column-index array, using stored (mesh-file) indices.
  /// @param weights     CSR nonzero-value array.
  PetscOperator(std::shared_ptr<const PetscIndexMapping> out_mapping,
                std::shared_ptr<const PetscIndexMapping> in_mapping, Array<int> rows,
                Array<int> cols, Array<BoutReal> weights)
      : out_mapping(std::move(out_mapping)), in_mapping(std::move(in_mapping)) {
    assembleFromCSR(rows, cols, weights);
  }

  /// @brief Construct an operator from an already-assembled PETSc Mat.
  ///
  /// Takes ownership of @p mat.
  ///
  /// @param out_mapping Shared mapping for the output (row) space.
  /// @param in_mapping  Shared mapping for the input (column) space.
  /// @param mat         Pre-assembled PETSc Mat. Ownership is transferred.
  PetscOperator(std::shared_ptr<const PetscIndexMapping> out_mapping,
                std::shared_ptr<const PetscIndexMapping> in_mapping, UniqueMat mat)
      : out_mapping(std::move(out_mapping)), in_mapping(std::move(in_mapping)),
        mat_operator(std::move(mat)) {}

  // Operators are non-copyable (PETSc matrices are expensive to copy) but movable.
  PetscOperator(const PetscOperator&) = delete;
  PetscOperator& operator=(const PetscOperator&) = delete;
  PetscOperator(PetscOperator&&) noexcept = default;
  PetscOperator& operator=(PetscOperator&&) noexcept = default;

  /// @brief Apply the operator to a typed vector: @c result = A * rhs.
  ///
  /// @param rhs Input vector; must share the same mapping as this operator's
  ///            input space.
  /// @returns A new output-space vector holding the result.
  PetscSpaceVector<OutSpaceTag>
  operator()(const PetscSpaceVector<InSpaceTag>& rhs) const {
    ASSERT0(in_mapping == rhs.getMapping());
    PetscSpaceVector<OutSpaceTag> result(out_mapping);
    BOUT_DO_PETSC(MatMult(*this->mat_operator, rhs.raw(), result.raw()));
    return result;
  }

  /// @brief Apply the operator to a Field3D: converts to a cell vector then applies.
  ///
  /// Only available when @p InSpaceTag is @ref CellSpaceTag and @p OutSpaceTag is
  /// not @ref CellSpaceTag. The Field3D is packed into a PetscCellVector via
  /// @ref makePetscCellVector before applying the matrix.
  ///
  /// @param rhs Source field; must have parallel slices split and allocated.
  /// @returns A new output-space vector holding the result.
  template <typename T = InSpaceTag, typename U = OutSpaceTag,
            typename = std::enable_if_t<std::is_same_v<T, CellSpaceTag>
                                        && !std::is_same_v<U, CellSpaceTag>>>
  PetscSpaceVector<OutSpaceTag> operator()(const Field3D& rhs) const {
    auto in = makePetscCellVector(
        std::static_pointer_cast<const PetscCellMapping>(in_mapping), rhs);
    return (*this)(in);
  }

  /// Apply operator to a Field3D, returning a Field3D.
  /// Enable if both input and output mapppings are in Cell Space
  template <typename T = InSpaceTag, typename U = OutSpaceTag,
            typename = std::enable_if_t<std::is_same_v<T, CellSpaceTag>
                                        && std::is_same_v<U, CellSpaceTag>>>
  Field3D operator()(const Field3D& rhs) const {
    auto in = makePetscCellVector(
        std::static_pointer_cast<const PetscCellMapping>(in_mapping), rhs);
    return toField3D((*this)(in), rhs);
  }

  /// @brief Return the transpose operator mapping from OutSpaceTag to InSpaceTag.
  ///
  /// Allocates a new PETSc Mat via MatTranspose (MAT_INITIAL_MATRIX). The mappings
  /// are swapped: the transposed operator's input mapping is this operator's output
  /// mapping, and vice versa.
  PetscOperator<InSpaceTag, OutSpaceTag> transpose() const {
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatTranspose(*mat_operator, MAT_INITIAL_MATRIX, out.get()));
    return PetscOperator<InSpaceTag, OutSpaceTag>(in_mapping, out_mapping,
                                                  std::move(out));
  }

  /// @brief Return a new operator equal to @c (*this) + @p rhs.
  ///
  /// Both operands must share the same input and output mappings. Uses
  /// UNKNOWN_NONZERO_PATTERN, so the sparsity structure of the result may differ
  /// from either operand.
  ///
  /// @param rhs Operator to add; must have the same in/out mappings.
  PetscOperator operator+(const PetscOperator& rhs) const {
    ASSERT0(out_mapping == rhs.out_mapping);
    ASSERT0(in_mapping == rhs.in_mapping);
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatDuplicate(*this->mat_operator, MAT_COPY_VALUES, out.get()));
    BOUT_DO_PETSC(MatAXPY(*out, 1.0, *rhs.mat_operator, UNKNOWN_NONZERO_PATTERN));
    return PetscOperator(out_mapping, in_mapping, std::move(out));
  }

  /// @brief Return a new operator equal to @c (*this) - @p rhs.
  ///
  /// Both operands must share the same input and output mappings. Uses
  /// UNKNOWN_NONZERO_PATTERN; see @ref operator+.
  ///
  /// @param rhs Operator to subtract; must have the same in/out mappings.
  PetscOperator operator-(const PetscOperator& rhs) const {
    ASSERT0(out_mapping == rhs.out_mapping);
    ASSERT0(in_mapping == rhs.in_mapping);
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatDuplicate(*this->mat_operator, MAT_COPY_VALUES, out.get()));
    BOUT_DO_PETSC(MatAXPY(*out, -1.0, *rhs.mat_operator, UNKNOWN_NONZERO_PATTERN));
    return PetscOperator(out_mapping, in_mapping, std::move(out));
  }

  /// @brief Return a new operator equal to @c (*this) * @p scalar (lvalue overload).
  ///
  /// Duplicates the matrix before scaling, leaving @c *this unchanged.
  ///
  /// @param scalar Scalar multiplier.
  PetscOperator operator*(BoutReal scalar) const& {
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatDuplicate(*this->mat_operator, MAT_COPY_VALUES, out.get()));
    BOUT_DO_PETSC(MatScale(*out, scalar));
    return PetscOperator(out_mapping, in_mapping, std::move(out));
  }

  /// @brief Scale this operator in place and return it (rvalue overload).
  ///
  /// More efficient than the lvalue overload when @c *this is a temporary,
  /// as no matrix duplication is needed.
  ///
  /// @param scalar Scalar multiplier.
  PetscOperator operator*(BoutReal scalar) && {
    BOUT_DO_PETSC(MatScale(*this->mat_operator, scalar));
    return std::move(*this);
  }

  /// @brief Return the underlying PETSc Mat handle (non-owning).
  ///
  /// Intended for use by the free @c operator* (matrix composition) and other
  /// low-level callers that need direct PETSc access.
  Mat raw() const { return *this->mat_operator; }

  /// @brief Construct a diagonal operator from a vector of diagonal entries.
  ///
  /// Only available when @p OutSpaceTag == @p InSpaceTag (square, same-space
  /// operator). Builds an MPIAIJ matrix with one nonzero per row via
  /// MatDiagonalSet.
  ///
  /// @param mapping Shared mapping for both the row and column space.
  /// @param diag    Vector of diagonal values; must share @p mapping.
  /// @returns A new diagonal operator.
  static PetscOperator diagonal(std::shared_ptr<const PetscIndexMapping> mapping,
                                const PetscSpaceVector<OutSpaceTag>& diag) {
    static_assert(std::is_same_v<OutSpaceTag, InSpaceTag>,
                  "Diagonal operators require the same input and output space");
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatCreate(BoutComm::get(), out.get()));
    BOUT_DO_PETSC(MatSetType(*out, MATMPIAIJ));
    BOUT_DO_PETSC(MatSetSizes(*out, mapping->localSize(), mapping->localSize(),
                              mapping->globalSize(), mapping->globalSize()));
    // Note: Off-diagonal blocks are empty because this is a square matrix
    BOUT_DO_PETSC(MatMPIAIJSetPreallocation(*out, 1, nullptr, 0, nullptr));
    BOUT_DO_PETSC(MatDiagonalSet(*out, diag.raw(), INSERT_VALUES));
    BOUT_DO_PETSC(MatAssemblyBegin(*out, MAT_FINAL_ASSEMBLY));
    BOUT_DO_PETSC(MatAssemblyEnd(*out, MAT_FINAL_ASSEMBLY));
    return PetscOperator(mapping, mapping, std::move(out));
  }

private:
  /// @brief Assemble the operator matrix from a mesh-file CSR representation.
  ///
  /// Iterates over the locally owned output rows (identified via
  /// @c out_mapping->localStoredIndices()), looks up each row's nonzeros in the
  /// CSR arrays, and inserts them using global PETSc row indices. Column indices
  /// are given in stored (mesh-file) numbering; the final matrix is obtained by
  /// post-multiplying the assembled temporary by
  /// @c in_mapping->getPetscToStored() to convert to PETSc column ordering.
  ///
  /// @param rows    CSR row-pointer array (size = total_stored_rows + 1).
  /// @param cols    CSR column-index array in stored (mesh-file) numbering.
  /// @param weights CSR nonzero-value array.
  void assembleFromCSR(const Array<int>& rows, const Array<int>& cols,
                       const Array<BoutReal>& weights) {
    const UniqueMat temp{new Mat{nullptr}};
    BOUT_DO_PETSC(MatCreate(BoutComm::get(), temp.get()));
    BOUT_DO_PETSC(MatSetType(*temp, MATMPIAIJ));
    BOUT_DO_PETSC(MatSetSizes(*temp, out_mapping->localSize(), PETSC_DECIDE,
                              out_mapping->globalSize(), in_mapping->globalSize()));

    const auto& local_rows = out_mapping->localStoredIndices();
    for (PetscInt local_row = 0; local_row < static_cast<PetscInt>(local_rows.size());
         ++local_row) {
      const int stored_row = local_rows[local_row];
      if (stored_row < 0 || stored_row + 1 >= static_cast<int>(rows.size())) {
        continue;
      }
      const int start = rows[stored_row];
      const int end = rows[stored_row + 1];
      if (end <= start) {
        continue;
      }
      const PetscInt petsc_row = out_mapping->rowStart() + local_row;
      BOUT_DO_PETSC(MatSetValues(*temp, 1, &petsc_row, end - start, &cols[start],
                                 &weights[start], INSERT_VALUES));
    }
    BOUT_DO_PETSC(MatAssemblyBegin(*temp, MAT_FINAL_ASSEMBLY));
    BOUT_DO_PETSC(MatAssemblyEnd(*temp, MAT_FINAL_ASSEMBLY));

    // Post-multiply by the stored→PETSc permutation to convert column numbering.
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatMatMult(*temp, in_mapping->getPetscToStored(), MAT_INITIAL_MATRIX,
                             PETSC_DETERMINE, out.get()));
    mat_operator = std::move(out);
  }

  std::shared_ptr<const PetscIndexMapping> out_mapping; ///< Row-space mapping.
  std::shared_ptr<const PetscIndexMapping> in_mapping;  ///< Column-space mapping.
  UniqueMat mat_operator{new Mat{nullptr}};             ///< Underlying PETSc Mat.

  template <typename O, typename M, typename I>
  friend PetscOperator<O, I> operator*(const PetscOperator<O, M>& lhs,
                                       const PetscOperator<M, I>& rhs);
};

/// @brief Compose two compatible operators: @c result = lhs * rhs.
///
/// The output space of @p rhs must match the input space of @p lhs, enforced at
/// compile time by the @p MidSpaceTag parameter. The resulting operator maps from
/// @p rhs's input space to @p lhs's output space.
///
/// @tparam OutSpaceTag Output space of @p lhs and of the result.
/// @tparam MidSpaceTag Input space of @p lhs and output space of @p rhs.
/// @tparam InSpaceTag  Input space of @p rhs and of the result.
/// @param lhs Left (outer) operator.
/// @param rhs Right (inner) operator.
/// @returns A new operator representing the composition.
template <typename OutSpaceTag, typename MidSpaceTag, typename InSpaceTag>
PetscOperator<OutSpaceTag, InSpaceTag>
operator*(const PetscOperator<OutSpaceTag, MidSpaceTag>& lhs,
          const PetscOperator<MidSpaceTag, InSpaceTag>& rhs) {
  ASSERT0(lhs.in_mapping == rhs.out_mapping);
  typename PetscOperator<OutSpaceTag, InSpaceTag>::UniqueMat out{new Mat{nullptr}};
  BOUT_DO_PETSC(
      MatMatMult(lhs.raw(), rhs.raw(), MAT_INITIAL_MATRIX, PETSC_DETERMINE, out.get()));
  return PetscOperator<OutSpaceTag, InSpaceTag>(lhs.out_mapping, rhs.in_mapping,
                                                std::move(out));
}

/// @brief Operator mapping from cell space C to forward leg space L+.
using PetscForwardOperator = PetscOperator<ForwardLegSpaceTag, CellSpaceTag>;
/// @brief Operator mapping from cell space C to backward leg space L-.
using PetscBackwardOperator = PetscOperator<BackwardLegSpaceTag, CellSpaceTag>;
/// @brief Operator mapping from forward leg space L+ to cell space C.
using PetscForwardToCellOperator = PetscOperator<CellSpaceTag, ForwardLegSpaceTag>;
/// @brief Operator mapping from backward leg space L- to cell space C.
using PetscBackwardToCellOperator = PetscOperator<CellSpaceTag, BackwardLegSpaceTag>;
/// @brief Operator mapping within the cell space C.
using PetscCellOperator = PetscOperator<CellSpaceTag, CellSpaceTag>;
/// @brief Operator mapping within the forward leg space L+.
using PetscForwardLegOperator = PetscOperator<ForwardLegSpaceTag, ForwardLegSpaceTag>;
/// @brief Operator mapping within the backward leg space L-.
using PetscBackwardLegOperator = PetscOperator<BackwardLegSpaceTag, BackwardLegSpaceTag>;

/// @brief Factory that constructs cell-space and leg-space operators from mesh metadata.
///
/// On construction, reads all required fields and integer scalars from the mesh and
/// sets up the @ref PetscCellMapping, forward @ref PetscLegMapping, and backward
/// @ref PetscLegMapping. Individual primitive operators (interpolation, injection,
/// leg weights) are constructed on demand. The high-level set of parallel operators
/// needed for Support-Operator-Method diffusion is assembled by @ref getParallel.
class PetscOperators {
public:
  /// @brief Construct the operator factory from a BOUT++ mesh.
  ///
  /// Reads the following fields from @p mesh:
  ///   - @c cell_number, @c forward_cell_number, @c backward_cell_number (Field3D)
  ///   - @c total_cells, @c n_forward_legs, @c n_backward_legs (int)
  ///   - @c forward_leg_interior_number, @c forward_leg_boundary_number (Field3D)
  ///   - @c backward_leg_interior_number, @c backward_leg_boundary_number (Field3D)
  ///
  /// @param mesh Non-owning pointer to the BOUT++ mesh. Must outlive this object.
  /// @throws BoutException if any required field or scalar is missing from the mesh.
  explicit PetscOperators(Mesh* mesh);

  /// @brief Build the forward interpolation operator F: L+ ← C.
  ///
  /// Maps cell-space values to forward leg space using the CSR arrays stored in the
  /// mesh (@c forward_rows, @c forward_columns, @c forward_weights).
  PetscForwardOperator forward() const;

  /// @brief Build the backward interpolation operator B: L- ← C.
  ///
  /// Maps cell-space values to backward leg space using the CSR arrays stored in the
  /// mesh (@c backward_rows, @c backward_columns, @c backward_weights).
  PetscBackwardOperator backward() const;

  /// @brief Build the forward injection operator I+: L+ ← C.
  ///
  /// Places the value of each interior source cell into its corresponding forward leg
  /// row(s), with weight 1. Boundary leg rows point to the virtual boundary cell in C.
  PetscForwardOperator injectForward() const;

  /// @brief Build the backward injection operator I-: L- ← C.
  ///
  /// Analogous to @ref injectForward for the backward leg space.
  PetscBackwardOperator injectBackward() const;

  /// @brief Load the forward leg-weight vector w+: L+.
  ///
  /// Each entry holds the quadrature weight of the corresponding forward leg row,
  /// equal to the interior or boundary fraction of that leg's sub-sample integral.
  PetscForwardLegVector forwardLegWeights() const;

  /// @brief Load the backward leg-weight vector w-: L-.
  ///
  /// Analogous to @ref forwardLegWeights for the backward leg space.
  PetscBackwardLegVector backwardLegWeights() const;

  /// @brief Build a diagonal cell-space operator from a Field3D of values.
  ///
  /// @param f Field3D of diagonal entries; must have parallel slices allocated.
  /// @returns A diagonal PetscCellOperator.
  PetscCellOperator diagonal(const Field3D& f) const;

  /// @brief Build a diagonal forward-leg-space operator from a leg vector.
  ///
  /// @param v Forward leg vector of diagonal entries.
  /// @returns A diagonal PetscForwardLegOperator.
  PetscForwardLegOperator diagonal(const PetscForwardLegVector& v) const;

  /// @brief Build a diagonal backward-leg-space operator from a leg vector.
  ///
  /// @param v Backward leg vector of diagonal entries.
  /// @returns A diagonal PetscBackwardLegOperator.
  PetscBackwardLegOperator diagonal(const PetscBackwardLegVector& v) const;

  /// @brief Collection of parallel differential operators assembled for the
  ///        Support Operator Method (SOM).
  ///
  /// All operators are constructed together in @ref getParallel so that intermediate
  /// quantities (leg weights, interpolated dl, dV) are shared correctly.
  ///
  /// Naming convention: the subscript on Grad/Div/Restrict refers to which side of
  /// the cell face the operator acts on, not the leg space it lives in.
  ///   - @c Grad_plus / @c Div_plus: forward-side half-step operators (C to L+).
  ///   - @c Grad_minus / @c Div_minus: backward-side half-step operators (C to L-).
  ///   - @c Restrict_minus = I+^T * W+: weighted restriction from L+ back to C,
  ///     paired with the minus-side gradient.
  ///   - @c Restrict_plus  = I-^T * W-: weighted restriction from L- back to C,
  ///     paired with the plus-side gradient.
  struct Parallel {
    /// Central parallel gradient: C ← C.
    PetscCellOperator Grad_par;
    /// Support-operator parallel divergence paired with Grad_par: C ← C.
    PetscCellOperator Div_par;
    /// Combined SOM parallel Laplacian (Div_par ∘ Grad_par): C ← C.
    PetscCellOperator Div_par_Grad_par;
    /// Cell volume J * dx * dy * dz.
    Field3D dV;

    /// Half-step gradient on the backward (y-1) side: L- ← C.
    PetscBackwardOperator Grad_minus;
    /// Half-step gradient on the forward (y+1) side: L+ ← C.
    PetscForwardOperator Grad_plus;
    /// SOM divergence paired with Grad_plus: C ← L+.
    PetscForwardToCellOperator Div_minus;
    /// SOM divergence paired with Grad_minus: C ← L-.
    PetscBackwardToCellOperator Div_plus;

    /// Raw backward interpolation operator B: L- ← C.
    PetscBackwardOperator Backward;
    /// Raw forward interpolation operator F: L+ ← C.
    PetscForwardOperator Forward;

    /// Injection operator into backward leg space I-: L- ← C.
    PetscBackwardOperator Inject_minus;
    /// Injection operator into forward leg space I+: L+ ← C.
    PetscForwardOperator Inject_plus;

    /// Average interpolation to backward leg positions: L- ← C.
    PetscBackwardOperator Interp_minus;
    /// Average interpolation to forward leg positions: L+ ← C.
    PetscForwardOperator Interp_plus;

    /// Weighted restriction from L+ to C: I+^T * W+.
    PetscForwardToCellOperator Restrict_minus;
    /// Weighted restriction from L- to C: I-^T * W-.
    PetscBackwardToCellOperator Restrict_plus;

    /// @brief Compute the SOM parallel diffusion Div_par(K * Grad_par(f)).
    ///
    /// Evaluates the half-step fluxes separately on each side, multiplies by the
    /// interpolated diffusion coefficient K, and averages the two divergence
    /// contributions.
    ///
    /// @param K Diffusion coefficient field; must have parallel slices allocated.
    /// @param f Field to differentiate; must have parallel slices allocated.
    /// @returns Field3D containing Div_par(K * Grad_par(f)) at interior cells.
    Field3D Div_par_K_Grad_par(const Field3D& K, const Field3D& f) const;
  };

  /// @brief Assemble and return the full set of parallel SOM operators.
  ///
  /// Constructs all primitive operators (Forward, Backward, Inject_plus,
  /// Inject_minus, leg weights), derives the interpolation, gradient, restriction,
  /// and divergence operators, and returns them packaged in a @ref Parallel struct.
  /// Intermediate Field3D quantities (dl, dV) are computed from the mesh
  /// coordinates and boundary conditions are applied before use.
  ///
  /// @returns A @ref Parallel struct containing all assembled operators.
  Parallel getParallel() const;

private:
  /// @brief Read a Field3D from the mesh, throwing if the field is absent.
  ///
  /// @param mesh Non-owning mesh pointer.
  /// @param name Field name as stored in the mesh file.
  /// @returns The requested Field3D.
  /// @throws BoutException if the field is not found.
  static Field3D meshGetField3D(Mesh* mesh, const std::string& name);

  /// @brief Read an Array<T> from the mesh, throwing if absent.
  ///
  /// @tparam T Element type of the array.
  /// @param name Array name as stored in the mesh file.
  /// @returns The requested Array<T>.
  /// @throws BoutException if the array is not found.
  template <typename T>
  Array<T> meshGetArray(const std::string& name) const {
    Array<T> result;
    if (mesh->get(result, name) != 0) {
      throw BoutException("PetscOperators requires array '{}'", name);
    }
    return result;
  }

  /// @brief Read a scalar integer from the mesh, throwing if absent.
  ///
  /// @param name Scalar name as stored in the mesh file.
  /// @returns The requested integer value.
  /// @throws BoutException if the scalar is not found.
  int meshGetInt(const std::string& name) const;

  /// @brief Build a leg-space injection operator for the given direction.
  ///
  /// Iterates owned interior cells and inserts a unit entry for each interior or
  /// boundary leg row, mapping to the corresponding cell-space column. Column
  /// indices come from @p interior_leg_number and @p boundary_leg_number.
  ///
  /// @tparam LegTag     @ref ForwardLegSpaceTag or @ref BackwardLegSpaceTag.
  /// @param interior_leg_number Field3D of stored interior leg row numbers.
  /// @param boundary_leg_number Field3D of stored boundary leg row numbers.
  /// @param leg_mapping         Shared mapping for the leg space.
  /// @returns The assembled injection operator LegTag ← CellSpaceTag.
  template <typename LegTag>
  PetscOperator<LegTag, CellSpaceTag>
  buildInjection(const Field3D& interior_leg_number, const Field3D& boundary_leg_number,
                 std::shared_ptr<const PetscLegMapping> leg_mapping) const;

  /// @brief Load a leg-weight vector from a named mesh array.
  ///
  /// Reads the flat array @p name from the mesh, then scatters values into the
  /// local PETSc vector using the mapping's stored index list. Stored indices that
  /// are out of range are assigned weight 0.
  ///
  /// @tparam LegVector  @ref PetscForwardLegVector or @ref PetscBackwardLegVector.
  /// @param name    Array name in the mesh file (e.g. "forward_leg_weights").
  /// @param mapping Shared mapping for the leg space.
  /// @returns A fully populated leg-weight vector.
  template <typename LegVector>
  LegVector loadLegWeights(const std::string& name,
                           std::shared_ptr<const PetscLegMapping> mapping) const {
    auto values = meshGetArray<BoutReal>(name);
    LegVector out(mapping);
    PetscScalar* x{nullptr};
    BOUT_DO_PETSC(VecGetArray(out.raw(), &x));
    const auto& local = mapping->localStoredIndices();
    for (PetscInt i = 0; i < static_cast<PetscInt>(local.size()); ++i) {
      const int stored = local[i];
      x[i] = (stored >= 0 && stored < static_cast<int>(values.size())) ? values[stored]
                                                                       : 0.0;
    }
    BOUT_DO_PETSC(VecRestoreArray(out.raw(), &x));
    return out;
  }

  Mesh* mesh{nullptr};                     ///< Non-owning pointer to the BOUT++ mesh.
  PetscCellMappingPtr cell_mapping;        ///< Cell-space index mapping.
  PetscLegMappingPtr forward_leg_mapping;  ///< Forward leg-space index mapping.
  PetscLegMappingPtr backward_leg_mapping; ///< Backward leg-space index mapping.
  Field3D forward_leg_interior_number;     ///< Stored interior leg numbers for L+.
  Field3D forward_leg_boundary_number;     ///< Stored boundary leg numbers for L+.
  Field3D backward_leg_interior_number;    ///< Stored interior leg numbers for L-.
  Field3D backward_leg_boundary_number;    ///< Stored boundary leg numbers for L-.
};

/// @brief Build the Neumann boundary operator N : C → C.
///
/// For each cell that has a forward or backward boundary virtual cell,
/// copies the adjacent interior cell value into the virtual cell row.
/// All other rows are identity. Used to ensure Restrict_minus has non-zero
/// rows for boundary cells so that boundary conditions propagate into
/// Grad_par and Div_par.
///
/// @param direction   Which Y boundary cells will be set
///
/// @returns A PetscCellOperator representing the Neumann boundary mapping.
PetscCellOperator makeNeumannOperator(const PetscCellMappingPtr& mapping,
                                      BoundaryDirection direction);

#else

#warning PETSc not available. No PetscOperators.

#endif

#endif // BOUT_PETSC_OPERATORS
