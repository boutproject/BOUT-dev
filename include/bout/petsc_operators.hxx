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
#include <petscvec.h>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

/// Tag type for the cell space C.
struct CellSpaceTag {};

/// Tag type for the forward leg space L+.
struct ForwardLegSpaceTag {};

/// Tag type for the backward leg space L-.
struct BackwardLegSpaceTag {};

/// Mapping between stored indices in a mesh file and PETSc row ownership.
///
/// Every space is defined by a global stored numbering written to the mesh file.
/// PETSc distributes rows contiguously across MPI ranks, so this class stores the
/// local owned stored indices in row order and a permutation matrix that maps
/// from stored numbering to PETSc numbering.
class PetscIndexMapping {
public:
  PetscIndexMapping() = default;
  virtual ~PetscIndexMapping();

  PetscIndexMapping(const PetscIndexMapping&) = delete;
  PetscIndexMapping& operator=(const PetscIndexMapping&) = delete;
  PetscIndexMapping(PetscIndexMapping&&) = delete;
  PetscIndexMapping& operator=(PetscIndexMapping&&) = delete;

  /// Number of rows owned locally by this MPI rank.
  PetscInt localSize() const {
    return static_cast<PetscInt>(local_stored_indices.size());
  }

  /// Global number of rows in this space.
  PetscInt globalSize() const { return global_size; }

  /// First global PETSc row owned by this MPI rank.
  PetscInt rowStart() const { return row_start; }

  /// One-past-last global PETSc row owned by this MPI rank.
  PetscInt rowEnd() const { return row_end; }

  /// Stored indices owned locally in PETSc row order.
  const std::vector<int>& localStoredIndices() const { return local_stored_indices; }

  /// Map a locally owned stored index to a PETSc global row.
  PetscInt storedToPetsc(int stored_index) const;

  /// Permutation matrix mapping PETSc ordering to stored numbering.
  Mat getPetscToStored() const { return mat_petsc_to_stored; }

protected:
  PetscLib lib;

  void buildPermutation(PetscInt nlocal, PetscInt nglobal,
                        const std::vector<int>& stored_indices);

  PetscInt global_size{0};
  PetscInt row_start{0};
  PetscInt row_end{0};
  std::vector<int> local_stored_indices;
  std::unordered_map<int, PetscInt> local_stored_to_petsc;
  Mat mat_stored_to_petsc{nullptr};
  Mat mat_petsc_to_stored{nullptr};
};

/// Cell-space mapping between stored cell numbers, PETSc numbering, and Field3D.
class PetscCellMapping : public PetscIndexMapping {
public:
  PetscCellMapping(const Field3D& cell_number, const Field3D& forward_cell_number,
                   const Field3D& backward_cell_number, int total_cells);

  static Region<Ind3D> create_region(const Field3D& cell_number);
  static Region<Ind3D> create_region_xin(const Field3D& cell_number);
  static Region<Ind3D> create_region_xout(const Field3D& cell_number);

  template <typename Function>
  void mapOwnedInteriorCells(Function func) const {
    PetscInt row = rowStart();
    BOUT_FOR_SERIAL(i, evolving_region) {
      func(row, i, ROUND(cell_number[i]));
      ++row;
    }
  }

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

  template <typename Function>
  void mapLocalYup(Function func) const {
    PetscInt row = evolving_region.size() + xin_region.size() + xout_region.size();
    BOUT_FOR_SERIAL(i, yup_region) {
      func(row, i.yp());
      ++row;
    }
  }

  template <typename Function>
  void mapLocalYdown(Function func) const {
    PetscInt row = evolving_region.size() + xin_region.size() + xout_region.size()
                   + yup_region.size();
    BOUT_FOR_SERIAL(i, ydown_region) {
      func(row, i.ym());
      ++row;
    }
  }

private:
  Field3D cell_number;
  Field3D forward_cell_number;
  Field3D backward_cell_number;
  Region<Ind3D> evolving_region;
  Region<Ind3D> xin_region;
  Region<Ind3D> xout_region;
  Region<Ind3D> yup_region;
  Region<Ind3D> ydown_region;
};

/// Mapping for forward or backward leg spaces.
class PetscLegMapping : public PetscIndexMapping {
public:
  PetscLegMapping(int total_legs, std::vector<int> local_leg_indices);
};

using PetscCellMappingPtr = std::shared_ptr<const PetscCellMapping>;
using PetscLegMappingPtr = std::shared_ptr<const PetscLegMapping>;

/// Wrapper around a PETSc Vec with a compile-time space tag.
template <typename SpaceTag>
class PetscSpaceVector {
  using UniqueVec = bout::petsc::UniqueVec;

public:
  PetscSpaceVector() = default;

  explicit PetscSpaceVector(std::shared_ptr<const PetscIndexMapping> mapping)
      : mapping(std::move(mapping)),
        vec(createVec(this->mapping->localSize(), this->mapping->globalSize())) {}

  PetscSpaceVector(std::shared_ptr<const PetscIndexMapping> mapping, UniqueVec vec)
      : mapping(std::move(mapping)), vec(std::move(vec)) {}

  PetscSpaceVector(const PetscSpaceVector&) = delete;
  PetscSpaceVector& operator=(const PetscSpaceVector&) = delete;
  PetscSpaceVector(PetscSpaceVector&&) noexcept = default;
  PetscSpaceVector& operator=(PetscSpaceVector&&) noexcept = default;

  Vec raw() const { return *this->vec; }
  const std::shared_ptr<const PetscIndexMapping>& getMapping() const { return mapping; }

  PetscSpaceVector duplicate() const {
    UniqueVec out{new Vec{nullptr}};
    BOUT_DO_PETSC(VecDuplicate(*this->vec, out.get()));
    return PetscSpaceVector(this->mapping, std::move(out));
  }

  PetscSpaceVector operator*(BoutReal scalar) const {
    auto out = this->duplicate();
    BOUT_DO_PETSC(VecScale(out.raw(), scalar));
    return out;
  }

  PetscSpaceVector operator+(const PetscSpaceVector& rhs) const {
    ASSERT0(mapping == rhs.mapping);
    auto out = this->duplicate();
    BOUT_DO_PETSC(VecAXPY(out.raw(), 1.0, rhs.raw()));
    return out;
  }

  PetscSpaceVector operator-(const PetscSpaceVector& rhs) const {
    ASSERT0(mapping == rhs.mapping);
    auto out = this->duplicate();
    BOUT_DO_PETSC(VecAXPY(out.raw(), -1.0, rhs.raw()));
    return out;
  }

  static PetscSpaceVector pointwiseMultiply(const PetscSpaceVector& lhs,
                                            const PetscSpaceVector& rhs) {
    ASSERT0(lhs.mapping == rhs.mapping);
    auto out = lhs.duplicate();
    BOUT_DO_PETSC(VecPointwiseMult(out.raw(), lhs.raw(), rhs.raw()));
    return out;
  }

  static PetscSpaceVector reciprocal(const PetscSpaceVector& in) {
    auto out = in.duplicate();
    BOUT_DO_PETSC(VecReciprocal(out.raw()));
    return out;
  }

private:
  static UniqueVec createVec(PetscInt local_size, PetscInt global_size) {
    UniqueVec out{new Vec{nullptr}};
    BOUT_DO_PETSC(VecCreate(BoutComm::get(), out.get()));
    BOUT_DO_PETSC(VecSetType(*out, VECMPI));
    BOUT_DO_PETSC(VecSetSizes(*out, local_size, global_size));
    return out;
  }

  std::shared_ptr<const PetscIndexMapping> mapping;
  UniqueVec vec;
};

// Tagged vectors provide compile-time check of compatible operations
using PetscCellVector = PetscSpaceVector<CellSpaceTag>;
using PetscForwardLegVector = PetscSpaceVector<ForwardLegSpaceTag>;
using PetscBackwardLegVector = PetscSpaceVector<BackwardLegSpaceTag>;

PetscCellVector makePetscCellVector(const PetscCellMappingPtr& mapping,
                                    const Field3D& field);
Field3D toField3D(const PetscCellVector& vec, const Field3D& prototype);

/// Typed linear operator between two spaces.
template <typename OutSpaceTag, typename InSpaceTag>
class PetscOperator {
  using UniqueMat = bout::petsc::UniqueMat;

public:
  PetscOperator() = default;
  ~PetscOperator() = default;

  PetscOperator(std::shared_ptr<const PetscIndexMapping> out_mapping,
                std::shared_ptr<const PetscIndexMapping> in_mapping, Array<int> rows,
                Array<int> cols, Array<BoutReal> weights)
      : out_mapping(std::move(out_mapping)), in_mapping(std::move(in_mapping)) {
    assembleFromCSR(rows, cols, weights);
  }

  PetscOperator(std::shared_ptr<const PetscIndexMapping> out_mapping,
                std::shared_ptr<const PetscIndexMapping> in_mapping, UniqueMat mat)
      : out_mapping(std::move(out_mapping)), in_mapping(std::move(in_mapping)),
        mat_operator(std::move(mat)) {}

  // Cannot be copied but can be moved
  PetscOperator(const PetscOperator&) = delete;
  PetscOperator& operator=(const PetscOperator&) = delete;
  PetscOperator(PetscOperator&&) noexcept = default;
  PetscOperator& operator=(PetscOperator&&) noexcept = default;

  PetscSpaceVector<OutSpaceTag>
  operator()(const PetscSpaceVector<InSpaceTag>& rhs) const {
    ASSERT0(in_mapping == rhs.getMapping());
    PetscSpaceVector<OutSpaceTag> result(out_mapping);
    BOUT_DO_PETSC(MatMult(*this->mat_operator, rhs.raw(), result.raw()));
    return result;
  }

  template <typename T = InSpaceTag,
            typename = std::enable_if_t<std::is_same_v<T, CellSpaceTag>>>
  PetscSpaceVector<OutSpaceTag> operator()(const Field3D& rhs) const {
    auto in = makePetscCellVector(
        std::static_pointer_cast<const PetscCellMapping>(in_mapping), rhs);
    return (*this)(in);
  }

  template <typename T = OutSpaceTag,
            typename = std::enable_if_t<std::is_same_v<T, CellSpaceTag>>>
  Field3D applyToField(const PetscSpaceVector<InSpaceTag>& rhs,
                       const Field3D& prototype) const {
    auto out = (*this)(rhs);
    return toField3D(out, prototype);
  }

  template <typename T = InSpaceTag, typename U = OutSpaceTag,
            typename = std::enable_if_t<std::is_same_v<T, CellSpaceTag>
                                        && std::is_same_v<U, CellSpaceTag>>>
  Field3D operator()(const Field3D& rhs, const Field3D& prototype) const {
    return applyToField((*this)(rhs), prototype);
  }

  PetscOperator<InSpaceTag, OutSpaceTag> transpose() const {
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatTranspose(*mat_operator, MAT_INITIAL_MATRIX, out.get()));
    return PetscOperator<InSpaceTag, OutSpaceTag>(in_mapping, out_mapping,
                                                  std::move(out));
  }

  PetscOperator operator+(const PetscOperator& rhs) const {
    ASSERT0(out_mapping == rhs.out_mapping);
    ASSERT0(in_mapping == rhs.in_mapping);
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatDuplicate(*this->mat_operator, MAT_COPY_VALUES, out.get()));
    BOUT_DO_PETSC(MatAXPY(*out, 1.0, *rhs.mat_operator, UNKNOWN_NONZERO_PATTERN));
    return PetscOperator(out_mapping, in_mapping, std::move(out));
  }

  PetscOperator operator-(const PetscOperator& rhs) const {
    ASSERT0(out_mapping == rhs.out_mapping);
    ASSERT0(in_mapping == rhs.in_mapping);
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatDuplicate(*this->mat_operator, MAT_COPY_VALUES, out.get()));
    BOUT_DO_PETSC(MatAXPY(*out, -1.0, *rhs.mat_operator, UNKNOWN_NONZERO_PATTERN));
    return PetscOperator(out_mapping, in_mapping, std::move(out));
  }

  PetscOperator operator*(BoutReal scalar) const& {
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatDuplicate(*this->mat_operator, MAT_COPY_VALUES, out.get()));
    BOUT_DO_PETSC(MatScale(*out, scalar));
    return PetscOperator(out_mapping, in_mapping, std::move(out));
  }

  PetscOperator operator*(BoutReal scalar) && {
    BOUT_DO_PETSC(MatScale(*this->mat_operator, scalar));
    return std::move(*this);
  }

  Mat raw() const { return *this->mat_operator; }

  static PetscOperator diagonal(std::shared_ptr<const PetscIndexMapping> mapping,
                                const PetscSpaceVector<OutSpaceTag>& diag) {
    static_assert(std::is_same_v<OutSpaceTag, InSpaceTag>,
                  "Diagonal operators require the same input and output space");
    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatCreate(BoutComm::get(), out.get()));
    BOUT_DO_PETSC(MatSetType(*out, MATMPIAIJ));
    BOUT_DO_PETSC(MatSetSizes(*out, mapping->localSize(), mapping->localSize(),
                              mapping->globalSize(), mapping->globalSize()));
    BOUT_DO_PETSC(MatMPIAIJSetPreallocation(*out, 1, nullptr, 0, nullptr));
    BOUT_DO_PETSC(MatDiagonalSet(*out, diag.raw(), INSERT_VALUES));
    BOUT_DO_PETSC(MatAssemblyBegin(*out, MAT_FINAL_ASSEMBLY));
    BOUT_DO_PETSC(MatAssemblyEnd(*out, MAT_FINAL_ASSEMBLY));
    return PetscOperator(mapping, mapping, std::move(out));
  }

private:
  void assembleFromCSR(const Array<int>& rows, const Array<int>& cols,
                       const Array<BoutReal>& weights) {
    UniqueMat temp{new Mat{nullptr}};
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

    UniqueMat out{new Mat{nullptr}};
    BOUT_DO_PETSC(MatMatMult(*temp, in_mapping->getPetscToStored(), MAT_INITIAL_MATRIX,
                             PETSC_DETERMINE, out.get()));
    mat_operator = std::move(out);
  }

  std::shared_ptr<const PetscIndexMapping> out_mapping;
  std::shared_ptr<const PetscIndexMapping> in_mapping;
  UniqueMat mat_operator{new Mat{nullptr}};

  template <typename O, typename M, typename I>
  friend PetscOperator<O, I> operator*(const PetscOperator<O, M>& lhs,
                                       const PetscOperator<M, I>& rhs);
};

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

using PetscForwardOperator = PetscOperator<ForwardLegSpaceTag, CellSpaceTag>;
using PetscBackwardOperator = PetscOperator<BackwardLegSpaceTag, CellSpaceTag>;
using PetscForwardToCellOperator = PetscOperator<CellSpaceTag, ForwardLegSpaceTag>;
using PetscBackwardToCellOperator = PetscOperator<CellSpaceTag, BackwardLegSpaceTag>;
using PetscCellOperator = PetscOperator<CellSpaceTag, CellSpaceTag>;
using PetscForwardLegOperator = PetscOperator<ForwardLegSpaceTag, ForwardLegSpaceTag>;
using PetscBackwardLegOperator = PetscOperator<BackwardLegSpaceTag, BackwardLegSpaceTag>;

/// Factory for operators defined on cell and leg spaces.
class PetscOperators {
public:
  explicit PetscOperators(Mesh* mesh);

  PetscForwardOperator forward() const;
  PetscBackwardOperator backward() const;

  PetscForwardOperator injectForward() const;
  PetscBackwardOperator injectBackward() const;

  PetscForwardLegVector forwardLegWeights() const;
  PetscBackwardLegVector backwardLegWeights() const;

  PetscCellOperator diagonal(const Field3D& f) const;
  PetscForwardLegOperator diagonal(const PetscForwardLegVector& v) const;
  PetscBackwardLegOperator diagonal(const PetscBackwardLegVector& v) const;

  struct Parallel {
    PetscForwardOperator Forward;
    PetscBackwardOperator Backward;
    PetscForwardOperator Inject_plus;
    PetscBackwardOperator Inject_minus;
    PetscForwardOperator Interp_plus;
    PetscBackwardOperator Interp_minus;
    PetscForwardOperator Grad_plus;
    PetscBackwardOperator Grad_minus;
    PetscForwardToCellOperator Div_minus;
    PetscBackwardToCellOperator Div_plus;
    PetscCellOperator Div_par_Grad_par;
    Field3D dV;

    Field3D Div_par_K_Grad_par(const Field3D& K, const Field3D& f) const;
  };

  Parallel getParallel() const;

private:
  static Field3D meshGetField3D(Mesh* mesh, const std::string& name);

  template <typename T>
  Array<T> meshGetArray(const std::string& name) const {
    Array<T> result;
    if (mesh->get(result, name) != 0) {
      throw BoutException("PetscOperators requires array '{}'", name);
    }
    return result;
  }

  int meshGetInt(const std::string& name) const;

  template <typename LegTag>
  PetscOperator<LegTag, CellSpaceTag>
  buildInjection(const Field3D& interior_leg_number, const Field3D& boundary_leg_number,
                 std::shared_ptr<const PetscLegMapping> leg_mapping) const;

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

  Mesh* mesh{nullptr};
  PetscCellMappingPtr cell_mapping;
  PetscLegMappingPtr forward_leg_mapping;
  PetscLegMappingPtr backward_leg_mapping;
  Field3D forward_leg_interior_number;
  Field3D forward_leg_boundary_number;
  Field3D backward_leg_interior_number;
  Field3D backward_leg_boundary_number;
};

#else

#warning PETSc not available. No PetscOperators.

#endif
