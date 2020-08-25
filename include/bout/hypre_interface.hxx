#ifndef BOUT_HYPRE_INTERFACE_H
#define BOUT_HYPRE_INTERFACE_H

#ifdef BOUT_HAS_HYPRE

#include "boutcomm.hxx"
#include "field.hxx"
#include "utils.hxx"
#include "bout/globalindexer.hxx"

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_utilities.h"

#include <memory>

namespace bout {

template <class T>
class HypreVector;
template <class T>
void swap(HypreVector<T>& lhs, HypreVector<T>& rhs);

// TODO: check correctness of initialise/assemble
// TODO: set sizes
// TODO: set contiguous blocks at once

template <class T>
class HypreVector {
  HYPRE_IJVector hypre_vector{nullptr};
  HYPRE_ParVector parallel_vector{nullptr};
  IndexerPtr<T> indexConverter;
  CELL_LOC location;
  bool initialised{false};

public:
  static_assert(bout::utils::is_Field<T>::value, "HypreVector only works with Fields");
  using ind_type = typename T::ind_type;

  HypreVector() = default;
  ~HypreVector() {
    // Destroy may print an error if the vector is already null
    if (hypre_vector == nullptr) {
      return;
    }
    // Also destroys parallel_vector
    HYPRE_IJVectorDestroy(hypre_vector);
  }

  // Disable copy, at least for now: not clear that HYPRE_IJVector is
  // copy-able
  HypreVector(const HypreVector<T>&) = delete;
  auto operator=(const HypreVector<T>&) = delete;

  HypreVector(HypreVector<T>&& other) {
    std::swap(hypre_vector, other.hypre_vector);
    std::swap(parallel_vector, other.parallel_vector);
    indexConverter = other.indexConverter;
    location = other.location;
    initialised = other.initialised;
    other.initialised = false;
  }

  HypreVector<T>& operator=(HypreVector<T>&& other) {
    std::swap(hypre_vector, other.hypre_vector);
    std::swap(parallel_vector, other.parallel_vector);
    indexConverter = other.indexConverter;
    location = other.location;
    initialised = other.initialised;
    other.initialised = false;
    return *this;
  }
  
  /// Construct from a field, copying over the field values
  explicit HypreVector(const T& f, IndexerPtr<T> indConverter)
      : indexConverter(indConverter) {
    ASSERT1(indConverter->getMesh() == f.getMesh());
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? f.getMesh()->getXcomm() : BoutComm::get();

    const auto& region = f.getRegion("RGN_ALL_THIN");
    const auto jlower =
        static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*std::begin(region)));
    const auto jupper =
        static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*--std::end(region)));

    HYPRE_IJVectorCreate(comm, jlower, jupper, &hypre_vector);
    HYPRE_IJVectorSetObjectType(hypre_vector, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(hypre_vector);

    location = f.getLocation();
    initialised = true;

    BOUT_FOR_SERIAL(i, region) {
      const auto index = static_cast<HYPRE_BigInt>(indexConverter->getGlobal(i));
      if (index != -1) {
        HYPRE_IJVectorSetValues(hypre_vector, static_cast<HYPRE_Int>(1), &index, &f[i]);
      }
    }

    // FIXME: call HYPRE_IJVectorSetMaxOffProcElmts

    assemble();
  }

  HypreVector<T>& operator=(const T& f) {
    HypreVector<T> result(f, indexConverter);
    *this = std::move(result);
    return *this;
  }

  // Mesh can't be const yet, need const-correctness on Mesh;
  // GlobalIndexer ctor also modifies mesh -- FIXME
  HypreVector(Mesh& mesh, IndexerPtr<T> indConverter)
      : indexConverter(indConverter) {
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? mesh.getXcomm() : BoutComm::get();
    ASSERT1(indConverter->getMesh() == &mesh);

    const auto& region = mesh.getRegion<T>("RGN_ALL_THIN");

    const auto jlower =
        static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*std::begin(region)));
    const auto jupper =
        static_cast<HYPRE_BigInt>(indexConverter->getGlobal(*--std::end(region)));

    HYPRE_IJVectorCreate(comm, jlower, jupper, &hypre_vector);
    HYPRE_IJVectorSetObjectType(hypre_vector, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(hypre_vector);
    initialised = true;
    location = CELL_LOC::centre;
  }

  void assemble() {
    HYPRE_IJVectorAssemble(hypre_vector);
    HYPRE_IJVectorGetObject(hypre_vector, reinterpret_cast<void**>(&parallel_vector));
  }

  T toField() {
    T result(indexConverter->getMesh());
    result.allocate().setLocation(location);

    // Note that this only populates boundaries to a depth of 1
    BOUT_FOR_SERIAL(i, result.getRegion("RGN_ALL_THIN")) {
      const auto index = static_cast<HYPRE_BigInt>(indexConverter->getGlobal(i));
      if (index != -1) {
        // Yes, complex, but this is a HYPRE typedef for real
        HYPRE_Complex value;
        HYPRE_IJVectorGetValues(hypre_vector, 1, &index, &value);
        result[i] = static_cast<BoutReal>(value);
      }
    }

    return result;
  }

  HYPRE_IJVector get() { return hypre_vector; }
  const HYPRE_IJVector& get() const { return hypre_vector; }

  HYPRE_ParVector getParallel() { return parallel_vector; }
  const HYPRE_ParVector& getParallel() const { return parallel_vector; }

  class Element {
    HYPRE_IJVector vector{nullptr};
    HYPRE_BigInt index;
    HYPRE_Complex value;

  public:
    Element() = delete;
    Element(const Element&) = default;
    Element(Element&&) = default;
    Element(HYPRE_IJVector vector_, int index_) : vector(vector_), index(index_) {
      BOUT_OMP(critical) {
        if (HYPRE_IJVectorGetValues(vector, 1, &index, &value) != 0) {
          throw BoutException("Error when getting element of a HYPRE vector, index = {}",
                              index);
        }
      }
    }

    Element& operator=(const Element& other) {
      return *this = static_cast<BoutReal>(other);
    }
    Element& operator=(BoutReal value_) {
      value = value_;
      BOUT_OMP(critical) {
        if (HYPRE_IJVectorSetValues(vector, 1, &index, &value) != 0) {
          throw BoutException("Error when setting element of a HYPRE vector, index = {}",
                              index);
        }
      }
      return *this;
    }
    Element& operator+=(BoutReal value_) {
      value = value_;
      BOUT_OMP(critical) {
        if (HYPRE_IJVectorAddToValues(vector, 1, &index, &value) != 0) {
          throw BoutException(
              "Error when adding to element of a HYPRE vector, index = {}", index);
        }
      }
      return *this;
    }
    operator BoutReal() const { return value; }
  };

  Element operator()(ind_type& index) {
#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Can not return element of uninitialised vector");
    }
#endif
    const int global = indexConverter->getGlobal(index);
#if CHECKLEVEL >= 1
    if (global == -1) {
      throw BoutException("Request to return invalid vector element");
    }
#endif
    return Element(hypre_vector, global);
  }

  friend void swap(HypreVector<T>& lhs, HypreVector<T>& rhs) {
    using std::swap;
    swap(lhs.hypre_vector, rhs.hypre_vector);
    swap(lhs.parallel_vector, rhs.parallel_vector);
    swap(lhs.indexConverter, rhs.indexConverter);
    swap(lhs.location, rhs.location);
    swap(lhs.initialised, rhs.initialised);
  }
};

template <class T>
class HypreMatrix;
template <class T>
void swap(HypreMatrix<T>& lhs, HypreMatrix<T>& rhs);

template <class T>
class HypreMatrix {
  std::shared_ptr<HYPRE_IJMatrix> hypre_matrix{nullptr};
  HYPRE_ParCSRMatrix parallel_matrix{nullptr};
  IndexerPtr<T> index_converter;
  bool initialised{false};
  int yoffset{0};
  ParallelTransform* parallel_transform{nullptr};
  bool assembled{false};

  struct MatrixDeleter {
    void operator()(HYPRE_IJMatrix* matrix) const {
      if (*matrix == nullptr) {
        return;
      }
      HYPRE_IJMatrixDestroy(*matrix);
      delete matrix;
    }
  };

public:
  static_assert(bout::utils::is_Field<T>::value, "HypreMatrix only works with Fields");
  using ind_type = typename T::ind_type;

  HypreMatrix() = default;
  HypreMatrix(const HypreMatrix<T>&) = delete;
  HypreMatrix(HypreMatrix<T>&& other) {
    using std::swap;
    swap(hypre_matrix, other.hypre_matrix);
    swap(parallel_matrix, other.parallel_matrix);
    index_converter = other.index_converter;
    initialised = other.initialised;
    yoffset = other.yoffset;
    parallel_transform = other.parallel_transform;
    other.initialised = false;
  }
  HypreMatrix<T>& operator=(const HypreMatrix<T>&) = delete;
  HypreMatrix<T>& operator=(HypreMatrix<T>&& other) {
    using std::swap;
    swap(hypre_matrix, other.hypre_matrix);
    swap(parallel_matrix, other.parallel_matrix);
    index_converter = other.index_converter;
    initialised = other.initialised;
    yoffset = other.yoffset;
    parallel_transform = other.parallel_transform;
    other.initialised = false;
    return *this;
  }
  
  /// Construct a matrix capable of operating on the specified field,
  /// preallocating memory if requeted and possible.
  ///
  /// note: preallocate not currently used, but here to match PetscMatrix interface 
  explicit HypreMatrix(IndexerPtr<T> indConverter, bool preallocate = true)
      : hypre_matrix(new HYPRE_IJMatrix, MatrixDeleter{}), index_converter(indConverter) {
    Mesh *mesh = indConverter->getMesh();
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? mesh->getXcomm() : BoutComm::get();
    parallel_transform = &mesh->getCoordinates()->getParallelTransform();

    const auto& region = mesh->getRegion<T>("RGN_ALL_THIN");

    const auto ilower =
        static_cast<HYPRE_BigInt>(index_converter->getGlobal(*std::begin(region)));
    const auto iupper =
        static_cast<HYPRE_BigInt>(index_converter->getGlobal(*--std::end(region)));

    HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &*hypre_matrix);
    HYPRE_IJMatrixSetObjectType(*hypre_matrix, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(*hypre_matrix);
    // TODO: redirect printf into BoutException using freopen
    // Second argument here is meaningless
    // HYPRE_IJMatrixSetPrintLevel(hypre_matrix, 0);

    initialised = true;
  }

  class Element {
    HYPRE_IJMatrix hypre_matrix;
    HYPRE_BigInt row, column;
    HYPRE_Complex value;
    std::vector<HYPRE_BigInt> positions;
    std::vector<BoutReal> weights;

  public:
    Element() = delete;
    Element(const Element&) = default;
    Element(Element&&) = default;

    Element(HypreMatrix<T>& matrix_, HYPRE_BigInt row_, HYPRE_BigInt column_,
            std::vector<HYPRE_BigInt> positions_ = {},
            std::vector<BoutReal> weights_ = {})
        : hypre_matrix(matrix_.get()), row(row_), column(column_), positions(positions_),
          weights(weights_) {
      ASSERT2(positions.size() == weights.size());
      if (positions.empty()) {
        positions = {column};
        weights = {1.0};
      }
      if (matrix_.isAssembled()) {
        HYPRE_Int ncolumns{1};
        BOUT_OMP(critial)
        if (HYPRE_IJMatrixGetValues(hypre_matrix, 1, &ncolumns, &row, &column, &value)
            != 0) {
          // Not clearing the (global) error will break future calls!
          HYPRE_ClearAllErrors();
          throw BoutException("Error when getting value of matrix, row = {}, column = {}",
                              row, column);
        }
      } else {
        value = 0.0;
      }
    }
    Element& operator=(const Element& other) {
      return *this = static_cast<BoutReal>(other);
    }
    Element& operator=(BoutReal value_) {
      value = value_;
      setValues(value);
      return *this;
    }
    Element& operator+=(BoutReal value_) {
      auto column_position = std::find(cbegin(positions), cend(positions), column);
      if (column_position != cend(positions)) {
        const auto i = std::distance(cbegin(positions), column_position);
        value += weights[i] * value_;
      }
      setValues(value);
      return *this;
    }

    operator BoutReal() const { return value; }

  private:
    void setValues(BoutReal value_) {
      ASSERT3(!positions.empty());
      std::vector<HYPRE_Complex> values;
      std::transform(
          weights.begin(), weights.end(), std::back_inserter(values),
          [&value_](BoutReal weight) -> HYPRE_Complex { return weight * value_; });
      auto ncolumns = static_cast<HYPRE_Int>(positions.size());
      BOUT_OMP(critical)
      if (HYPRE_IJMatrixSetValues(hypre_matrix, 1, &ncolumns, &row, positions.data(),
                                  values.data())
          != 0) {
        HYPRE_ClearAllErrors();
        throw BoutException("Error when setting elements of a HYPRE matrix, row = {}, "
                            "column = {}",
                            row, column);
      }
    }
  };

  BoutReal operator()(const ind_type& row, const ind_type& column) const {
    const auto global_row = index_converter->getGlobal(row);
    const auto global_column = index_converter->getGlobal(column);
#if CHECKLEVEL >= 1
    if (global_row == -1 or global_column == -1) {
      throw BoutException(
          "Invalid HypreMatrix element at row = ({}, {}, {}), column = ({}, {}, {})",
          row.x(), row.y(), row.z(), column.x(), column.y(), column.z());
    }
#endif
    return this->operator()(global_row, global_column);
  }

  BoutReal operator()(const HYPRE_BigInt& row, const HYPRE_BigInt& column) const {
#if CHECKLEVEL >= 1
    if (not initialised) {
      throw BoutException("Cannot return element of uninitialised HypreMatrix");
    }
    if (not assembled) {
      throw BoutException("Cannot return element of unassembled HypreMatrix (call "
                          "HypreMatrix::assemble() first)");
    }
#endif
    HYPRE_Complex value;
    // The following need to be non-const for some reason
    HYPRE_Int ncolumns{1};
    HYPRE_BigInt row_{row};
    HYPRE_BigInt column_{column};
    BOUT_OMP(critial)
    if (HYPRE_IJMatrixGetValues(*hypre_matrix, 1, &ncolumns, &row_, &column_, &value)
        != 0) {
      // Not clearing the (global) error will break future calls!
      HYPRE_ClearAllErrors();
      throw BoutException("Error when getting value of matrix, row = {}, column = {}",
                          row, column);
    }
    return static_cast<BoutReal>(value);
  }

  Element operator()(const ind_type& row, const ind_type& column) {
    const auto global_row = index_converter->getGlobal(row);
    const auto global_column = index_converter->getGlobal(column);

#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Cannot return element of uninitialised HypreMatrix");
    }
    if (global_row == -1 or global_column == -1) {
      throw BoutException(
          "Invalid HypreMatrix element at row = ({}, {}, {}), column = ({}, {}, {})",
          row.x(), row.y(), row.z(), column.x(), column.y(), column.z());
    }
#endif

    std::vector<HYPRE_Int> positions{};
    std::vector<HYPRE_Complex> weights{};

    if (yoffset != 0) {
      ASSERT1(yoffset == (column.y() - row.y()));
      const auto pw = [this, &row, &column]() {
        if (this->yoffset == -1) {
          return parallel_transform->getWeightsForYDownApproximation(column.x(), row.y(),
                                                                     column.z());
        } else if (this->yoffset == 1) {
          return parallel_transform->getWeightsForYUpApproximation(column.x(), row.y(),
                                                                   column.z());
        } else {
          return parallel_transform->getWeightsForYApproximation(
              column.x(), row.y(), column.z(), this->yoffset);
        }
      }();

      const int ny =
          std::is_same<T, FieldPerp>::value ? 1 : index_converter->getMesh()->LocalNy;
      const int nz =
          std::is_same<T, Field2D>::value ? 1 : index_converter->getMesh()->LocalNz;
      std::transform(
          pw.begin(), pw.end(), std::back_inserter(positions),
          [this, ny, nz](ParallelTransform::PositionsAndWeights p) -> HYPRE_Int {
            return this->index_converter->getGlobal(
                ind_type(p.i * ny * nz + p.j * nz + p.k, ny, nz));
          });
      std::transform(pw.begin(), pw.end(), std::back_inserter(weights),
                     [](ParallelTransform::PositionsAndWeights p) -> HYPRE_Complex {
                       return p.weight;
                     });
    }
    return Element(*this, global_row, global_column, positions, weights);
  }

  void assemble() {
    HYPRE_IJMatrixAssemble(*hypre_matrix);
    HYPRE_IJMatrixGetObject(*hypre_matrix, reinterpret_cast<void**>(&parallel_matrix));
    assembled = true;
  }

  bool isAssembled() const { return assembled; }

  HypreMatrix<T> yup(int index = 0) { return ynext(index + 1); }
  HypreMatrix<T> ydown(int index = 0) { return ynext(-index - 1); }
  HypreMatrix<T> ynext(int dir) {
    if (std::is_same<T, FieldPerp>::value and ((yoffset + dir) != 0)) {
      throw BoutException("Can not get ynext for FieldPerp");
    }
    HypreMatrix<T> result;
    result.hypre_matrix = hypre_matrix;
    result.index_converter = index_converter;
    result.parallel_transform = parallel_transform;
    result.yoffset = std::is_same<T, Field2D>::value ? 0 : yoffset + dir;
    result.initialised = initialised;
    return result;
  }

  HYPRE_IJMatrix get() { return *hypre_matrix; }
  const HYPRE_IJMatrix& get() const { return *hypre_matrix; }

  HYPRE_ParCSRMatrix getParallel() { return parallel_matrix; }
  const HYPRE_ParCSRMatrix& getParallel() const { return parallel_matrix; }
};

} // namespace bout

#endif // BOUT_HAS_HYPRE

#endif // BOUT_HYPRE_INTERFACE_H
