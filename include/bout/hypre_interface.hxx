#ifndef BOUT_HYPRE_INTERFACE_H
#define BOUT_HYPRE_INTERFACE_H

#include "bout/build_config.hxx"

#ifdef BOUT_HAS_HYPRE

#include "bout/bout_enum_class.hxx"
#include "boutcomm.hxx"
#include "field.hxx"
#include "utils.hxx"
#include "bout/globalindexer.hxx"

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_utilities.h"
#include "_hypre_utilities.h"


#include <memory>

// BOUT_ENUM_CLASS does not work inside namespaces
BOUT_ENUM_CLASS(HYPRE_SOLVER_TYPE, gmres);

namespace bout {
#ifdef BOUT_USE_CUDA   //HYPRE with Cuda enabled 
#define HypreMalloc(P, SIZE)  cudaMallocManaged(&P, SIZE)
#define HypreFree(P)  cudaFree(P)
#else
#define HypreMalloc(P, SIZE) P = static_cast<decltype(P)>(malloc(SIZE))
#define HypreFree(P)  free(P)
#endif


namespace {
int checkHypreError(int error) {
  if (error) {
    // Aaron Fisher warns that it is possible that the error code is non-zero even
    // though a solve was successful. If this occurs, we might want to make this a
    // warning, or otherwise investigate further.
    char error_cstr[2048];
    HYPRE_DescribeError(error, &error_cstr[0]);
    if(error > 1)
      throw BoutException("A Hypre call failed with error code {:d}: {:s}",
                        error, &error_cstr[0]);
    else if(error == 1) {
       printf("Hypre returns %d %s\n",error,&error_cstr[0]);
       HYPRE_ClearAllErrors();
    }  

  }
  return error;
}
}


template <class T>
class HypreVector;
template <class T>
void swap(HypreVector<T>& lhs, HypreVector<T>& rhs);

// TODO: check correctness of initialise/assemble
// TODO: set sizes
// TODO: set contiguous blocks at once

template <class T>
class HypreVector {
  MPI_Comm comm;
  HYPRE_BigInt jlower, jupper, vsize;
  HYPRE_IJVector hypre_vector{nullptr};
  HYPRE_ParVector parallel_vector{nullptr};
  IndexerPtr<T> indexConverter;
  CELL_LOC location;
  bool initialised{false};
  bool have_indices{false};
  HYPRE_BigInt *I{nullptr};
  HYPRE_Complex *V{nullptr};


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
    checkHypreError(HYPRE_IJVectorDestroy(hypre_vector));
    HypreFree(I);  
    HypreFree(V);
  }

  // Disable copy, at least for now: not clear that HYPRE_IJVector is
  // copy-able
  HypreVector(const HypreVector<T>&) = delete;
  auto operator=(const HypreVector<T>&) = delete;

  HypreVector(HypreVector<T>&& other) {
    comm = other.comm;
    jlower = other.jlower;
    jupper = other.jupper;
    vsize = other.vsize;
    std::swap(hypre_vector, other.hypre_vector);
    std::swap(parallel_vector, other.parallel_vector);
    indexConverter = other.indexConverter;
    location = other.location;
    initialised = other.initialised;
    other.initialised = false;
    have_indices = other.have_indices;
    other.have_indices = false;
    I = other.I;
    V = other.V;
    other.I = nullptr;
    other.V = nullptr;
  }

  HypreVector<T>& operator=(HypreVector<T>&& other) {
    comm = other.comm;
    jlower = other.jlower;
    jupper = other.jupper;
    vsize = other.vsize;
    std::swap(hypre_vector, other.hypre_vector);
    std::swap(parallel_vector, other.parallel_vector);
    indexConverter = other.indexConverter;
    location = other.location;
    initialised = other.initialised;
    other.initialised = false;
    have_indices = other.have_indices;
    other.have_indices = false;
    I = other.I;
    V = other.V;
    other.I = nullptr;
    other.V = nullptr;
    return *this;
  }
  
  HypreVector<T>& operator=(const T& f) {
    ASSERT0(indexConverter); // Needs to have an index set
    HypreVector<T> result(f, indexConverter);
    *this = std::move(result);
    return *this;
  }

/// Construct from a field, copying over the field values
  explicit HypreVector(const T& f, IndexerPtr<T> indConverter)
      : indexConverter(indConverter) {
    ASSERT1(indConverter->getMesh() == f.getMesh());
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? f.getMesh()->getXcomm() : BoutComm::get();

    HYPRE_BigInt jlower = indConverter->getGlobalStart();
    HYPRE_BigInt jupper = jlower + indConverter->size() - 1; // inclusive end
    vsize = jupper - jlower + 1;

#ifdef BOUT_USE_CUDA
    hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_DEVICE;
    HYPRE_MemoryLocation memory_location = HYPRE_MEMORY_DEVICE;
    hypre_HandleMemoryLocation(hypre_handle())    = memory_location;
#else
    hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_HOST;
    HYPRE_MemoryLocation memory_location = HYPRE_MEMORY_HOST;
    hypre_HandleMemoryLocation(hypre_handle())    = memory_location;
#endif

    checkHypreError(HYPRE_IJVectorCreate(comm, jlower, jupper, &hypre_vector));
    checkHypreError(HYPRE_IJVectorSetObjectType(hypre_vector, HYPRE_PARCSR));
    checkHypreError(HYPRE_IJVectorInitialize(hypre_vector));

    location = f.getLocation();
    initialised = true;
    HypreMalloc(I, vsize*sizeof(HYPRE_BigInt)); 
    HypreMalloc(V, vsize*sizeof(HYPRE_Complex));
    importValuesFromField(f);
  }

  /// Construct a vector with given index set, but don't set any values
  HypreVector(IndexerPtr<T> indConverter) : indexConverter(indConverter) {
    Mesh& mesh = *indConverter->getMesh();
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? mesh.getXcomm() : BoutComm::get();
    
    HYPRE_BigInt jlower = indConverter->getGlobalStart();
    HYPRE_BigInt jupper = jlower + indConverter->size() - 1; // inclusive end
    vsize = jupper - jlower + 1;
#ifdef BOUT_USE_CUDA
    hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_DEVICE;
    HYPRE_MemoryLocation memory_location = HYPRE_MEMORY_DEVICE;
    hypre_HandleMemoryLocation(hypre_handle())    = memory_location;
#else
    hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_HOST;
    HYPRE_MemoryLocation memory_location = HYPRE_MEMORY_HOST;
    hypre_HandleMemoryLocation(hypre_handle())    = memory_location;
#endif

    checkHypreError(HYPRE_IJVectorCreate(comm, jlower, jupper, &hypre_vector));
    checkHypreError(HYPRE_IJVectorSetObjectType(hypre_vector, HYPRE_PARCSR));
    checkHypreError(HYPRE_IJVectorInitialize(hypre_vector));
    initialised = true;
    location = CELL_LOC::centre;
    HypreMalloc(I, vsize*sizeof(HYPRE_BigInt)); 
    HypreMalloc(V, vsize*sizeof(HYPRE_Complex));
  }

  void assemble() {
    writeCacheToHypre();
    checkHypreError(HYPRE_IJVectorAssemble(hypre_vector));
    checkHypreError(HYPRE_IJVectorGetObject(hypre_vector, reinterpret_cast<void**>(&parallel_vector)));
  }

  void writeCacheToHypre()
  {
    checkHypreError(HYPRE_IJVectorSetValues(hypre_vector, vsize, I, V));
  }

  void readCacheFromHypre()
  {
    checkHypreError(HYPRE_IJVectorGetValues(hypre_vector, vsize, I, V));
  }

  T toField() {
    T result(indexConverter->getMesh());
    result.allocate().setLocation(location);

    readCacheFromHypre();
    // Note that this only populates boundaries to a depth of 1
    int count = 0;
    BOUT_FOR_SERIAL(i, indexConverter->getRegionAll()) {
      HYPRE_BigInt index = static_cast<HYPRE_BigInt>(indexConverter->getGlobal(i));
      if (index != -1) { // Todo double check why index out of bounds does not return -1
        result[i] = static_cast<BoutReal>(V[count]);
        count++;
      }
    }

    return result;
  }

  void importValuesFromField(const T &f) {
    int vec_i = 0;
    BOUT_FOR_SERIAL(i, indexConverter->getRegionAll()) {
      HYPRE_BigInt index = static_cast<HYPRE_BigInt>(indexConverter->getGlobal(i));
      if (index != -1) {
        I[vec_i] = index;        
        V[vec_i] = f[i];
        vec_i ++;
      }
    }

    ASSERT2(vec_i == vsize);
    //writeCacheToHypre(); // redundant assemble already performs writeCacheToHypre
    assemble();
    have_indices = true;
  }

  HYPRE_IJVector get() { return hypre_vector; }
  const HYPRE_IJVector& get() const { return hypre_vector; }

  HYPRE_ParVector getParallel() { return parallel_vector; }
  const HYPRE_ParVector& getParallel() const { return parallel_vector; }

  class Element {
    HypreVector<T> *vector;
    HYPRE_IJVector hypre_vector{nullptr};
    HYPRE_BigInt index;
    int vec_i;
    HYPRE_Complex value;

  public:
    Element() = delete;
    Element(const Element&) = default;
    Element(Element&&) = default;
    Element(HypreVector<T> &vector_, int index_) : vector(&vector_), hypre_vector(vector_.get()), index(index_) {
      vector = &vector_;
      vec_i = -1;
      for (HYPRE_BigInt i = 0; i < vector->vsize; ++i) {
        if (vector->I[i] == index) {
          vec_i = i;
          break;
        }
      }
      if (vec_i < 0)
        throw BoutException("Error cannot find global index in HypreVector, index = {}",index);

      value = vector->V[vec_i];
    }

    Element& operator=(const Element& other) {
      ASSERT3(finite(static_cast<BoutReal>(other)));
      return *this = static_cast<BoutReal>(other);
    }
    Element& operator=(BoutReal value_) {
      ASSERT3(finite(value_));
      value = value_;
      vector->V[vec_i] = value_;
      return *this;
    }
    Element& operator+=(BoutReal value_) {
      ASSERT3(finite(value_));
      value += value_;
      ASSERT3(finite(value));
      vector->V[vec_i] += value_;
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
    return Element(*this, global);
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
  MPI_Comm comm;
  HYPRE_BigInt ilower, iupper;
  std::shared_ptr<HYPRE_IJMatrix> hypre_matrix{nullptr};
  HYPRE_ParCSRMatrix parallel_matrix{nullptr};
  IndexerPtr<T> index_converter;
  CELL_LOC location;
  bool initialised{false};
  int yoffset{0};
  ParallelTransform* parallel_transform{nullptr};
  bool assembled{false};
  HYPRE_BigInt num_rows;
  std::vector<HYPRE_BigInt> *I;
  std::vector<std::vector<HYPRE_BigInt> > *J;
  std::vector<std::vector<HYPRE_Complex> > *V;

  // todo also take care of I,J,V
  struct MatrixDeleter {
    void operator()(HYPRE_IJMatrix* matrix) const {
      if (*matrix == nullptr) {
        return;
      }
      checkHypreError(HYPRE_IJMatrixDestroy(*matrix));
    }
  };

public:
  static_assert(bout::utils::is_Field<T>::value, "HypreMatrix only works with Fields");
  using ind_type = typename T::ind_type;

  HypreMatrix() = default;
  HypreMatrix(const HypreMatrix<T>&) = delete;
  HypreMatrix(HypreMatrix<T>&& other) {
    comm = other.comm;
    ilower = other.ilower;
    iupper = other.iupper;
    std::swap(hypre_matrix, other.hypre_matrix);
    std::swap(parallel_matrix, other.parallel_matrix);
    index_converter = other.index_converter;
    location = other.location;
    initialised = other.initialised;
    yoffset = other.yoffset;
    parallel_transform = other.parallel_transform;
    assembled = other.assembled;
    num_rows = other.num_rows;
    I = other.I;
    J = other.J;
    V = other.V;
  }
  HypreMatrix<T>& operator=(const HypreMatrix<T>&) = delete;
  HypreMatrix<T>& operator=(HypreMatrix<T>&& other) {
    comm = other.comm;
    ilower = other.ilower;
    iupper = other.iupper;
    std::swap(hypre_matrix, other.hypre_matrix);
    std::swap(parallel_matrix, other.parallel_matrix);
    index_converter = other.index_converter;
    location = other.location;
    initialised = other.initialised;
    yoffset = other.yoffset;
    parallel_transform = other.parallel_transform;
    assembled = other.assembled;
    num_rows = other.num_rows;
    I = other.I;
    J = other.J;
    V = other.V;
    return *this;
  }
  
  /// Construct a matrix capable of operating on the specified field,
  /// preallocating memory if requeted and possible.
  ///
  /// note: preallocate not currently used, but here to match PetscMatrix interface 
  explicit HypreMatrix(IndexerPtr<T> indConverter, bool UNUSED(preallocate) = true)
      : hypre_matrix(new HYPRE_IJMatrix, MatrixDeleter{}), index_converter(indConverter) {
    Mesh *mesh = indConverter->getMesh();
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? mesh->getXcomm() : BoutComm::get();
    parallel_transform = &mesh->getCoordinates()->getParallelTransform();

    ilower = indConverter->getGlobalStart();
    iupper = ilower + indConverter->size() - 1; // inclusive end
    num_rows = iupper - ilower + 1;

    I = new std::vector<HYPRE_BigInt>;
    J = new std::vector<std::vector<HYPRE_BigInt> >;
    V = new std::vector<std::vector<HYPRE_Complex> >;
    (*I).resize(num_rows);
    (*J).resize(num_rows);
    (*V).resize(num_rows);
    for (HYPRE_BigInt i = 0; i < num_rows; ++i) {
      (*I)[i] = ilower + i;
      (*J)[i].resize(0);
      (*J)[i].reserve(10);
      (*V)[i].resize(0);
      (*V)[i].reserve(10);
    }

    checkHypreError(HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &*hypre_matrix));
    checkHypreError(HYPRE_IJMatrixSetObjectType(*hypre_matrix, HYPRE_PARCSR));
    checkHypreError(HYPRE_IJMatrixInitialize(*hypre_matrix));
    // TODO: redirect printf into BoutException using freopen
    // Second argument here is meaningless
    // HYPRE_IJMatrixSetPrintLevel(hypre_matrix, 0);

    initialised = true;
  }

  class Element {
    HypreMatrix *matrix;
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
#if CHECK > 2
      for (const auto val : weights) {
        ASSERT3(finite(val));
      }
#endif
      ASSERT2(positions.size() == weights.size());
      if (positions.empty()) {
        positions = {column};
        weights = {1.0};
      }
      matrix = &matrix_;
      value = matrix->getVal(row, column);
    }
    Element& operator=(const Element& other) {
      AUTO_TRACE();
      ASSERT3(finite(static_cast<BoutReal>(other)));
      return *this = static_cast<BoutReal>(other);
    }
    Element& operator=(BoutReal value_) {
      AUTO_TRACE();
      ASSERT3(finite(value_));
      value = value_;
      setValues(value);
      return *this;
    }
    Element& operator+=(BoutReal value_) {
      AUTO_TRACE();
      ASSERT3(finite(value_));
      auto column_position = std::find(cbegin(positions), cend(positions), column);
      if (column_position != cend(positions)) {
        const auto i = std::distance(cbegin(positions), column_position);
        value += weights[i] * value_;
      }
      addValues(value_);
      return *this;
    }

    operator BoutReal() const { return value; }

  private:
    void setValues(BoutReal value_) {
      TRACE("HypreMatrix setting values at ({}, {})", row, column);
      ASSERT3(!positions.empty());
      std::vector<HYPRE_Complex> values;
      std::transform(
          weights.begin(), weights.end(), std::back_inserter(values),
          [&value_](BoutReal weight) -> HYPRE_Complex { return weight * value_; });
      HYPRE_BigInt ncolumns = static_cast<HYPRE_BigInt>(positions.size());
      //BOUT_OMP(critical)
      for (HYPRE_BigInt i = 0; i < positions.size(); ++i) {
        matrix->setVal(row, positions[i], values[i]);
      }

    }

    void addValues(BoutReal value_) {
      TRACE("HypreMatrix setting values at ({}, {})", row, column);
      ASSERT3(!positions.empty());
      std::vector<HYPRE_Complex> values;
      std::transform(
          weights.begin(), weights.end(), std::back_inserter(values),
          [&value_](BoutReal weight) -> HYPRE_Complex { return weight * value_; });
      HYPRE_BigInt ncolumns = static_cast<HYPRE_BigInt>(positions.size());
      //BOUT_OMP(critical)
      for (HYPRE_BigInt i = 0; i < positions.size(); ++i) {
        matrix->addVal(row, positions[i], values[i]);
      }

    }
  };

  BoutReal getVal(const ind_type& row, const ind_type& column) const {
    const HYPRE_BigInt global_row = index_converter->getGlobal(row);
    const HYPRE_BigInt global_column = index_converter->getGlobal(column);
#if CHECKLEVEL >= 1
    if (global_row < 0 or global_column < 0) {
      throw BoutException(
          "Invalid HypreMatrix element at row = ({}, {}, {}), column = ({}, {}, {})",
          row.x(), row.y(), row.z(), column.x(), column.y(), column.z());
    }
#endif
    return getVal(global_row, global_column);
  }

  BoutReal getVal(const HYPRE_BigInt row, const HYPRE_BigInt column) const {
    HYPRE_Complex value = 0.0;
    HYPRE_BigInt i = row - ilower;
    ASSERT2(i >= 0 && i < num_rows);
    for(HYPRE_BigInt col_ind = 0; col_ind < (*J)[i].size(); ++col_ind) {
      if ((*J)[i][col_ind] == column) {
        value = (*V)[i][col_ind];
        break;
      }
    }
    return static_cast<BoutReal>(value);
  }

  void setVal(const ind_type& row, const ind_type& column, BoutReal value) {
    const HYPRE_BigInt global_row = index_converter->getGlobal(row);
    const HYPRE_BigInt global_column = index_converter->getGlobal(column);
#if CHECKLEVEL >= 1
    if (global_row < 0 or global_column < 0) {
      throw BoutException(
          "Invalid HypreMatrix element at row = ({}, {}, {}), column = ({}, {}, {})",
          row.x(), row.y(), row.z(), column.x(), column.y(), column.z());
    }
#endif
    return setVal(global_row, global_column, value);
  }

  void setVal(const HYPRE_BigInt row, const HYPRE_BigInt column, BoutReal value) {
    HYPRE_BigInt i = row - ilower;
    ASSERT2(i >= 0 && i < num_rows);
    bool value_set = false;
    for(HYPRE_BigInt col_ind = 0; col_ind < (*J)[i].size(); ++col_ind) {
      if ((*J)[i][col_ind] == column) {
        (*V)[i][col_ind] = value;
        value_set  = true;
        break;
      }
    }

    if (!value_set)
    {
      (*V)[i].push_back(value);
      (*J)[i].push_back(column);
    }
  }

  void addVal(const ind_type& row, const ind_type& column, BoutReal value) {
    const HYPRE_BigInt global_row = index_converter->getGlobal(row);
    const HYPRE_BigInt global_column = index_converter->getGlobal(column);
#if CHECKLEVEL >= 1
    if (global_row < 0 or global_column < 0) {
      throw BoutException(
          "Invalid HypreMatrix element at row = ({}, {}, {}), column = ({}, {}, {})",
          row.x(), row.y(), row.z(), column.x(), column.y(), column.z());
    }
#endif
    return addVal(global_row, global_column, value);
  }

  void addVal(const HYPRE_BigInt row, const HYPRE_BigInt column, BoutReal value) {
    HYPRE_BigInt i = row - ilower;
    ASSERT2(i >= 0 && i < num_rows);
    bool value_set = false;
    for(HYPRE_BigInt col_ind = 0; col_ind < (*J)[i].size(); ++col_ind) {
      if ((*J)[i][col_ind] == column) {
        (*V)[i][col_ind] += value;
        value_set  = true;
        break;
      }
    }

    if (!value_set)
    {
      (*V)[i].push_back(value);
      (*J)[i].push_back(column);
    }
  }

  BoutReal operator()(const ind_type& row, const ind_type& column) const {
    HYPRE_BigInt global_row = index_converter->getGlobal(row);
    HYPRE_BigInt global_column = index_converter->getGlobal(column);
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
    return getVal(row, column);
  }

  Element operator()(const ind_type& row, const ind_type& column) {
    HYPRE_BigInt global_row = index_converter->getGlobal(row);
    HYPRE_BigInt global_column = index_converter->getGlobal(column);

#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Cannot return element of uninitialised HypreMatrix");
    }
    if (global_row == -1 or global_column == -1) {
      throw BoutException(
          "Invalid HypreMatrix element at row = ({}, {}, {}), column = ({}, {}, {}), global: {} {}",
          row.x(), row.y(), row.z(), column.x(), column.y(), column.z(), global_row, global_column);
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
    HYPRE_BigInt num_entries = 0;
    HYPRE_BigInt *num_cols;
    HYPRE_BigInt *cols;
    HYPRE_BigInt *rawI;
    HYPRE_Complex *vals;

    HypreMalloc(num_cols, num_rows*sizeof(HYPRE_BigInt)); 
    for (HYPRE_BigInt i = 0; i < num_rows; ++i) {
      num_cols[i] = (*J)[i].size();;
      num_entries += (*J)[i].size();
    }

    HypreMalloc(rawI, num_rows*sizeof(HYPRE_BigInt));
    HypreMalloc(cols, num_entries*sizeof(HYPRE_BigInt)); 
    HypreMalloc(vals, num_entries*sizeof(HYPRE_Complex));
    HYPRE_BigInt entry = 0;
    for (HYPRE_BigInt i = 0; i < num_rows; ++i) {
      rawI[i] = (*I)[i];
      for(HYPRE_BigInt col_ind = 0; col_ind < num_cols[i]; ++col_ind) {
        cols[entry] = (*J)[i][col_ind];
        vals[entry] = (*V)[i][col_ind];
        entry ++;
      }
    }
    checkHypreError(HYPRE_IJMatrixSetValues(*hypre_matrix, num_rows, num_cols, rawI, cols, vals));
    checkHypreError(HYPRE_IJMatrixAssemble(*hypre_matrix));
    checkHypreError(HYPRE_IJMatrixGetObject(*hypre_matrix, reinterpret_cast<void**>(&parallel_matrix)));
    assembled = true;

    HypreFree(rawI);
    HypreFree(num_cols);
    HypreFree(cols);
    HypreFree(vals);
  }

  bool isAssembled() const { return assembled; }

  HypreMatrix<T> yup(int index = 0) { return ynext(index + 1); }
  HypreMatrix<T> ydown(int index = 0) { return ynext(-index - 1); }
  HypreMatrix<T> ynext(int dir) {
    if (std::is_same<T, FieldPerp>::value and ((yoffset + dir) != 0)) {
      throw BoutException("Can not get ynext for FieldPerp");
    }
    HypreMatrix<T> result;
    result.comm = comm;
    result.ilower = ilower;
    result.iupper = iupper;
    result.hypre_matrix = hypre_matrix;
    result.parallel_matrix = parallel_matrix;    
    result.index_converter = index_converter;
    result.location = location;
    result.initialised = initialised;
    result.yoffset = std::is_same<T, Field2D>::value ? 0 : yoffset + dir;
    result.parallel_transform = parallel_transform;
    result.assembled = assembled;
    result.num_rows = num_rows;
    result.I = I;   //We want the pointer to transfer so this works like a view
    result.J = J;
    result.V = V;

    return result;
  }

  HYPRE_IJMatrix get() { return *hypre_matrix; }
  const HYPRE_IJMatrix& get() const { return *hypre_matrix; }

  HYPRE_ParCSRMatrix getParallel() { return parallel_matrix; }
  const HYPRE_ParCSRMatrix& getParallel() const { return parallel_matrix; }

  // y = alpha*A*x + beta*y
  // Note result is returned in 'y' argument
  void computeAxpby(double alpha, HypreVector<T> &x, double beta, HypreVector<T> &y)
  {
    checkHypreError(HYPRE_ParCSRMatrixMatvec(alpha, parallel_matrix, x.getParallel(), beta, y.getParallel()));
  }

  // y = A*x
  // Note result is returned in 'y' argument
  void computeAx(HypreVector<T> &x, HypreVector<T> &y)
  {
    checkHypreError(HYPRE_ParCSRMatrixMatvec(1.0, parallel_matrix, x.getParallel(), 0.0, y.getParallel()));
  }

};

// This is a specific Hypre matrix system.  here we will be setting is the system
// AMG(Px=b) A x = b
// and we will be solvingthe system using the preconditioned conjugate gradient 
// method utilizing the AMG preconditioner on the system Px = b.  In some cases
// it will be advantage to us a P that is an approximation to A instead of P = A
template <class T>
class HypreSystem
{
private:
  HypreMatrix<T> *P{nullptr};
  HypreMatrix<T> *A{nullptr};
  HypreVector<T> *x{nullptr};
  HypreVector<T> *b{nullptr};
  MPI_Comm comm;
  HYPRE_Solver solver;
  HYPRE_Solver precon;
  bool solver_setup;
  HYPRE_SOLVER_TYPE solver_type = HYPRE_SOLVER_TYPE::gmres;
public:
  HypreSystem(Mesh& mesh, Options& options)
  {
    HYPRE_Init(); // each processor inits Hypre ; still need to work in interaction with hyprelib init/finalize

    solver_type = options["hypre_solver_type"]
                  .doc("Type of solver to use when solving Hypre system. Possible "
                       "values are: gmres")
                  .withDefault(HYPRE_SOLVER_TYPE::gmres);

#ifdef BOUT_USE_CUDA
    hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_DEVICE;
    HYPRE_MemoryLocation memory_location = HYPRE_MEMORY_DEVICE;
    hypre_HandleMemoryLocation(hypre_handle())    = memory_location;
#else
    hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_HOST;
    HYPRE_MemoryLocation memory_location = HYPRE_MEMORY_HOST;
    hypre_HandleMemoryLocation(hypre_handle())    = memory_location;
#endif

    comm = std::is_same<T, FieldPerp>::value ? mesh.getXcomm() : BoutComm::get();

    if (solver_type == HYPRE_SOLVER_TYPE::gmres) {
      HYPRE_ParCSRGMRESCreate(comm, &solver);
      /* set the GMRES parameters */
      //HYPRE_GMRESSetKDim(solver, 5); // commented out to let Hypre pick default
    } else {
      throw BoutException("Unsupported hypre_solver_type {}", solver_type);
    }

    setRelTol(options["rtol"].doc("Relative tolerance for Hypre solver")
                             .withDefault(1.0e-7));
    setAbsTol(options["atol"].doc("Relative tolerance for Hypre solver")
                             .withDefault(1.0e-12));
    setMaxIter(options["maxits"].doc("Maximum iterations for Hypre solver")
                                .withDefault(10000));

    HYPRE_BoomerAMGCreate(&precon);
    HYPRE_BoomerAMGSetOldDefault(precon);
#ifdef BOUT_USE_CUDA
    HYPRE_BoomerAMGSetRelaxType(precon, 18);  // 18 or 7 for GPU implementation // 7 throws error code 256 did not converge
    HYPRE_BoomerAMGSetRelaxOrder(precon, false); // must be false for GPU
    //HYPRE_BoomerAMGSetCoarsenType(precon, 8); // must be PMIS (8) for GPU // currently causes solver_err = 1 Generic Error and non-convergence
    HYPRE_BoomerAMGSetInterpType(precon, 15); // must be 3 or 15 for GPU 
#endif
    HYPRE_BoomerAMGSetNumSweeps(precon, 1);
    HYPRE_BoomerAMGSetMaxLevels(precon, 20);
    HYPRE_BoomerAMGSetKeepTranspose(precon, 1);
    HYPRE_BoomerAMGSetTol(precon, 0.0);
    HYPRE_BoomerAMGSetPrintLevel(
      precon,
      options["hypre_print_level"]
      .doc("Verbosity for Hypre solver. Integer from 0 (silent) to 4 (most verbose).")
      .withDefault(0)
    );

    solver_setup = false;
  }

  ~HypreSystem() {
    //HYPRE_Finalize(); // now handled by hyprelib 
  }

  void setRelTol(double tol)
  {
    if (solver_type == HYPRE_SOLVER_TYPE::gmres) {
      checkHypreError(HYPRE_ParCSRGMRESSetTol(solver, tol));
    } else {
      throw BoutException("Unsupported hypre_solver_type {}", solver_type);
    }
  }

  void setAbsTol(double tol)
  {
    if (solver_type == HYPRE_SOLVER_TYPE::gmres) {
      checkHypreError(HYPRE_ParCSRGMRESSetAbsoluteTol(solver, tol));
    } else {
      throw BoutException("Unsupported hypre_solver_type {}", solver_type);
    }
  }

  void setMaxIter(int max_iter)
  {
    if (solver_type == HYPRE_SOLVER_TYPE::gmres) {
      checkHypreError(HYPRE_ParCSRGMRESSetMaxIter(solver, max_iter));
    } else {
      throw BoutException("Unsupported hypre_solver_type {}", solver_type);
    }
  }

  double getFinalRelResNorm()
  {
    HYPRE_Real resnorm;
    if (solver_type == HYPRE_SOLVER_TYPE::gmres) {
      checkHypreError(HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &resnorm));
    } else {
      throw BoutException("Unsupported hypre_solver_type {}", solver_type);
    }
    return resnorm;
  }

  int getNumItersTaken()
  {
    HYPRE_Int iters;
    if (solver_type == HYPRE_SOLVER_TYPE::gmres) {
      checkHypreError(HYPRE_ParCSRGMRESGetNumIterations(solver, &iters));
    } else {
      throw BoutException("Unsupported hypre_solver_type {}", solver_type);
    }
    return iters;
  }

  void setMatrix(HypreMatrix<T>* A_) {
    A=A_;
  }

  int setupAMG(HypreMatrix<T>* P_) {
    P = P_;
    if (solver_type == HYPRE_SOLVER_TYPE::gmres) {
      checkHypreError(HYPRE_ParCSRGMRESSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precon));
    } else {
      throw BoutException("Unsupported hypre_solver_type {}", solver_type);
    }
    int setup_err = checkHypreError(HYPRE_BoomerAMGSetup(precon, P->getParallel(), nullptr, nullptr));

    return setup_err;
  }

  void setSolutionVector(HypreVector<T>* x_) {
    x = x_;
  }

  void setRHSVector(HypreVector<T>* b_) {
    b = b_;
  }

  HypreMatrix<T>* getMatrix() {
    return A;
  }

  HypreMatrix<T>* getAMGMatrix() {
    return P;
  }  

  HypreVector<T>* getRHS() {
    return b;
  }

  HypreVector<T>* getSolution() {
    return x;
  }

  HYPRE_PtrToParSolverFcn SetupFcn() const { 
    return (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup; 
  }
  
  HYPRE_PtrToParSolverFcn SolveFcn() const { 
    return (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve; 
  }

  int solve()
  {
    int solve_err;
    ASSERT2(A != nullptr);
    ASSERT2(x != nullptr);
    ASSERT2(b != nullptr);
    if (solver_type == HYPRE_SOLVER_TYPE::gmres) {
      if (!solver_setup) {
        checkHypreError(HYPRE_ParCSRGMRESSetup(solver, A->getParallel(), b->getParallel(), x->getParallel()));
        solver_setup = true;
      }

      solve_err = checkHypreError(HYPRE_ParCSRGMRESSolve(solver, A->getParallel(), b->getParallel(), x->getParallel()));

    } else {
      throw BoutException("Unsupported hypre_solver_type {}", solver_type);
    }

    return solve_err;
  }


};


} // namespace bout

#endif // BOUT_HAS_HYPRE

#endif // BOUT_HYPRE_INTERFACE_H
