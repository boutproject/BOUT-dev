/*!***********************************************************************
 * \file petsc_interface.hxx
 * Classes to wrap PETSc matrices and vectors, providing a convenient
 * interface to them. In particular, they will internally convert
 * between BOUT++ indices and PETSc ones, making it far easier to set
 * up a linear system.
 *
 **************************************************************************
 * Copyright 2019 C. MacMackin
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef BOUT_PETSC_INTERFACE_H
#define BOUT_PETSC_INTERFACE_H

#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include <algorithm>
#include <iterator>
#include <memory>
#include <type_traits>
#include <vector>

#include <bout/bout_types.hxx>
#include <bout/boutcomm.hxx>
#include <bout/globalindexer.hxx>
#include <bout/mesh.hxx>
#include <bout/operatorstencil.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/petsclib.hxx>
#include <bout/region.hxx>
#include <bout/traits.hxx>

#include <petscsystypes.h>
#include <petscvec.h>

/*!
 * A class which wraps PETSc vector objects, allowing them to be
 * indexed using the BOUT++ scheme. Note that boundaries are only
 * included in the vector to a depth of 1.
 */
template <class T>
class PetscVector;
template <class T>
void swap(PetscVector<T>& first, PetscVector<T>& second);

template <class T>
inline MPI_Comm getComm([[maybe_unused]] const T& field) {
  return BoutComm::get();
}

template <>
inline MPI_Comm getComm([[maybe_unused]] const FieldPerp& field) {
  return field.getMesh()->getXcomm();
}

template <class T>
class PetscVector {
public:
  static_assert(bout::utils::is_Field_v<T>, "PetscVector only works with Fields");
  using ind_type = typename T::ind_type;

  struct VectorDeleter {
    void operator()(Vec* vec) const {
      VecDestroy(vec);
      delete vec;
    }
  };

  PetscVector() : vector(new Vec{}) {}
  ~PetscVector() = default;

  PetscVector(const PetscVector<T>& vec)
      : vector(new Vec()), indexConverter(vec.indexConverter), location(vec.location),
        initialised(vec.initialised), vector_values(vec.vector_values) {
    VecDuplicate(*vec.vector, vector.get());
    VecCopy(*vec.vector, *vector);
  }

  PetscVector(PetscVector<T>&& vec) noexcept
      : indexConverter(vec.indexConverter), location(vec.location),
        initialised(vec.initialised), vector_values(std::move(vec.vector_values)) {
    std::swap(vector, vec.vector);
    vec.initialised = false;
  }

  /// Construct from a field, copying over the field values
  PetscVector(const T& field, IndexerPtr<T> indConverter)
      : vector(new Vec()), indexConverter(indConverter), location(field.getLocation()),
        initialised(true), vector_values(field.size()) {
    ASSERT1(indConverter->getMesh() == field.getMesh());
    MPI_Comm comm = getComm(field);

    const int size = indexConverter->size();
    VecCreateMPI(comm, size, PETSC_DECIDE, vector.get());
    // This allows us to pass negative indices
    VecSetOption(*vector, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
    *this = field;
  }

  /// Construct a vector like v, but using data from a raw PETSc
  /// Vec. That Vec (not a copy) will then be owned by the new object.
  PetscVector(const PetscVector<T>& other, Vec* vec) {
#if CHECKLEVEL >= 2
    int fsize = other.indexConverter->size();
    int msize = 0;
    VecGetLocalSize(*vec, &msize);
    ASSERT2(fsize == msize);
#endif
    vector.reset(vec);
    indexConverter = other.indexConverter;
    location = other.location;
    initialised = true;
  }

  /// Copy assignment
  PetscVector<T>& operator=(const PetscVector<T>& rhs) {
    *this = PetscVector<T>(rhs);
    return *this;
  }

  /// Move assignment
  PetscVector<T>& operator=(PetscVector<T>&& rhs) noexcept {
    swap(*this, rhs);
    return *this;
  }

  PetscVector<T>& operator=(const T& field) {
    ASSERT1(indexConverter); // Needs to have index set
    BOUT_FOR_SERIAL(i, indexConverter->getRegionAll()) { (*this)(i) = field[i]; }
    assemble();
    return *this;
  }

  friend void swap<T>(PetscVector<T>& first, PetscVector<T>& second);

  BoutReal& operator()(const ind_type& index) {
#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Can not return element of uninitialised vector");
    }
#endif
    return vector_values[index.ind];
  }

  BoutReal operator()(const ind_type& index) const {
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
    BoutReal value = BoutNaN;
    int status = 0;
    BOUT_OMP_SAFE(critical)
    status = VecGetValues(*get(), 1, &global, &value);
    if (status != 0) {
      throw BoutException("Error when getting element of a PETSc vector.");
    }
    return value;
  }

  void assemble() {
    auto ierr = VecSetValues(*vector, vector_values.size(),
                             indexConverter->getIntIndices().begin(),
                             vector_values.begin(), INSERT_VALUES);
    if (ierr < 0) {
      throw BoutException("PETSc error");
    }

    ierr = VecAssemblyBegin(*vector);
    if (ierr < 0) {
      throw BoutException("PETSc error");
    }
    ierr = VecAssemblyEnd(*vector);
    if (ierr < 0) {
      throw BoutException("PETSc error");
    }
  }

  void destroy() {
    if (*vector != nullptr && initialised) {
      VecDestroy(vector.get());
      *vector = nullptr;
      initialised = false;
    }
  }

  /// Returns a field constructed from the contents of this vector
  T toField() const {
    T result(indexConverter->getMesh());
    result.allocate();
    result.setLocation(location);
    // Note that this only populates boundaries to a depth of 1
    BOUT_FOR_SERIAL(i, indexConverter->getRegionAll()) {
      const PetscInt ind = indexConverter->getGlobal(i);
      PetscScalar val = BoutNaN;
      VecGetValues(*vector, 1, &ind, &val);
      result[i] = val;
    }
    return result;
  }

  /// Provides a reference to the raw PETSc Vec object.
  Vec* get() { return vector.get(); }
  const Vec* get() const { return vector.get(); }

private:
  PetscLib lib{};
  std::unique_ptr<Vec, VectorDeleter> vector = nullptr;
  IndexerPtr<T> indexConverter{};
  CELL_LOC location = CELL_LOC::deflt;
  bool initialised = false;
  Array<BoutReal> vector_values{};
};

/*!
 * A class which wraps PETSc vector objects, allowing them to be
 * indexed using the BOUT++ scheme. It provides the option of setting
 * a y-offset that interpolates onto field lines.
 */
template <class T>
class PetscMatrix;
template <class T>
void swap(PetscMatrix<T>& first, PetscMatrix<T>& second);

template <class T>
class PetscMatrix {
public:
  static_assert(bout::utils::is_Field_v<T>, "PetscMatrix only works with Fields");
  using ind_type = typename T::ind_type;

  struct MatrixDeleter {
    void operator()(Mat* mat) const {
      MatDestroy(mat);
      delete mat;
    }
  };

  /// Default constructor does nothing
  PetscMatrix() : matrix(new Mat()) {}
  ~PetscMatrix() = default;

  /// Copy constructor
  PetscMatrix(const PetscMatrix<T>& mat)
      : matrix(new Mat()), indexConverter(mat.indexConverter), pt(mat.pt),
        yoffset(mat.yoffset), initialised(mat.initialised) {
    MatDuplicate(*mat.matrix, MAT_COPY_VALUES, matrix.get());
  }

  /// Move constrcutor
  PetscMatrix(PetscMatrix<T>&& mat) noexcept
      : matrix(std::move(mat.matrix)), indexConverter(std::move(mat.indexConverter)),
        pt(mat.pt), yoffset(mat.yoffset), initialised(mat.initialised) {
    mat.initialised = false;
  }

  // Construct a matrix capable of operating on the specified field,
  // preallocating memory if requeted and possible.
  PetscMatrix(IndexerPtr<T> indConverter, bool preallocate = true)
      : matrix(new Mat()), indexConverter(indConverter),
        pt(&indConverter->getMesh()->getCoordinates()->getParallelTransform()) {
    MPI_Comm comm = std::is_same_v<T, FieldPerp> ? indConverter->getMesh()->getXcomm()
                                                 : BoutComm::get();
    const int size = indexConverter->size();

    MatCreate(comm, matrix.get());
    MatSetSizes(*matrix, size, size, PETSC_DECIDE, PETSC_DECIDE);
    MatSetType(*matrix, MATMPIAIJ);

    // If a stencil has been provided, preallocate memory
    if (preallocate && indexConverter->sparsityPatternAvailable()) {
      MatMPIAIJSetPreallocation(*matrix, 0, indexConverter->getNumDiagonal().data(), 0,
                                indexConverter->getNumOffDiagonal().data());
    }

    MatSetUp(*matrix);
    yoffset = 0;
    initialised = true;
  }

  /// Copy assignment
  PetscMatrix<T>& operator=(PetscMatrix<T> rhs) {
    swap(*this, rhs);
    return *this;
  }
  /// Move assignment
  PetscMatrix<T>& operator=(PetscMatrix<T>&& rhs) noexcept {
    matrix = rhs.matrix;
    indexConverter = rhs.indexConverter;
    pt = rhs.pt;
    yoffset = rhs.yoffset;
    initialised = rhs.initialised;
    rhs.initialised = false;
    return *this;
  }
  friend void swap<T>(PetscMatrix<T>& first, PetscMatrix<T>& second);

  /*!
   * A class which is used to assign to a particular element of a PETSc
   * matrix, potentially with a y-offset. It is meant to be transient
   * and will be destroyed immediately after use. In general you should
   * not try to assign an instance to a variable.
   *
   * The Element object will store a copy of the value which has been
   * assigned to it. This is because mixing calls to the PETSc
   * function MatGetValues() and MatSetValues() without intermediate
   * matrix assembly will cause errors. Thus, if the user wishes to
   * get the value of the matrix, it must be stored here.
   */
  class Element {
  public:
    Element() = delete;
    ~Element() = default;
    Element(Element&&) noexcept = default;
    Element& operator=(Element&&) noexcept = default;
    Element(const Element& other) = default;
    Element(Mat* matrix, PetscInt row, PetscInt col, std::vector<PetscInt> position = {},
            std::vector<BoutReal> weight = {})
        : petscMatrix(matrix), petscRow(row), petscCol(col),
          positions(std::move(position)), weights(std::move(weight)) {
      ASSERT2(positions.size() == weights.size());
#if CHECK > 2
      for (const auto val : weights) {
        ASSERT3(finite(val));
      }
#endif
      if (positions.empty()) {
        positions = {col};
        weights = {1.0};
      }
      PetscBool assembled = PETSC_FALSE;
      MatAssembled(*petscMatrix, &assembled);
      if (assembled == PETSC_TRUE) {
        BOUT_OMP_SAFE(critical)
        MatGetValues(*petscMatrix, 1, &petscRow, 1, &petscCol, &value);
      } else {
        value = 0.;
      }
    }
    Element& operator=(const Element& other) {
      AUTO_TRACE();
      if (this == &other) {
        return *this;
      }
      ASSERT3(finite(static_cast<BoutReal>(other)));
      *this = static_cast<BoutReal>(other);
      return *this;
    }
    Element& operator=(BoutReal val) {
      AUTO_TRACE();
      ASSERT3(finite(val));
      value = val;
      setValues(val, INSERT_VALUES);
      return *this;
    }
    Element& operator+=(BoutReal val) {
      AUTO_TRACE();
      ASSERT3(finite(val));
      auto columnPosition = std::find(positions.begin(), positions.end(), petscCol);
      if (columnPosition != positions.end()) {
        const int index = std::distance(positions.begin(), columnPosition);
        value += weights[index] * val;
        ASSERT3(finite(value));
      }
      setValues(val, ADD_VALUES);
      return *this;
    }
    operator BoutReal() const { return value; }

  private:
    void setValues(BoutReal val, InsertMode mode) {
      TRACE("PetscMatrix setting values at ({}, {})", petscRow, petscCol);
      ASSERT3(positions.size() > 0);
      std::vector<PetscScalar> values;
      std::transform(weights.begin(), weights.end(), std::back_inserter(values),
                     [&val](BoutReal weight) -> PetscScalar { return weight * val; });

      int status = 0;
      BOUT_OMP_SAFE(critical)
      status = MatSetValues(*petscMatrix, 1, &petscRow, positions.size(),
                            positions.data(), values.data(), mode);
      if (status != 0) {
        throw BoutException("Error when setting elements of a PETSc matrix.");
      }
    }
    Mat* petscMatrix;
    PetscInt petscRow, petscCol;
    PetscScalar value;
    std::vector<PetscInt> positions;
    std::vector<BoutReal> weights;
  };

  Element operator()(const ind_type& index1, const ind_type& index2) {
    const int global1 = indexConverter->getGlobal(index1);
    const int global2 = indexConverter->getGlobal(index2);
#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Can not return element of uninitialised matrix");
    }
    if (global1 == -1 || global2 == -1) {
      throw BoutException(
          "Request to return invalid matrix element: (({}, {}, {}), ({}, {}, {}))",
          index1.x(), index1.y(), index1.z(), index2.x(), index2.y(), index2.z());
    }
#endif
    std::vector<PetscInt> positions;
    std::vector<PetscScalar> weights;
    if (yoffset != 0) {
      ASSERT1(yoffset == index2.y() - index1.y());
      const auto pws =
          pt->getWeightsForYApproximation(index2.x(), index1.y(), index2.z(), yoffset);
      const int ny =
          std::is_same_v<T, FieldPerp> ? 1 : indexConverter->getMesh()->LocalNy;
      const int nz = std::is_same_v<T, Field2D> ? 1 : indexConverter->getMesh()->LocalNz;

      std::transform(
          pws.begin(), pws.end(), std::back_inserter(positions),
          [this, ny, nz](ParallelTransform::PositionsAndWeights pos) -> PetscInt {
            return this->indexConverter->getGlobal(
                ind_type(pos.i * ny * nz + pos.j * nz + pos.k, ny, nz));
          });
      std::transform(pws.begin(), pws.end(), std::back_inserter(weights),
                     [](ParallelTransform::PositionsAndWeights pos) -> PetscScalar {
                       return pos.weight;
                     });
    }
    return Element(matrix.get(), global1, global2, positions, weights);
  }

  BoutReal operator()(const ind_type& index1, const ind_type& index2) const {
    ASSERT2(yoffset == 0);
    const int global1 = indexConverter->getGlobal(index1);
    const int global2 = indexConverter->getGlobal(index2);
#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Can not return element of uninitialised matrix");
    }
    if (global1 == -1 || global2 == -1) {
      throw BoutException(
          "Request to return invalid matrix element: (({}, {}, {}), ({}, {}, {}))",
          index1.x(), index1.y(), index1.z(), index2.x(), index2.y(), index2.z());
    }
#endif
    BoutReal value = BoutNaN;
    int status = 0;
    BOUT_OMP_SAFE(critical)
    status = MatGetValues(*get(), 1, &global1, 1, &global2, &value);
    if (status != 0) {
      throw BoutException("Error when getting elements of a PETSc matrix.");
    }
    return value;
  }

  // Assemble the matrix prior to use
  void assemble() {
    MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
  }

  // Partially assemble the matrix so you can switch between adding
  // and inserting values
  void partialAssemble() {
    MatAssemblyBegin(*matrix, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*matrix, MAT_FLUSH_ASSEMBLY);
  }

  void destroy() {
    if (*matrix != nullptr && initialised) {
      MatDestroy(matrix.get());
      *matrix = nullptr;
      initialised = false;
    }
  }

  PetscMatrix<T> yup(int index = 0) { return ynext(index + 1); }
  PetscMatrix<T> ydown(int index = 0) { return ynext(-index - 1); }
  PetscMatrix<T> ynext(int dir) {
    if (std::is_same_v<T, FieldPerp> && yoffset + dir != 0) {
      throw BoutException("Can not get ynext for FieldPerp");
    }
    PetscMatrix<T> result; // Can't use copy constructor because don't
                           // want to duplicate the matrix
    result.matrix = matrix;
    result.indexConverter = indexConverter;
    result.pt = pt;
    result.yoffset = std::is_same_v<T, Field2D> ? 0 : yoffset + dir;
    result.initialised = initialised;
    return result;
  }

  /// Provides a reference to the raw PETSc Mat object.
  Mat* get() { return matrix.get(); }
  const Mat* get() const { return matrix.get(); }

private:
  PetscLib lib;
  std::shared_ptr<Mat> matrix = nullptr;
  IndexerPtr<T> indexConverter{};
  ParallelTransform* pt{};
  int yoffset = 0;
  bool initialised = false;
};

/*!
 * Move reference to one Vec from \p first to \p second and
 * vice versa.
 */
template <class T>
void swap(PetscVector<T>& first, PetscVector<T>& second) {
  std::swap(first.vector, second.vector);
  std::swap(first.indexConverter, second.indexConverter);
  std::swap(first.location, second.location);
  std::swap(first.initialised, second.initialised);
}

/*!
 * Move reference to one Mat from \p first to \p second and
 * vice versa.
 */
template <class T>
void swap(PetscMatrix<T>& first, PetscMatrix<T>& second) {
  std::swap(first.matrix, second.matrix);
  std::swap(first.indexConverter, second.indexConverter);
  std::swap(first.pt, second.pt);
  std::swap(first.yoffset, second.yoffset);
  std::swap(first.initialised, second.initialised);
}

/*!
 *  Performs matrix-multiplication on the supplied vector
 */
template <class T>
PetscVector<T> operator*(const PetscMatrix<T>& mat, const PetscVector<T>& vec) {
  Vec rhs = *vec.get();
  Vec* result = new Vec();
  VecDuplicate(rhs, result);
  VecAssemblyBegin(*result);
  VecAssemblyEnd(*result);
  const int err = MatMult(*mat.get(), rhs, *result);
  ASSERT2(err == 0);
  return PetscVector<T>(vec, result);
}

// Compatibility wrappers
// For < 3.24
#if PETSC_VERSION_GE(3, 24, 0) \
    || (PETSC_VERSION_GE(3, 23, 0) && PETSC_VERSION_RELEASE == 0)
namespace bout {
template <class T>
constexpr auto cast_MatFDColoringFn(T func) {
  return func;
}
} // namespace bout
#else
using MatFDColoringFn = PetscErrorCode (*)();
namespace bout {
template <class T>
constexpr auto cast_MatFDColoringFn(T func) {
  return reinterpret_cast<MatFDColoringFn>(func); // NOLINT(*-reinterpret-cast)
}
} // namespace bout
#endif

#endif // BOUT_HAS_PETSC

#endif // BOUT_PETSC_INTERFACE_H
