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

#ifndef __PETSC_INTERFACE_H__
#define __PETSC_INTERFACE_H__

#include <algorithm>
#include <memory>
#include <type_traits>
#include <vector>

#include <bout/mesh.hxx>
#include <bout/operatorstencil.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/petsclib.hxx>
#include <bout/region.hxx>
#include <bout/traits.hxx>
#include <bout_types.hxx>
#include <boutcomm.hxx>

#ifdef BOUT_HAS_PETSC
template <class T>
class GlobalIndexer;

template <class T>
using IndexerPtr = std::shared_ptr<GlobalIndexer<T>>;

using InterpolationWeights = std::vector<ParallelTransform::PositionsAndWeights>;

/*!
 * A singleton which accepts index objects produced by iterating over
 * fields and returns a global index. This index can be used when
 * constructing PETSc arrays. Guard regions used for communication
 * between processes will have the indices of the part of the interior
 * region they are mirroring. Boundaries required by the stencil will
 * have unique indices. If no stencil is provided then only the
 * interior region will be assigned global indices. By default, the
 * indexer is fully initialised so that guard cells are communicated
 * (ensuring they hold the appropriate global indices). However, by
 * passing `autoInitialise = false` this behaviour can be turned off
 * and the user can then manually call the `initialise()` method
 * later. This can be useful for mocking/faking the class when
 * testing.
 */
template <class T>
class GlobalIndexer {
public:
  static_assert(bout::utils::is_Field<T>::value, "GlobalIndexer only works with Fields");
  using ind_type = typename T::ind_type;

  GlobalIndexer() = default;

  GlobalIndexer(Mesh* localmesh,
                OperatorStencil<ind_type> stencil = OperatorStencil<ind_type>(),
                bool autoInitialise = true)
      : fieldmesh(localmesh), indices(-1., localmesh), stencils(std::move(stencil)) {

    Region<ind_type> allCandidate, bndryCandidate;
    if (stencils.getNumParts() > 0) {
      std::set<ind_type> allIndices(getRegionNobndry().getIndices().begin(),
                                    getRegionNobndry().getIndices().end()),
          newIndices;
      BOUT_FOR_SERIAL(i, getRegionNobndry()) {
        for (const IndexOffset<ind_type> j : stencils.getStencilPart(i)) {
          insertIndex(i + j, allIndices, newIndices);
        }
      }
      std::set<ind_type> candidateIndices = newIndices;
      while (candidateIndices.size() > 0) {
        newIndices.clear();
        for (const ind_type i : candidateIndices) {
          insertIndex(i, allIndices, newIndices);
        }
        candidateIndices = newIndices;
      }
      std::vector<ind_type> allIndicesVec(allIndices.begin(), allIndices.end());
      allCandidate = Region<ind_type>(allIndicesVec);
    }

    bndryCandidate = mask(allCandidate, getRegionNobndry());

    regionInnerX = getUnion(bndryCandidate, indices.getRegion("RGN_INNER_X"));
    regionOuterX = getUnion(bndryCandidate, indices.getRegion("RGN_OUTER_X"));
    if (std::is_same<T, FieldPerp>::value) {
      regionLowerY = Region<ind_type>({});
      regionUpperY = Region<ind_type>({});
    } else {
      regionLowerY = getUnion(bndryCandidate, indices.getRegion("RGN_LOWER_Y"));
      regionUpperY = getUnion(bndryCandidate, indices.getRegion("RGN_UPPER_Y"));
    }
    regionBndry = regionLowerY + regionInnerX + regionOuterX + regionUpperY;
    regionAll = getRegionNobndry() + regionBndry;
    regionBndry.sort();
    regionAll.sort();

    int localSize = size();
    MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? fieldmesh->getXcomm() : BoutComm::get();
    fieldmesh->getMpi().MPI_Scan(&localSize, &globalEnd, 1, MPI_INT, MPI_SUM, comm);
    globalEnd--;
    int counter = globalStart = globalEnd - size() + 1;

    BOUT_FOR_SERIAL(i, regionAll) { indices[i] = counter++; }

    if (autoInitialise) {
      initialise();
    }
  }

  virtual ~GlobalIndexer() {}
  
  /// Call this immediately after construction when running unit tests.
  void initialiseTest() {}

  /// Finish setting up the indexer, communicating indices across
  /// processes and, if possible, calculating the sparsity pattern of
  /// any matrices.
  void initialise() {
    fieldmesh->communicate(indices);
  }

  Mesh* getMesh() const { return fieldmesh; }

  /// Convert the local index object to a global index which can be
  /// used in PETSc vectors and matrices.
  PetscInt getGlobal(const ind_type& ind) const {
    return static_cast<PetscInt>(std::round(indices[ind]));
  }

  /// Check whether the local index corresponds to an element which is
  /// stored locally.
  bool isLocal(const ind_type& ind) const {
    if (ind.ind < 0) {
      return false;
    }
    PetscInt index = getGlobal(ind);
    return (globalStart <= index) && (index <= globalEnd);
  }

  PetscInt getGlobalStart() const { return globalStart; }

  const Region<ind_type>& getRegionAll() const { return regionAll; }
  const Region<ind_type>& getRegionNobndry() const {
    return indices.getRegion("RGN_NOBNDRY");
  }
  const Region<ind_type>& getRegionBndry() const { return regionBndry; }
  const Region<ind_type>& getRegionLowerY() const { return regionLowerY; }
  const Region<ind_type>& getRegionUpperY() const { return regionUpperY; }
  const Region<ind_type>& getRegionInnerX() const { return regionInnerX; }
  const Region<ind_type>& getRegionOuterX() const { return regionOuterX; }

  bool sparsityPatternAvailable() const { return stencils.getNumParts() > 0; }

  const std::vector<int>& getNumDiagonal() const {
    ASSERT2(sparsityPatternAvailable());
    if (!sparsityCalculated) {
      calculateSparsity();
    }
    return numDiagonal;
  }

  const std::vector<int>& getNumOffDiagonal() const {
    ASSERT2(sparsityPatternAvailable());
    if (!sparsityCalculated) {
      calculateSparsity();
    }
    return numOffDiagonal;
  }

  int size() const { return regionAll.size(); }

protected:
  // Must not be const as the index field needs to be mutable in order
  // to fake parallel communication in the unit tests.
  T& getIndices() { return indices; }

private:
  static void insertIndex(const ind_type i, std::set<ind_type>& allInds,
                          std::set<ind_type>& newInds) {
    auto result = allInds.insert(i);
    if (result.second) {
      newInds.insert(i);
    }
  }

  /// This gets called by initialiseTest and is used to register
  /// fields with fake parallel meshes.
  virtual void registerFieldForTest(T& UNUSED(f)) {
    // This is a place-holder which does nothing. It can be overridden
    // by descendent classes if necessary to set up testing.
    return;
  }

  void calculateSparsity() const {
    numDiagonal = std::vector<int>(size());
    numOffDiagonal = std::vector<int>(size(), 0);

    // Set initial guess for number of on-diagonal elements
    BOUT_FOR(i, regionAll) {
      numDiagonal[getGlobal(i) - globalStart] = stencils.getStencilSize(i);
    }

    BOUT_FOR_SERIAL(i, indices.getRegion("RGN_GUARDS")) {
      if (getGlobal(i) >= 0 && !isLocal(i)) {
        for (const auto& j : stencils.getIndicesWithStencilIncluding(i)) {
          if (isLocal(j)) {
            const int n = getGlobal(j) - globalStart;
            numDiagonal[n] -= 1;
            numOffDiagonal[n] += 1;
          }
        }
      }
    }

    sparsityCalculated = true;
  }

  Mesh* fieldmesh;

  /// Fields containing the indices for each element (as reals)
  T indices;
  /// The first and last global index on this processor (inclusive in
  /// both cases)
  PetscInt globalStart, globalEnd;

  /// Stencil for which this indexer has been configured
  OperatorStencil<ind_type> stencils;

  /// Regions containing the elements for which there are global indices
  Region<ind_type> regionAll, regionLowerY, regionUpperY, regionInnerX, regionOuterX,
      regionBndry;

  mutable bool sparsityCalculated = false;
  mutable std::vector<PetscInt> numDiagonal, numOffDiagonal;
};

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
class PetscVector {
public:
  static_assert(bout::utils::is_Field<T>::value, "PetscVector only works with Fields");
  using ind_type = typename T::ind_type;

  struct VectorDeleter {
    void operator()(Vec* v) const {
      VecDestroy(v);
      delete v;
    }
  };

  /// Default constructor does nothing
  PetscVector() : vector(new Vec(), VectorDeleter()) {}

  /// Copy constructor
  PetscVector(const PetscVector<T>& v) : vector(new Vec(), VectorDeleter()) {
    VecDuplicate(*v.vector, vector.get());
    VecCopy(*v.vector, *vector);
    indexConverter = v.indexConverter;
    location = v.location;
    initialised = v.initialised;
  }

  /// Move constrcutor
  PetscVector(PetscVector<T>&& v) {
    std::swap(vector, v.vector);
    indexConverter = v.indexConverter;
    location = v.location;
    initialised = v.initialised;
    v.initialised = false;
  }

  /// Construct from a field, copying over the field values
  PetscVector(const T& f, IndexerPtr<T> indConverter)
      : vector(new Vec(), VectorDeleter()), indexConverter(indConverter) {
    ASSERT1(indConverter->getMesh() == f.getMesh());
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? f.getMesh()->getXcomm() : BoutComm::get();
    const int size = indexConverter->size();
    VecCreateMPI(comm, size, PETSC_DECIDE, vector.get());
    location = f.getLocation();
    initialised = true;
    *this = f;
  }

  /// Construct a vector like v, but using data from a raw PETSc
  /// Vec. That Vec (not a copy) will then be owned by the new object.
  PetscVector(const PetscVector<T>& v, Vec* vec) {
#if CHECKLEVEL >= 2
    int fsize = v.indexConverter->size(), msize;
    VecGetLocalSize(*vec, &msize);
    ASSERT2(fsize == msize);
#endif
    vector.reset(vec);
    indexConverter = v.indexConverter;
    location = v.location;
    initialised = true;
  }

  /// Copy assignment
  PetscVector<T>& operator=(const PetscVector<T>& rhs) {
    return *this = PetscVector<T>(rhs);
  }

  /// Move assignment
  PetscVector<T>& operator=(PetscVector<T>&& rhs) {
    swap(*this, rhs);
    return *this;
  }

  PetscVector<T>& operator=(const T& f) {
    BOUT_FOR_SERIAL(i, indexConverter->getRegionAll()) {
      const PetscInt ind = indexConverter->getGlobal(i);
      if (ind != -1) {
        // TODO: consider how VecSetValues() could be used where there
        // are continuous stretches of field data which should be
        // copied into the vector.
        VecSetValue(*vector, ind, f[i], INSERT_VALUES);
      }
    }
    assemble();
    return *this;
  }

  friend void swap<T>(PetscVector<T>& first, PetscVector<T>& second);

  /*!
   * A class which is used to assign to a particular element of a PETSc
   * vector. It is meant to be transient and will be destroyed immediately
   * after use. In general you should not try to assign an instance to a
   * variable.
   *
   * The Element object will store a copy of the value which has been
   * assigned to it. This is because mixing calls to the PETSc
   * function VecGetValues() and VecSetValues() without intermediate
   * vector assembly will cause errors. Thus, if the user wishes to
   * get the value of the vector, it must be stored here.
   */
  class Element {
  public:
    Element() = delete;
    Element(const Element& other) = default;
    Element(Vec* vector, int index) : petscVector(vector), petscIndex(index) {
      int status;
      BOUT_OMP(critical) status = VecGetValues(*petscVector, 1, &petscIndex, &value);
      if (status != 0) {
        value = 0.;
      }
    }
    Element& operator=(Element& other) { return *this = static_cast<BoutReal>(other); }
    Element& operator=(BoutReal val) {
      value = val;
      int status;
      BOUT_OMP(critical)
      status = VecSetValue(*petscVector, petscIndex, val, INSERT_VALUES);
      if (status != 0) {
        throw BoutException("Error when setting elements of a PETSc vector.");
      }
      return *this;
    }
    Element& operator+=(BoutReal val) {
      value += val;
      int status;
      BOUT_OMP(critical)
      status = VecSetValue(*petscVector, petscIndex, val, ADD_VALUES);
      if (status != 0) {
        throw BoutException("Error when setting elements of a PETSc vector.");
      }
      return *this;
    }
    operator BoutReal() const { return value; }

  private:
    Vec* petscVector = nullptr;
    PetscInt petscIndex;
    PetscScalar value;
  };

  Element operator()(const ind_type& index) {
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
    return Element(vector.get(), global);
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
    BoutReal value;
    int status;
    BOUT_OMP(critical) status = VecGetValues(*get(), 1, &global, &value);
    if (status != 0) {
      throw BoutException("Error when getting element of a PETSc vector.");
    }
    return value;
  }

  void assemble() {
    VecAssemblyBegin(*vector);
    VecAssemblyEnd(*vector);
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
      PetscInt ind = indexConverter->getGlobal(i);
      PetscScalar val;
      VecGetValues(*vector, 1, &ind, &val);
      result[i] = val;
    }
    return result;
  }

  /// Provides a reference to the raw PETSc Vec object.
  Vec* get() { return vector.get(); }
  const Vec* get() const { return vector.get(); }

private:
  PetscLib lib;
  std::unique_ptr<Vec, VectorDeleter> vector = nullptr;
  IndexerPtr<T> indexConverter;
  CELL_LOC location;
  bool initialised = false;
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
  static_assert(bout::utils::is_Field<T>::value, "PetscMatrix only works with Fields");
  using ind_type = typename T::ind_type;

  struct MatrixDeleter {
    void operator()(Mat* m) const {
      MatDestroy(m);
      delete m;
    }
  };

  PetscMatrix() : matrix(new Mat(), MatrixDeleter()) {}

  PetscMatrix(const PetscMatrix<T>& m) : matrix(new Mat(), MatrixDeleter()), pt(m.pt) {
    MatDuplicate(*m.matrix, MAT_COPY_VALUES, matrix.get());
    indexConverter = m.indexConverter;
    yoffset = m.yoffset;
    initialised = m.initialised;
  }

  PetscMatrix(PetscMatrix<T>&& m)
      : matrix(std::move(m.matrix)), indexConverter(std::move(m.indexConverter)),
        pt(std::move(m.pt)), yoffset(m.yoffset), initialised(m.initialised) {
    m.initialised = false;
  }

  // Construct a matrix capable of operating on the specified field,
  // preallocating memory if requeted and possible.
  PetscMatrix(IndexerPtr<T> indConverter, bool preallocate = true)
      : matrix(new Mat(), MatrixDeleter()), indexConverter(indConverter) {
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? indConverter->getMesh()->getXcomm() : BoutComm::get();
    pt = &indConverter->getMesh()->getCoordinates()->getParallelTransform();
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

  PetscMatrix<T>& operator=(PetscMatrix<T>&& rhs) {
    matrix = std::move(rhs.matrix);
    indexConverter = std::move(rhs.indexConverter);
    pt = std::move(rhs.pt);
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
    Element(const Element& other) = default;
    Element(Mat* matrix, PetscInt row, PetscInt col, std::vector<PetscInt> p = {},
            std::vector<BoutReal> w = {})
        : petscMatrix(matrix), petscRow(row), petscCol(col), positions(p), weights(w) {
      ASSERT2(positions.size() == weights.size());
      if (positions.size() == 0) {
        positions = {col};
        weights = {1.0};
      }
      PetscBool assembled;
      MatAssembled(*petscMatrix, &assembled);
      if (assembled == PETSC_TRUE) {
        BOUT_OMP(critical)
        MatGetValues(*petscMatrix, 1, &petscRow, 1, &petscCol, &value);
      } else {
        value = 0.;
      }
    }
    Element& operator=(Element& other) { return *this = static_cast<BoutReal>(other); }
    Element& operator=(BoutReal val) {
      value = val;
      setValues(val, INSERT_VALUES);
      return *this;
    }
    Element& operator+=(BoutReal val) {
      auto columnPosition = std::find(positions.begin(), positions.end(), petscCol);
      if (columnPosition != positions.end()) {
        int i = std::distance(positions.begin(), columnPosition);
        value += weights[i] * val;
      }
      setValues(val, ADD_VALUES);
      return *this;
    }
    operator BoutReal() const { return value; }

  private:
    void setValues(BoutReal val, InsertMode mode) {
      ASSERT3(positions.size() > 0);
      std::vector<PetscScalar> values;
      std::transform(weights.begin(), weights.end(), std::back_inserter(values),
                     [&val](BoutReal weight) -> PetscScalar { return weight * val; });
      int status;
      BOUT_OMP(critical)
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
    const int global1 = indexConverter->getGlobal(index1),
              global2 = indexConverter->getGlobal(index2);
#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Can not return element of uninitialised matrix");
    } else if (global1 == -1 || global2 == -1) {
      std::cout << "(" << index1.x() << ", " << index1.y() << ", " << index1.z() << ")  ";
      std::cout << "(" << index2.x() << ", " << index2.y() << ", " << index2.z() << ")\n";
      throw BoutException("Request to return invalid matrix element");
    }
#endif
    std::vector<PetscInt> positions;
    std::vector<PetscScalar> weights;
    if (yoffset != 0) {
      ASSERT1(yoffset == index2.y() - index1.y());
      const auto pw = [this, &index1, &index2]() {
        if (this->yoffset == -1) {
          return pt->getWeightsForYDownApproximation(index2.x(), index1.y(), index2.z());
        } else if (this->yoffset == 1) {
          return pt->getWeightsForYUpApproximation(index2.x(), index1.y(), index2.z());
        } else {
          return pt->getWeightsForYApproximation(index2.x(), index1.y(), index2.z(),
                                                 this->yoffset);
        }
      }();

      const int ny =
          std::is_same<T, FieldPerp>::value ? 1 : indexConverter->getMesh()->LocalNy;
      const int nz =
          std::is_same<T, Field2D>::value ? 1 : indexConverter->getMesh()->LocalNz;

      std::transform(
          pw.begin(), pw.end(), std::back_inserter(positions),
          [this, ny, nz](ParallelTransform::PositionsAndWeights p) -> PetscInt {
            return this->indexConverter->getGlobal(
                ind_type(p.i * ny * nz + p.j * nz + p.k, ny, nz));
          });
      std::transform(pw.begin(), pw.end(), std::back_inserter(weights),
                     [](ParallelTransform::PositionsAndWeights p) -> PetscScalar {
                       return p.weight;
                     });
    }
    return Element(matrix.get(), global1, global2, positions, weights);
  }

  BoutReal operator()(const ind_type& index1, const ind_type& index2) const {
    ASSERT2(yoffset == 0);
    const int global1 = indexConverter->getGlobal(index1),
              global2 = indexConverter->getGlobal(index2);
#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Can not return element of uninitialised matrix");
    } else if (global1 == -1 || global2 == -1) {
      std::cout << "(" << index1.x() << ", " << index1.y() << ", " << index1.z() << ")  ";
      std::cout << "(" << index2.x() << ", " << index2.y() << ", " << index2.z() << ")\n";
      throw BoutException("Request to return invalid matrix element");
    }
#endif
    BoutReal value;
    int status;
    BOUT_OMP(critical)
    status = MatGetValues(*get(), 1, &global1, 1, &global2, &value);
    if (status != 0) {
      throw BoutException("Error when setting elements of a PETSc matrix.");
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
    if (std::is_same<T, FieldPerp>::value && yoffset + dir != 0) {
      throw BoutException("Can not get ynext for FieldPerp");
    }
    PetscMatrix<T> result; // Can't use copy constructor because don't
                           // want to duplicate the matrix
    result.matrix = matrix;
    result.indexConverter = indexConverter;
    result.pt = pt;
    result.yoffset = std::is_same<T, Field2D>::value ? 0 : yoffset + dir;
    result.initialised = initialised;
    return result;
  }

  /// Provides a reference to the raw PETSc Mat object.
  Mat* get() { return matrix.get(); }
  const Mat* get() const { return matrix.get(); }

private:
  PetscLib lib;
  std::shared_ptr<Mat> matrix = nullptr;
  IndexerPtr<T> indexConverter;
  ParallelTransform* pt;
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
  const Vec rhs = *vec.get();
  Vec* result = new Vec();
  VecDuplicate(rhs, result);
  VecAssemblyBegin(*result);
  VecAssemblyEnd(*result);
  const int err = MatMult(*mat.get(), rhs, *result);
  ASSERT2(err == 0);
  return PetscVector<T>(vec, result);
}

#endif // BOUT_HAS_PETSC

#endif // __PETSC_INTERFACE_H__
