/*!***********************************************************************
 * \file petsc_interface.hxx
 * Classes to wrap PETSc matrices and vectors, providing a convenient
 * interface to them. In particular, they will internally convert
 * between BOUT++ indices and PETSc ones, making it far easier to set
 * up a linear system.
 *
 **************************************************************************
 * Copyright 2013 J. Buchanan, J.Omotani
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

#include <vector>
#include <memory>
#include <type_traits>

#include <bout_types.hxx>
#include <bout/petsclib.hxx>
#include <bout/region.hxx>
#include <bout/mesh.hxx>
#include <boutcomm.hxx>
#include <bout/paralleltransform.hxx>

using namespace std;

#ifdef BOUT_HAS_PETSC
class GlobalIndexer;
using IndexerPtr = shared_ptr<GlobalIndexer>;
using InterpolationWeights = std::vector<ParallelTransform::positionsAndWeights>;


/*!
 * A singleton which accepts index objects produced by iterating over
 * fields and returns a global index. This index can be used when
 * constructing PETSc arrays. Guard regions used for communication
 * between processes will have the indices of the part of the interior
 * region they are mirroring.
 */
class GlobalIndexer {
public:
  /// If \p localmesh is the same as the global one, return a pointer
  /// to the global instance. Otherwise create a new one.
  static IndexerPtr getInstance(Mesh* localmesh);
  /// Call this immediately after construction when running unit tests.
  void initialiseTest();
  /// Finish setting up the indexer, communicating indices across processes.
  void initialise();
  
  /// Convert the local index object to a global index which can be
  /// used in PETSc vectors and matrices.
  PetscInt getGlobal(Ind2D ind);
  PetscInt getGlobal(Ind3D ind);
  PetscInt getGlobal(IndPerp ind);

protected:
  GlobalIndexer(Mesh* localmesh);
  Mesh* fieldmesh;

private:
  /// This gets called by initialiseTest and is used to register
  /// fields with fake parallel meshes.
  virtual void registerFieldForTest(FieldData& f);
  virtual void registerFieldForTest(FieldPerp& f);

  PetscLib lib;

  /// Fields containing the indices for each element (as reals)
  Field3D indices3D;
  Field2D indices2D;
  FieldPerp indicesPerp;
  
  /// The only instance of this class acting on the global Mesh
  static IndexerPtr globalInstance;
  static bool initialisedGlobal;
  static Mesh* globalmesh;
  bool initialised;
};


/*!
 * A class which is used to assign to a particular element of a PETSc
 * vector. It is meant to be transient and will be destroyed immediately
 * after use. In general you should not try to assign an instance to a
 * variable.
 */
class PetscVectorElement {
public:
  BoutReal operator=(BoutReal val);
  BoutReal operator+=(BoutReal val);

  /// This is the only valid method for constructing new instances,
  /// guaranteeing they can be safely deleted once used.
  static PetscVectorElement& newElement(Vec* vector, PetscInt index);

  // Prevent non-transient copies from being created
  PetscVectorElement() = delete;
  PetscVectorElement(const PetscVectorElement& p) = delete;
  PetscVectorElement& operator=(const PetscVectorElement& rhs) = delete;
  PetscVectorElement& operator=(PetscVectorElement&& rhs) = delete;

private:
  PetscVectorElement(Vec* vector, int index);
  Vec* petscVector;
  PetscInt petscIndex;
};


/*!
 * A class which is used to assign to a particular element of a PETSc
 * matrix, potentially with a y-offset. It is meant to be transient
 * and will be destroyed immediately after use. In general you should
 * not try to assign an instance to a variable.
 */
class PetscMatrixElement {
public:
  BoutReal operator=(BoutReal val);
  BoutReal operator+=(BoutReal val);

  /// This is the only valid method for constructing new instances,
  /// guaranteeing they can be safely deleted once used.
  static PetscMatrixElement& newElement(Mat* matrix, PetscInt row, PetscInt col,
					std::vector<PetscInt> p = std::vector<PetscInt>(),
					std::vector<BoutReal> w = std::vector<BoutReal>());

  // Prevent non-transient copies from being created
  PetscMatrixElement() = delete;
  PetscMatrixElement(const PetscMatrixElement& p) = delete;
  PetscMatrixElement& operator=(const PetscMatrixElement& rhs) = delete;
  PetscMatrixElement& operator=(PetscMatrixElement&& rhs) = delete;

private:
  PetscMatrixElement(Mat* matrix, PetscInt row,
		     std::vector<PetscInt> p, std::vector<BoutReal> w);
  void setValues(BoutReal val, InsertMode mode);
  Mat* petscMatrix;
  PetscInt petscRow;
  std::vector<PetscInt> positions;
  std::vector<BoutReal> weights;
};


/*!
 * A class which wraps PETSc vector objects, allowing them to be
 * indexed using the BOUT++ scheme.
 */
template<class F> class PetscVector;
template<class F> void swap(PetscVector<F>& first, PetscVector<F>& second);

template <class F>
class PetscVector {
public:
  using ind_type = typename F::ind_type;
  /// Default constructor does nothing
  PetscVector() {
    initialised = false;
  }
  
  /// Copy constructor
  PetscVector(const PetscVector<F>& v) {
    VecCopy(v.vector, vector);
    indexConverter = v.indexConverter;
    location = v.location;
    initialised = v.initialised;
  }
  /// Move constrcutor
  PetscVector(PetscVector<F>&& v) {
    vector = v.vector;
    indexConverter = v.indexConverter;
    location = v.location;
    initialised = v.initialised;
  }

  /// Construct from a field, copying over the field values
  PetscVector(F f) {
    MPI_Comm comm;
    if (std::is_same<F, FieldPerp>::value) {
      comm = f.getMesh()->getXcomm();
    } else {
      comm = BoutComm::get();
    }
    indexConverter = GlobalIndexer::getInstance(f.getMesh());
    int size;
    if (std::is_same<F, FieldPerp>::value) {
      size = f.getMesh()->localSizePerp();
    } else if (std::is_same<F, Field2D>::value) {
      size = f.getMesh()->localSize2D();
    } else if (std::is_same<F, Field3D>::value) {
      size = f.getMesh()->localSize3D();
    } else {
      throw BoutException("PetscVector initialised for non-field type.");
    }
    VecCreateMPI(comm, size, PETSC_DECIDE, &vector);
    location = f.getLocation();
    initialised = true;
    *this = f;
  }
  
  ~PetscVector() {
    // FIXME: Should I add a check to ensure the vector has actually been created in the first place? Is that possible in Petsc?
    VecDestroy(&vector);
  }
  
  /// Copy assignment
  PetscVector<F>& operator=(const PetscVector<F>& rhs) {
    swap(*this, rhs);
    return *this;
  }
  
  /// Move assignment
  PetscVector<F>& operator=(PetscVector<F>&& rhs) {
    vector = rhs.vector;
    indexConverter = rhs.indexConverter;
    location = rhs.location;
    initialised = rhs.initialised;
  }
  
  friend void swap<F>(PetscVector<F>& first, PetscVector<F>& second);

  /// Assign from field, copying over field values. The vector must
  /// already have been created for this to work.
  PetscVector<F>& operator=(F& rhs) {
    PetscInt ind;
    ASSERT1(initialised);
    ASSERT2(location == rhs.getLocation());
    // Note that physical boundaries will not be set here. Not sure if that matters.
    BOUT_FOR(i, rhs.getRegion(RGN_NOBNDRY)) {
      ind = indexConverter->getGlobal(i);
      VecSetValues(vector, 1, &ind, &rhs[i], INSERT_VALUES);
    }
    assemble();
    return *this;
  }

  PetscVectorElement& operator()(ind_type& index) {
    return PetscVectorElement::newElement(&vector, indexConverter->getGlobal(index));
  }
  
  void assemble() {
    VecAssemblyBegin(vector);
    VecAssemblyEnd(vector);
  }
  
  void destroy() {
    VecDestroy(&vector);
  }

  /// Returns a field constructed from the contents of this vector
  const F toField() {
    F result;
    result.allocate();
    result.setLocation(location);
    PetscScalar val;
    PetscInt ind;
    // Need to set the physical boundaries, I think
    BOUT_FOR(i, result.getRegion(RGN_NOBNDRY)) {
      ind = indexConverter->getGlobal(i);
      VecGetValues(vector, 1, &ind, &val);
      result[i] = val;
    }
    return result;
  }

  /// Provides a reference to the raw PETSc Vec object.
  Vec* getVectorPointer() {
    return &vector;
  }

private:
  Vec vector;
  IndexerPtr indexConverter;
  CELL_LOC location;
  bool initialised;
  PetscLib lib;
};


/*!
 * A class which wraps PETSc vector objects, allowing them to be
 * indexed using the BOUT++ scheme. It provides the option of setting
 * a y-offset that interpolates onto field lines.
 */
template<class F> class PetscMatrix;
template<class F> void swap(PetscMatrix<F>& first, PetscMatrix<F>& second);

template <class F>
class PetscMatrix {
public:
  using ind_type = typename F::ind_type;
  /// Default constructor does nothing
  PetscMatrix() {  }

  /// Copy constructor
  PetscMatrix(const PetscMatrix<F>& v) : pt(v.pt) {
    throw BoutException("Not implemented");
  }
  /// Move constrcutor
  PetscMatrix(PetscMatrix<F>&& v) : pt(v.pt) {
    throw BoutException("Not implemented");
  }

  PetscMatrix(F &f) {
    throw BoutException("Not implemented");    
  }
  
  ~PetscMatrix() {
    throw BoutException("Not implemented");
  }

  /// Copy assignment
  PetscMatrix<F>& operator=(const PetscMatrix<F>& rhs) {
    swap(*this, rhs);
    return *this;
  }
  /// Move assignment
  PetscMatrix<F>& operator=(PetscMatrix<F>&& rhs);
  friend void swap<F>(PetscMatrix<F>& first, PetscMatrix<F>& second);

  PetscMatrixElement& operator()(ind_type& index1, ind_type& index2) {
    throw BoutException("Not implemented");
  }
  void assemble() {
    throw BoutException("Not implemented");
  }
  void destroy() {
    throw BoutException("Not implemented");
  }

  PetscMatrix<F>& yup(int index = 0) {
    return ynext(index + 1);
  }
  PetscMatrix<F>& ydown(int index = 0) {
    return ynext(-index - 1);
  }
  PetscMatrix<F>& ynext(int dir) {
    throw BoutException("Not implemented");
  }

  /// Provides a reference to the raw PETSc Mat object.
  Mat* getMatrixPointer() {
    return matrix.get();
  }

private:
  shared_ptr<Mat> matrix;
  ParallelTransform* pt;
  int yoffset;
  PetscLib lib;
};


/*!
 * Move reference to one Vec from \p first to \p second and
 * vice versa.
 */
template <class F>
void swap(PetscVector<F>& first, PetscVector<F>& second) {
  VecSwap(first.vector, second.vector);
  swap(first.indexConverter, second.indexConverter);
  swap(first.location, second.location);
  swap(first.initialised, second.initialised);
}

/*!
 * Move reference to one Mat from \p first to \p second and
 * vice versa.
 */
template <class F>
void swap(PetscMatrix<F>& first, PetscMatrix<F>& second) {
    throw BoutException("Not implemented");
}

/*!
 *  Performs matrix-multiplication on the supplied vector
 */
template <class F>
PetscVector<F> operator*(const PetscMatrix<F>& mat, PetscVector<F>& vec) {
  throw BoutException("Not implemented");
}
  

#endif // BOUT_HAS_PETSC

#endif // __PETSC_INTERFACE_H__
