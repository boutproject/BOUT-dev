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

#include <bout_types.hxx>
#include <bout/petsclib.hxx>
#include <bout/region.hxx>
#include <bout/mesh.hxx>
#include <bout/paralleltransform.hxx>

using namespace std;

#ifdef BOUT_HAS_PETSC
class GlobalIndexer;
using IndexerPtr = shared_ptr<GlobalIndexer>;
using InterpolationWeights = std::vector<ParallelTransform::positionsAndWeights>;

enum PetscObjectMode {unset_mode, insert_mode, add_mode};

template<class I> I makeGlobalIndexer(Mesh* localmesh);
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

private:
  /// This gets called by initialiseTest and is used to register
  /// fields with fake parallel meshes.
  virtual void registerFieldsForTest(Mesh* localmesh, FieldData& f);
  virtual void registerFieldsForTest(Mesh* localmesh, FieldPerp& f);

  /// The only instance of this class acting on the global Mesh
  static IndexerPtr globalInstance;
  Mesh* fieldmesh;
  bool initialised;
};


/*!
 * An abstract class used when assigning to an element of a PETSc
 * vector or matrix.
 */
class PetscElement {
public:
  virtual BoutReal operator=(BoutReal val);
  virtual BoutReal operator+=(BoutReal val);
};


/*!
 * A class which is used to assign to a particular element of a PETSc
 * vector. It is meant to be transient and will be destroyed immediately
 * after use. In general you should not try to assign an instance to a
 * variable.
 */
class PetscVectorElement : public PetscElement {
public:
  BoutReal operator=(BoutReal val) override;
  BoutReal operator+=(BoutReal val) override;

  /// This is the only valid method for constructing new instances,
  /// guaranteeing they can be safely deleted once used.
  static PetscVectorElement& newElement(Vec* vector, PetscInt index);

  // Prevent non-transient copies from being created
  PetscVectorElement() = delete;
  PetscVectorElement(const PetscVectorElement& p) = delete;
  PetscVectorElement& operator=(const PetscVectorElement& rhs) = delete;
  PetscVectorElement& operator=(PetscVectorElement&& rhs) = delete;

private:
  PetscVectorElement(Vec* vector, int index) : petscVector(vector),
					       petscIndex(index) {};
  Vec* petscVector;
  PetscInt petscIndex;
};


/*!
 * A class which is used to assign to a particular element of a PETSc
 * matrix, potentially with a y-offset. It is meant to be transient
 * and will be destroyed immediately after use. In general you should
 * not try to assign an instance to a variable.
 */
class PetscMatrixElement : public PetscElement {
public:
  BoutReal operator=(BoutReal val) override;
  BoutReal operator+=(BoutReal val) override;

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
  PetscMatrixElement(Mat* matrix, PetscInt row, PetscInt col,
		     std::vector<PetscInt> p,
		     std::vector<BoutReal> w)
    : petscMatrix(matrix), petscRow(row), petscCol(col), positions(p),
      weights(w) {};
  Mat* petscMatrix;
  PetscInt petscRow, petscCol;
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
  PetscVector() {}
  
  /// Copy constructor
  PetscVector(const PetscVector<F>& v) {
    throw BoutException("Not implemented");
  }
  /// Move constrcutor
  PetscVector(PetscVector<F>&& v) {
    throw BoutException("Not implemented");
  }

  /// Construct from a field, copying over the field values
  PetscVector(F f) {
    throw BoutException("Not implemented");    
  }
  
  ~PetscVector() {
    throw BoutException("Not implemented");
  }
  
  /// Copy assignment
  PetscVector<F>& operator=(const PetscVector<F>& rhs) {
    throw BoutException("Not implemented");
  }
  /// Move assignment
  PetscVector<F>& operator=(PetscVector<F>&& rhs) {
    throw BoutException("Not implemented");
  }
  friend void swap<F>(PetscVector<F>& first, PetscVector<F>& second);
  /// Assign from field, copying over field values
  PetscVector<F>& operator=(F& rhs) {
    throw BoutException("Not Implemented");
  }

  PetscVectorElement& operator()(ind_type& index) {
    throw BoutException("Not implemented");
  }
  void assemble() {
    throw BoutException("Not implemented");
  }
  void destroy() {
    throw BoutException("Not implemented");
  }

  /// Returns a field constructed from the contents of thsi vector
  const F toField() {
    throw BoutException("Not implemented");
  }

  /// Provides a reference to the raw PETSc Vec object.
  Vec* getVectorPointer() {
    return &vector;
  }

private:
  Vec vector;
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
  PetscMatrix() {}

  /// Copy constructor
  PetscMatrix(const PetscMatrix<F>& v) : pt(v.pt) {
    throw BoutException("Not implemented");
  }
  /// Move constrcutor
  PetscMatrix(PetscMatrix<F>&& v) : pt(v.pt) {
    throw BoutException("Not implemented");
  }

  PetscMatrix(F &f);
  
  ~PetscMatrix() {
    throw BoutException("Not implemented");
  }

  /// Copy assignment
  PetscMatrix<F>& operator=(const PetscMatrix<F>& rhs) {
    throw BoutException("Not implemented");
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
    throw BoutException("Not implemented");
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
