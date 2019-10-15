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
#include <algorithm>

#include <bout_types.hxx>
#include <bout/petsclib.hxx>
#include <bout/region.hxx>
#include <bout/mesh.hxx>
#include <boutcomm.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/traits.hxx>

#ifdef BOUT_HAS_PETSC
class GlobalIndexer;
using IndexerPtr = std::shared_ptr<GlobalIndexer>;
using InterpolationWeights = std::vector<ParallelTransform::PositionsAndWeights>;

/*!
 * A singleton which accepts index objects produced by iterating over
 * fields and returns a global index. This index can be used when
 * constructing PETSc arrays. Guard regions used for communication
 * between processes will have the indices of the part of the interior
 * region they are mirroring. Boundaries will have unique indices, but
 * are only included to a depth of 1.
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
  Mesh* getMesh();
  
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
 * A class which wraps PETSc vector objects, allowing them to be
 * indexed using the BOUT++ scheme. Note that boundaries are only
 * included in the vector to a depth of 1.
 */
template<class T> class PetscVector;
template<class T> void swap(PetscVector<T>& first, PetscVector<T>& second);

template <class T>
class PetscVector {
public:
  static_assert(bout::utils::is_Field<T>::value, "PetscVector only works with Fields");
  using ind_type = typename T::ind_type;
  /// Default constructor does nothing
  PetscVector() = default;
  
  /// Copy constructor
  PetscVector(const PetscVector<T>& v) {
    VecDuplicate(v.vector, &vector);
    VecCopy(v.vector, vector);
    indexConverter = v.indexConverter;
    location = v.location;
    initialised = v.initialised;
  }

  /// Move constrcutor
  PetscVector(PetscVector<T>&& v) {
    vector = v.vector;
    indexConverter = v.indexConverter;
    location = v.location;
    initialised = v.initialised;
    v.vector = nullptr;
    v.initialised = false;
  }

  /// Construct from a field, copying over the field values
  PetscVector(const T &f) {
    MPI_Comm comm;
    if (std::is_same<T, FieldPerp>::value) {
      comm = f.getMesh()->getXcomm();
    } else {
      comm = BoutComm::get();
    }
    indexConverter = GlobalIndexer::getInstance(f.getMesh());
    int size;
    if (std::is_same<T, FieldPerp>::value) {
      size = f.getMesh()->localSizePerp();
    } else if (std::is_same<T, Field2D>::value) {
      size = f.getMesh()->localSize2D();
    } else if (std::is_same<T, Field3D>::value) {
      size = f.getMesh()->localSize3D();
    } else {
      throw BoutException("PetscVector initialised for non-field type.");
    }
    VecCreateMPI(comm, size, PETSC_DECIDE, &vector);
    location = f.getLocation();
    initialised = true;
    PetscInt ind;
    BOUT_FOR(i, f.getRegion(RGN_ALL)) {
      ind = indexConverter->getGlobal(i);
      if (ind != -1) {
	VecSetValues(vector, 1, &ind, &f[i], INSERT_VALUES);
      }
    }
    assemble();
  }

  /// Construct a vector like v, but using data from a raw PETSc
  /// Vec. That Vec (not a copy) will then be owned by the new object.
  PetscVector(const PetscVector<T> v, Vec vec) {
#if CHECKLEVEL >= 2
    int fsize, msize;
    if (std::is_same<T, FieldPerp>::value) {
      fsize = v.indexConverter->getMesh()->localSizePerp();
    } else if (std::is_same<T, Field2D>::value) {
      fsize = v.indexConverter->getMesh()->localSize2D();
    } else if (std::is_same<T, Field3D>::value) {
      fsize = v.indexConverter->getMesh()->localSize3D();
    } else {
      throw BoutException("PetscVector initialised for non-field type.");
    }
    VecGetSize(vec, &msize);
    ASSERT2(fsize == msize);
#endif
    vector = vec;
    indexConverter = v.indexConverter;
    location = v.location;
    initialised = true;
}

  ~PetscVector() {
    // FIXME: Should I add a check to ensure the vector has actually been created in the first place? Is that possible in Petsc?
    if (vector != nullptr && initialised) {
      VecDestroy(&vector);
    }
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

  friend void swap<T>(PetscVector<T>& first, PetscVector<T>& second);

  /*!
   * A class which is used to assign to a particular element of a PETSc
   * vector. It is meant to be transient and will be destroyed immediately
   * after use. In general you should not try to assign an instance to a
   * variable.
   */
  class Element {
  public:
    Element() = delete;
    Element(Vec* vector, int index) : petscVector(vector), petscIndex(index) {}
    Element operator=(BoutReal val) {
      auto status = VecSetValues(*petscVector, 1, &petscIndex, &val, INSERT_VALUES);
      if (status != 0) {
	throw BoutException("Error when setting elements of a PETSc vector.");
      }
      return *this;
    }
    Element operator+=(BoutReal val) {
      auto status = VecSetValues(*petscVector, 1, &petscIndex, &val, ADD_VALUES);
      if (status != 0) {
	throw BoutException("Error when setting elements of a PETSc vector.");
      }
      return *this;
    }
    operator BoutReal() const {
      BoutReal value;
      auto status = VecGetValues(*petscVector, 1, &petscIndex, &value);
      if (status != 0) {
	throw BoutException("Error when getting elements of a PETSc vector.");
      }  
      return value;
    }

  private:
    Vec* petscVector;
    PetscInt petscIndex;
  };
  
  Element operator()(ind_type& index) {
#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Can not return element of uninitialised vector");
    } 
#endif
    int global = indexConverter->getGlobal(index);
#if CHECKLEVEL >= 1
    if (global == -1) {
      throw BoutException("Request to return invalid vector element");
    }
#endif
    return Element(&vector, global);
  }
  
  void assemble() {
    VecAssemblyBegin(vector);
    VecAssemblyEnd(vector);
  }
  
  void destroy() {
    if (vector != nullptr && initialised) {
      VecDestroy(&vector);
      vector = nullptr;
      initialised = false;
    }
  }

  /// Returns a field constructed from the contents of this vector
  T toField() {
    T result(indexConverter->getMesh());
    result.allocate();
    result.setLocation(location);
    PetscScalar val;
    PetscInt ind;
    // Note that this only populates boundaries to a depth of 1
    BOUT_FOR(i, result.getRegion(RGN_ALL)) {
      ind = indexConverter->getGlobal(i);
      if (ind == -1) {
	result[i] = -1.0;
      } else {
	VecGetValues(vector, 1, &ind, &val);
	result[i] = val;
      }
    }
    return result;
  }

  /// Provides a reference to the raw PETSc Vec object.
  Vec* getVectorPointer() {
    return &vector;
  }

private:
  Vec vector = nullptr;
  IndexerPtr indexConverter;
  CELL_LOC location;
  bool initialised = false;
  PetscLib lib;
};


/*!
 * A class which wraps PETSc vector objects, allowing them to be
 * indexed using the BOUT++ scheme. It provides the option of setting
 * a y-offset that interpolates onto field lines.
 */
template<class T> class PetscMatrix;
template<class T> void swap(PetscMatrix<T>& first, PetscMatrix<T>& second);

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
  
  /// Default constructor does nothing
  PetscMatrix() : matrix(new Mat(), MatrixDeleter()) {}

  /// Copy constructor
  PetscMatrix(const PetscMatrix<T>& m) : matrix(new Mat(), MatrixDeleter()), pt(m.pt) {
    MatDuplicate(*m.matrix, MAT_COPY_VALUES, matrix.get());
    // MatCopy(*m.matrix, *matrix);
    indexConverter = m.indexConverter;
    yoffset = m.yoffset;
    initialised = m.initialised;
  }
  
  /// Move constrcutor
  PetscMatrix(PetscMatrix<T>&& m) : pt(m.pt) {
    matrix = m.matrix;
    indexConverter = m.indexConverter;
    yoffset = m.yoffset;
    initialised = m.initialised;
    m.initialised = false;
  }

  // Construct a matrix capable of operating on the specified field
  PetscMatrix(T &f) : matrix(new Mat(), MatrixDeleter()) {
    MPI_Comm comm;
    if (std::is_same<T, FieldPerp>::value) {
      comm = f.getMesh()->getXcomm();
    } else {
      comm = BoutComm::get();
    }
    indexConverter = GlobalIndexer::getInstance(f.getMesh());
    pt = &f.getMesh()->getCoordinates()->getParallelTransform();
    int size;
    if (std::is_same<T, FieldPerp>::value) {
      size = f.getMesh()->localSizePerp();
    } else if (std::is_same<T, Field2D>::value) {
      size = f.getMesh()->localSize2D();
    } else if (std::is_same<T, Field3D>::value) {
      size = f.getMesh()->localSize3D();
    } else {
      throw BoutException("PetscVector initialised for non-field type.");
    }
    MatCreate(comm, matrix.get());
    MatSetSizes(*matrix, size, size, PETSC_DECIDE, PETSC_DECIDE);
    MatSetType(*matrix, MATMPIAIJ);
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
  PetscMatrix<T>& operator=(PetscMatrix<T>&& rhs) {
    matrix = rhs.matrix;
    indexConverter = rhs.indexConverter;
    pt = rhs.pt;
    yoffset = rhs.yoffset;
    initialised = rhs.initialised;
    rhs.initialised = false;
  }
  friend void swap<T>(PetscMatrix<T>& first, PetscMatrix<T>& second);

  /*!
   * A class which is used to assign to a particular element of a PETSc
   * matrix, potentially with a y-offset. It is meant to be transient
   * and will be destroyed immediately after use. In general you should
   * not try to assign an instance to a variable.
   */
  class Element {
  public:
    Element() = delete;
    Element(Mat* matrix, PetscInt row, PetscInt col, std::vector<PetscInt> p = {},
	    std::vector<BoutReal> w = {}) : petscMatrix(matrix),
					    petscRow(row), positions(p),
					    weights(w) {
      ASSERT2(positions.size() == weights.size());
      if (positions.size() == 0) {
	positions = { col };
	weights = { 1.0 };
      }
    }
    Element operator=(BoutReal val) {
      setValues(val, INSERT_VALUES);
      return *this;
    }
    Element operator+=(BoutReal val) {
      setValues(val, ADD_VALUES);
      return *this;
    }
    //operator BoutReal() const {
      // Tricky...
    //}
    
  private:
    void setValues(BoutReal val, InsertMode mode) {
      ASSERT3(positions.size() > 0);
      std::vector<PetscScalar> values;
      std::transform(weights.begin(), weights.end(), std::back_inserter(values),
		     [&val](BoutReal weight) -> PetscScalar {return weight * val;});
      auto status = MatSetValues(*petscMatrix, 1, &petscRow, positions.size(),
				 positions.data(), values.data(), mode);
      if (status != 0) {
	throw BoutException("Error when setting elements of a PETSc matrix.");
      }
    }
    Mat* petscMatrix;
    PetscInt petscRow;
    std::vector<PetscInt> positions;
    std::vector<BoutReal> weights;
  };
  
  Element operator()(ind_type& index1, ind_type& index2) {
    int global1 = indexConverter->getGlobal(index1),
        global2 = indexConverter->getGlobal(index2);
#if CHECKLEVEL >= 1
    if (!initialised) {
      throw BoutException("Can not return element of uninitialised matrix");
    } else if (global1 == -1 || global2 == -1) {
      throw BoutException("Request to return invalid matrix element");
    }
#endif
    std::vector<PetscInt> positions;
    std::vector<PetscScalar> weights;
    if (yoffset != 0) {
      const auto pw = [this, &index2]() {
	    if (this->yoffset == -1) {
	      return pt->getWeightsForYDownApproximation(index2.x(), index2.y(), index2.z());
	    } else if (this->yoffset == 1) {
	      return pt->getWeightsForYUpApproximation(index2.x(), index2.y(), index2.z());
	    } else {
	      return pt->getWeightsForYApproximation(index2.x(), index2.y(), index2.z(),
						     this->yoffset);
	    }
      }();
      
      int ny = indexConverter->getMesh()->LocalNy, nz = indexConverter->getMesh()->LocalNz;
      if (std::is_same<T, FieldPerp>::value) {
	ny = 1;
      } else if (std::is_same<T, Field2D>::value) {
	nz = 1;
      }
      std::transform(pw.begin(), pw.end(), std::back_inserter(positions),
		     [this, &ny, &nz](ParallelTransform::PositionsAndWeights p) -> PetscInt
		     {return this->indexConverter->getGlobal(ind_type(p.i*ny*nz + p.j*nz
								      + p.k, ny, nz));});
      std::transform(pw.begin(), pw.end(), std::back_inserter(weights),
		     [](ParallelTransform::PositionsAndWeights p) ->
		     PetscScalar {return p.weight;});
    }
    return Element(matrix.get(), global1, global2, positions, weights);
  }

  void assemble() {
    MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
  }

  void destroy() {
    if (*matrix != nullptr && initialised) {
      MatDestroy(matrix.get());
      *matrix = nullptr;
      initialised = false;
    }
  }

  PetscMatrix<T> yup(int index = 0) {
    return ynext(index + 1);
  }
  PetscMatrix<T> ydown(int index = 0) {
    return ynext(-index - 1);
  }
  PetscMatrix<T> ynext(int dir) {
    if (std::is_same<T, FieldPerp>::value && yoffset + dir != 0) {
      throw BoutException("Can not get ynext for FieldPerp");
    } else {
      PetscMatrix<T> result; // Can't use copy constructor because don't
                             // want to duplicate the matrix
      result.matrix = matrix;
      result.indexConverter = indexConverter;
      result.pt = pt;
      result.yoffset = std::is_same<T, Field2D>::value ? 0 : yoffset + dir;
      result.initialised = initialised;
      return result;
    }
  }

  /// Provides a reference to the raw PETSc Mat object.
  Mat* getMatrixPointer() {
    return matrix.get();
  }

private:
  std::shared_ptr<Mat> matrix = nullptr;
  IndexerPtr indexConverter;
  ParallelTransform* pt;
  int yoffset = 0;
  bool initialised = false;
  PetscLib lib;
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
PetscVector<T> operator*(PetscMatrix<T>& mat, PetscVector<T>& vec) {
  Vec rhs = *vec.getVectorPointer(), result;
  VecDuplicate(rhs, &result);
  VecAssemblyBegin(result);
  VecAssemblyEnd(result);  
  int err = MatMult(*mat.getMatrixPointer(), rhs, result);
  ASSERT2(err == 0);
  return PetscVector<T>(vec, result);
}
  

#endif // BOUT_HAS_PETSC

#endif // __PETSC_INTERFACE_H__
