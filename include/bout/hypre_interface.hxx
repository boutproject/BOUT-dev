#ifndef BOUT_HYPRE_INTERFACE_H
#define BOUT_HYPRE_INTERFACE_H

#ifdef BOUT_HAS_HYPRE

#include "bout/globalindexer.hxx"
#include "field.hxx"
#include "utils.hxx"
#include "boutcomm.hxx"

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"

namespace bout {

template <class T>
class HypreVector;
template <class T>
void swap(HypreVector<T>& lhs, HypreVector<T>& rhs);

template <class T>
class HypreVector {
  HYPRE_IJVector hypre_vector{nullptr};
  HYPRE_ParVector parallel_vector{nullptr};
  IndexerPtr indexConverter;
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

  explicit HypreVector(const T& f) {
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? f.getMesh()->getXcomm() : BoutComm::get();
    indexConverter = GlobalIndexer::getInstance(f.getMesh());

    const auto region = f.getRegion("RGN_ALL_THIN");
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
    assemble();
  }

  HypreVector<T>& operator=(const T& f) {
    HypreVector<T> result(f);
    *this = std::move(result);
    return *this;
  }

  // Mesh can't be const yet, need const-correctness on Mesh;
  // GlobalIndexer ctor also modifies mesh -- FIXME
  explicit HypreVector(Mesh& mesh) {
    const MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? mesh.getXcomm() : BoutComm::get();
    indexConverter = GlobalIndexer::getInstance(&mesh);

    const auto region = mesh.getRegion<T>("RGN_ALL_THIN");

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
  const HYPRE_IJVector get() const { return hypre_vector; }

  class Element {
    HYPRE_IJVector vector{nullptr};
    HYPRE_BigInt index;
    HYPRE_Complex value;

  public:
    Element() = delete;
    Element(const Element&) = default;
    Element(Element&&) = default;
    Element(HYPRE_IJVector vector_, int index_) : vector(vector_), index(index_) {
      int status;
      BOUT_OMP(critical) {
        if((status = HYPRE_IJVectorGetValues(vector, 1, &index, &value)) != 0) {
          throw BoutException("Error when getting elements of a HYPRE vector, index = {}", index);
        }
      }
    }

    Element& operator=(Element& other) { return *this = static_cast<BoutReal>(other); }
    Element& operator=(BoutReal value_) {
      value = value_;
      int status;
      BOUT_OMP(critical) {
        if((status = HYPRE_IJVectorSetValues(vector, 1, &index, &value)) != 0) {
          throw BoutException("Error when setting elements of a HYPRE vector, index = {}", index);
        }
      }
      return *this;
    }
    Element& operator+=(BoutReal value_) {
      value = value_;
      int status;
      BOUT_OMP(critical) {
        if((status = HYPRE_IJVectorAddToValues(vector, 1, &index, &value)) != 0) {
          throw BoutException("Error when setting elements of a HYPRE vector, index = {}", index);
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

} // namespace bout

#endif // BOUT_HAS_HYPRE

#endif // BOUT_HYPRE_INTERFACE_H
