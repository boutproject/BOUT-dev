/*
 * Handle arrays of data
 *
 * Provides an interface to create, iterate over and release
 * arrays of templated types.
 *
 * Each type and array size has an object store, so when
 * arrays are released they are put into a store. Rather
 * than allocating memory, objects are retrieved from the
 * store. This minimises new and delete operations.
 *
 *
 * Ben Dudson, University of York, 2015
 *
 *
 * Changelog
 * ---------
 *
 * 2015-03-04  Ben Dudson <bd512@york.ac.uk>
 *     o Initial version
 *
 * 2021 Holger Jones, Ben Dudson
 *     o Added Umpire support, in multiple iterations/variations
 */

#ifndef BOUT_ARRAY_H
#define BOUT_ARRAY_H

#include <algorithm>
#include <map>
#include <memory>
#include <vector>

#if BOUT_USE_OPENMP
#include <omp.h>
#endif

#include "bout/build_defines.hxx"

#if BOUT_HAS_UMPIRE
#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"
#endif

#include <bout/assert.hxx>
#include <bout/openmpwrap.hxx>

namespace {
template <typename T>
using iterator = T*;
template <typename T>
using const_iterator = const T*;
} // namespace

/*!
 * ArrayData holds the actual data
 * Handles the allocation and deletion of data
 */
template <typename T>
struct ArrayData {
  ArrayData(int size) : len(size) {
    // Allocate memory array
    // Note: By allocating with Umpire, the raw data
    //       array can be accessed on GPU devices
    //       even though the Array object itself can't.
#if BOUT_HAS_UMPIRE
    auto& rm = umpire::ResourceManager::getInstance();
#if BOUT_HAS_CUDA
    auto allocator = rm.getAllocator(umpire::resource::Pinned);
#else
    auto allocator = rm.getAllocator("HOST");
#endif
    data = static_cast<T*>(allocator.allocate(size * sizeof(T)));
#else // BOUT_HAS_UMPIRE
    data = new T[len];
#endif
  }

  /// Move constructor
  ArrayData(ArrayData&& in) noexcept : len(in.len), data(in.data) {
    in.len = 0;
    in.data = nullptr;
  }

  /// Disable copy, since ArrayData takes ownership of data
  ArrayData(const ArrayData& in) = delete;

  ~ArrayData() {
#if BOUT_HAS_UMPIRE
    auto& rm = umpire::ResourceManager::getInstance();
    rm.deallocate(data);
#else
    delete[] data;
#endif
  }
  iterator<T> begin() const { return data; }
  iterator<T> end() const { return data + len; }
  int size() const { return len; }

  /// Copy assignment
  /// Copy the underlying data from one array to the other
  ///
  /// @param in  ArrayData with the same len as this.
  ArrayData<T>& operator=(const ArrayData<T>& in) {
    ASSERT1(in.len == len); // Ensure that the array lengths are the same
    std::copy(std::begin(in), std::end(in), begin());
    return *this;
  }

  /// Move assignment
  ArrayData<T>& operator=(ArrayData<T>&& in) noexcept {
    if (this != &in) {
      // Free resources
#if BOUT_HAS_UMPIRE
      auto& rm = umpire::ResourceManager::getInstance();
      rm.deallocate(data);
#else
      delete[] data;
#endif
      // Copy pointers
      len = in.len;
      data = in.data;

      // Remove pointer from input so that it is
      // not freed multiple times
      in.len = 0;
      in.data = nullptr;
    }
    return *this;
  }

  inline T& operator[](int ind) { return data[ind]; }
  inline const T& operator[](int ind) const { return data[ind]; }

private:
  int len; ///< Size of the array
  T* data; ///< Array of data
};

/*!
 * Data array type with automatic memory management
 *
 * This implements a container similar to std::vector
 * but with reference counting like a smart pointer
 * and custom memory management to minimise new and delete calls
 *
 * This can be used as an alternative to static arrays
 *
 * Array<dcomplex> vals(100); // 100 complex numbers
 *
 * vals[10] = 1.0;  // ok
 *
 * When an Array goes out of scope or is deleted,
 * the underlying memory (dataBlock/Backing) is put into
 * a map, rather than being freed.
 * If the same size arrays are used repeatedly then this
 * avoids the need to use new and delete.
 *
 * This behaviour can be disabled by calling the static function useStore:
 *
 * Array<dcomplex>::useStore(false); // Disables memory store
 *
 * The second template argument determines what type of container to use to
 * store data. This defaults to a custom struct but can be std::valarray (
 * provided T is a compatible type), std::vector etc. Must provide the following :
 *  size, operator=, operator[], begin, end
 *
 * Notes:
 *  - Arrays can't be used in GPU code. To access Array data
 *    inside a RAJA loop, first extract the raw pointer
 */
template <typename T, typename Backing = ArrayData<T>>
class Array {
public:
  using data_type = T;
  using backing_type = Backing;
  using size_type = int;

  /*!
   * Create an empty array
   *
   * Array a();
   * a.empty(); // True
   *
   */
  Array() noexcept : ptr(nullptr) {}

  /*!
   * Create an array of given length
   */
  Array(size_type len) { ptr = get(len); }

  /*!
   * Destructor. Releases the underlying dataBlock
   */
  ~Array() noexcept { release(ptr); }

  /*!
   * Copy constructor
   */
  Array(const Array& other) noexcept { ptr = other.ptr; }

  /*!
   * Assignment operator
   * After this both Arrays share the same dataBlock
   *
   * Uses copy-and-swap idiom
   */
  Array& operator=(Array other) noexcept {
    swap(*this, other);
    return *this;
  }

  /*!
   * Move constructor
   */
  Array(Array&& other) noexcept { swap(*this, other); }

  /*!
   * Reallocate the array with size = \p new_size
   *
   * Note that this invalidates the existing data!
   */
  void reallocate(size_type new_size) {
    release(ptr);
    ptr = get(new_size);
  }

  /*!
   * Holds a static variable which controls whether
   * memory blocks (dataBlock) are put into a store
   * or new/deleted each time.
   *
   * The variable is initialised to true on first use,
   * but can be set to false by passing "false" as input.
   * Once set to false it can't be changed back to true.
   */
  static bool useStore(bool keep_using = true) noexcept {
    static bool value = true;
    if (keep_using) {
      return value;
    }
    // Change to false
    value = false;
    return value;
  }

  /*!
   * Release data. After this the Array is empty and any data access
   * will be invalid
   */
  void clear() noexcept { release(ptr); }

  /*!
   * Delete all data from the store and disable the store
   *
   * Note: After this is called the store cannot be re-enabled
   */
  static void cleanup() {
    // Clean the store, deleting data
    store(true);
    // Don't use the store anymore. This is so that array releases
    // after cleanup() get deleted rather than put into the store
    useStore(false);
  }

  /*!
   * Returns true if the Array is empty
   */
  bool empty() const noexcept { return ptr == nullptr; }

  /*!
   * Return size of the array. Zero if the array is empty.
   */
  size_type size() const noexcept {
    if (!ptr) {
      return 0;
    }

    // Note: std::valarray::size is technically not noexcept, so
    // Array::size shouldn't be either if we're using valarrays -- in
    // practice, it is so this shouldn't matter
    return ptr->size();
  }

  /*!
   * Returns true if the data is unique to this Array.
   *
   */
  bool unique() const noexcept { return ptr.use_count() == 1; }

  /*!
   * Ensures that this Array does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() {
    if (!ptr || unique()) {
      return;
    }

    // Get a new (unique) block of data
    dataPtrType p = get(size());

    // Make copy of the underlying data
    p->operator=((*ptr));

    // Update the local pointer and release old
    // Can't just do ptr=p as need to try to add to store.
    release(ptr);
    ptr = std::move(p);
  }

  //////////////////////////////////////////////////////////
  // Iterators

  iterator<T> begin() noexcept { return (ptr) ? std::begin(*ptr) : nullptr; }

  iterator<T> end() noexcept { return (ptr) ? std::end(*ptr) : nullptr; }

  // Const iterators
  const_iterator<T> begin() const noexcept { return (ptr) ? std::begin(*ptr) : nullptr; }

  const_iterator<T> end() const noexcept { return (ptr) ? std::end(*ptr) : nullptr; }

  //////////////////////////////////////////////////////////
  // Element access

  /*!
   * Access a data element. This will fail if the Array is empty (so ptr is null),
   * or if ind is out of bounds. For efficiency no checking is performed,
   * so the user should perform checks.
   */
  inline T& operator[](size_type ind) {
    ASSERT3(0 <= ind && ind < size());
    return ptr->operator[](ind);
  }
  inline const T& operator[](size_type ind) const {
    ASSERT3(0 <= ind && ind < size());
    return ptr->operator[](ind);
  }

  /*!
   * Exchange contents with another Array of the same type.
   * Sizes of the arrays may differ.
   */
  friend void swap(Array<T, Backing>& first, Array<T, Backing>& second) noexcept {
    using std::swap;
    swap(first.ptr, second.ptr);
  }

private:
  // Type defs to help keep things brief -- which backing do we use
  using dataBlock = Backing;
  using dataPtrType = std::shared_ptr<dataBlock>;

  /*!
   * Pointer to the data container object owned by this Array.
   * May be null
   */
  dataPtrType ptr;

  using storeType = std::map<size_type, std::vector<dataPtrType>>;
  using arenaType = std::vector<storeType>;

  /*!
   * This maps from array size (size_type) to vectors of pointers to dataBlock objects
   *
   * By putting the static store inside a function it is initialised on first use,
   * and doesn't need to be separately declared for each type T
   *
   * Inputs
   * ------
   *
   * @param[in] cleanup   If set to true, deletes all dataBlock and clears the store
   */
  static storeType& store(bool cleanup = false) {
    static arenaType arena(omp_get_max_threads());
    if (!cleanup) {
      return arena[omp_get_thread_num()];
    }

    // Clean by deleting all data -- possible that just stores.clear() is
    // sufficient rather than looping over each entry.
    BOUT_OMP_SAFE(single)
    {
      for (auto& stores : arena) {
        for (auto& p : stores) {
          auto& v = p.second;
          for (dataPtrType a : v) {
            a.reset();
          }
          v.clear();
        }
        stores.clear();
      }
      // Here we ensure there is exactly one empty map still
      // left in the arena as we have to return one such item
      arena.resize(1);
    }

    // Store should now be empty but we need to return something,
    // so return an empty storeType from the arena.
    return arena[0];
  }

  /*!
   * Returns a pointer to a dataBlock object of size \p len with no
   * references. This is either from the store, or newly allocated
   *
   * Expects \p len >= 0
   */
  dataPtrType get(size_type len) {
    ASSERT3(len >= 0);

    dataPtrType p;

    auto& st = store()[len];

    if (!st.empty()) {
      p = st.back();
      st.pop_back();
    } else {
      // Ensure that when we release the data block later we'll have
      // enough space to put it in the store so that `release` can be
      // noexcept
      st.reserve(1);
      p = std::make_shared<dataBlock>(len);
    }

    return p;
  }

  /*!
   * Release an dataBlock object, reducing its reference count by one.
   * If no more references, then put back into the store.
   * It's important to pass a reference to the pointer, otherwise we get
   * a copy of the shared_ptr, which therefore increases the use count
   * and doesn't allow us to free the pass pointer directly
   *
   * Note that this is noexcept only because we've ensure that both a)
   * store()[<size>] already exists, and b) it has space for at least
   * one data block. Of course, store() could throw -- in which case
   * we're doomed anyway, so the only thing we can do is abort
   */
  void release(dataPtrType& d) noexcept {
    if (!d) {
      return;
    }

    // Reduce reference count, and if zero return to store
    if (d.use_count() == 1) {
      if (useStore()) {
        // Put back into store
        store()[d->size()].push_back(std::move(d));
        // Could return here but seems to slow things down a lot
      }
    }

    // Finish by setting pointer to nullptr if not putting on store
    d = nullptr;
  }
};

/*!
 * Create a copy of an Array, which does not share data
 */
template <typename T, typename Backing>
Array<T, Backing> copy(const Array<T, Backing>& other) {
  Array<T, Backing> a(other);
  a.ensureUnique();
  return a;
}

/// Compare arrays, which may have different backing
template <typename T, typename B1, typename B2>
bool operator==(const Array<T, B1>& lhs, const Array<T, B2>& rhs) {
  return std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

#endif // BOUT_ARRAY_H
