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
 */


#ifndef __ARRAY_H__
#define __ARRAY_H__

#include <algorithm>
#include <map>
#include <memory>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bout/build_config.hxx"

#if BOUT_USE_CUDA && defined(__CUDACC__)
#define BOUT_HOST_DEVICE __host__ __device__
#define BOUT_HOST __host__
#define BOUT_DEVICE __device__
#else
#define BOUT_HOST_DEVICE
#define BOUT_HOST
#define BOUT_DEVICE
#endif

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
}

/*!
 * ArrayData holds the actual data
 * Handles the allocation and deletion of data
 */
template <typename T>
struct ArrayData {
  ArrayData(int size) : len(size) { data = new T[len]; }
  ~ArrayData() { delete[] data; }
  iterator<T> begin() const { return data; }
  iterator<T> end() const { return data + len; }
  int size() const { return len; }
  void operator=(ArrayData<T>& in) { std::copy(std::begin(in), std::end(in), begin()); }
  T& operator[](int ind) { return data[ind]; };
private:
  int len; ///< Size of the array
  T* data; ///< Array of data
};

/*!
 * SharedArrayData holds data and a reference count
 * Handles the allocation and deletion of data
 *
 * This exists so that Arrays can be used in CUDA code, which
 * requires that all accessed memory is allocated with Umpire.
 *
 * When used on the CPU it behaves like a shared pointer (pointing to itself)
 * When used on the GPU (__CUDA_ARCH__ defined) behaves like raw pointer
 *
 * An alternative might be to work out how to use an Umpire
 * allocator inside std::make_shared
 *
 * Note: In CUDA we can't create a new Array, but we can copy
 *       and delete those copies. This disables reference counting.
 * Note: __CUDA_ARCH__ is only defined device side
 */
template <typename T>
struct SharedArrayData {
  /// Default allocation.
  /// Note: The result contains a null data pointer
  SharedArrayData() {}

  /// Create a new array. This will call `new`, or use Umpire if enabled.
  /// Note: Both data and counter need to be allocated, preferably in the same data block
  SharedArrayData(int size) : len(size) {
#if BOUT_HAS_UMPIRE
    auto& rm = umpire::ResourceManager::getInstance();
#if BOUT_USE_CUDA
    auto allocator = rm.getAllocator(umpire::resource::Pinned);
#else
    auto allocator = rm.getAllocator("HOST");
#endif
    // Allocate space for both counter and data
    auto raw_ptr = allocator.allocate(size * sizeof(T) + sizeof(int));
    // data and counter are both pointers into the block of memory
    data = static_cast<T*>(raw_ptr);
    counter = static_cast<int*>(raw_ptr + size * sizeof(T));
#else
    // Un-managed
    data = new T[len];
    counter = new int;
#endif
    *counter = 1; // Initialise with a single reference
  }

  /// Delete and deallocate/delete the data
  /// If called from CUDA doesn't do anything, because in CUDA we don't take any ownership
  BOUT_HOST_DEVICE ~SharedArrayData() {
#ifndef __CUDA_ARCH__
    if (!data) {
      return; // No data allocated
    }
    (*counter)--; // Decrement counter. Not done in CUDA

    if (*counter == 0) {
#if BOUT_HAS_UMPIRE
      auto& rm = umpire::ResourceManager::getInstance();
      rm.deallocate(data); // counter is part of the same memory allocation
#else
      delete[] data;
      delete counter;
#endif
    }
#endif
  }

  /// Make a shallow copy, referring to the same data
  /// Increments the shared counter
  BOUT_HOST_DEVICE SharedArrayData(const SharedArrayData& other) {
    len = other.len;
    data = other.data;
    counter = other.counter;
#ifndef __CUDA_ARCH__
    (*counter)++; // Increment the shared reference counter
    // Note: In CUDA don't take any ownership, treat as raw pointer
#endif
  }

  /// Move. No need to update counter
  BOUT_HOST_DEVICE SharedArrayData(SharedArrayData&& other) {
    len = other.len;
    data = other.data;
    counter = other.counter;
  }

  BOUT_HOST_DEVICE inline iterator<T> begin() const { return data; }
  BOUT_HOST_DEVICE inline iterator<T> end() const { return data + len; }
  BOUT_HOST_DEVICE inline int size() const { return len; }

  /// Assign, deep copying the underlying data
  BOUT_HOST_DEVICE inline SharedArrayData<T>& operator=(const SharedArrayData<T>& in) {
#ifndef __CUDA_ARCH__
    std::copy(std::begin(in), std::end(in), begin());
#else
    int threadId =
        threadIdx.x + blockDim.x * threadIdx.y + (blockDim.x * blockDim.y) * threadIdx.z;

    if (threadId < len) {
      data[threadId] = in[threadId];
    }
#endif
    return *this;
  }

  BOUT_HOST_DEVICE inline T& operator[](int ind) { return data[ind]; };
  BOUT_HOST_DEVICE inline const T& operator[](int ind) const { return data[ind]; };

  /// Dereference operator.
  /// Here to emulate std::shared_ptr<ArrayData>
  SharedArrayData& operator*() {
    return *this;
  }

  /// Check if the data pointer is null
  explicit operator bool() const noexcept {
    return data != nullptr;
  }

  /// Return the number of references to this data
  /// Name chosen to emulate std::shared_ptr<ArrayData>
  int use_count() {
    if (!counter) {
      return 0;
    }
    return *counter;
  }

private:
  T* data = nullptr;
  int len = 0;
  int* counter = nullptr;
};

// Creation Policies
//
// These determine how backing arrays are created. The result of `create(int len)` should
// have shared pointer semantics
//

/// Create a shared pointer to the backing array
struct CreateSharedPointer {
  /// The type returned by create(len) is a pointer to the dataBlock type
  template<class dataBlock>
  using dataPtrType = std::shared_ptr<dataBlock>;

  template<class dataBlock>
  std::shared_ptr<dataBlock> create(int len) {
    return std::make_shared<dataBlock>(len);
  }
};

/// Special case where array also has pointer semantics (SharedArrayData)
struct CreateDirectly {
  template<class dataBlock>
  using dataPtrType = dataBlock;

  template<class dataBlock>
  dataBlock create(int len) {
    return dataBlock(len);
  }
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
 */

#if BOUT_HAS_UMPIRE
template <typename T, typename Backing = SharedArrayData<T>, typename CreationPolicy = CreateDirectly>
#else
template <typename T, typename Backing = ArrayData<T>, typename CreationPolicy = CreateSharedPointer>
#endif
class Array : private CreationPolicy {
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
  Array() noexcept : ptr() {}

  /*!
   * Create an array of given length
   */
  Array(size_type len) {
    ptr = get(len);
  }

  /*!
   * Destructor. Releases the underlying dataBlock
   */
  BOUT_HOST_DEVICE inline ~Array() {
#ifndef __CUDA_ARCH__
    // Only free resources if not running on the device
    // If running on the device don't take ownership, treat like raw pointer
    release(ptr);
#endif
  }

  /*!
   * Copy constructor
   */
  BOUT_HOST_DEVICE Array(const Array& other) noexcept : ptr(other.ptr) {}

  /*!
   * Assignment operator
   * After this both Arrays share the same dataBlock
   *
   * Uses copy-and-swap idiom
   */
  BOUT_HOST_DEVICE inline Array& operator=(Array other) noexcept {
    swap(*this, other);
    return *this;
  }

  /*!
   * Move constructor
   */
  BOUT_HOST_DEVICE inline Array(Array&& other) noexcept { swap(*this, other); }

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
  static bool useStore( bool keep_using = true ) noexcept {
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
  BOUT_HOST_DEVICE inline bool empty() const noexcept { return ptr == nullptr; }

  /*!
   * Return size of the array. Zero if the array is empty.
   */
  BOUT_HOST_DEVICE inline size_type size() const noexcept {
    if (!ptr)
      return 0;

    // Note: std::valarray::size is technically not noexcept, so
    // Array::size shouldn't be either if we're using valarrays -- in
    // practice, it is so this shouldn't matter
    return ptr->size();
  }

  /*!
   * Returns true if the data is unique to this Array.
   *
   */
  bool unique() const noexcept {
    return ptr.use_count() == 1;
  }

  /*!
   * Ensures that this Array does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() {
    if (!ptr || unique())
      return;

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

  BOUT_HOST_DEVICE inline iterator<T> begin() noexcept { return (ptr) ? std::begin(*ptr) : nullptr; }

  BOUT_HOST_DEVICE inline iterator<T> end() noexcept { return (ptr) ? std::end(*ptr) : nullptr; }

  // Const iterators
  BOUT_HOST_DEVICE inline const_iterator<T> begin() const noexcept { return (ptr) ? std::begin(*ptr) : nullptr; }

  BOUT_HOST_DEVICE inline const_iterator<T> end() const noexcept { return (ptr) ? std::end(*ptr) : nullptr; }

  //////////////////////////////////////////////////////////
  // Element access

  /*!
   * Access a data element. This will fail if the Array is empty (so ptr is null),
   * or if ind is out of bounds. For efficiency no checking is performed,
   * so the user should perform checks.
   */
  BOUT_HOST_DEVICE inline T& operator[](size_type ind) {
#ifndef __CUDA_ARCH__
    ASSERT3(0 <= ind && ind < size());
#endif
    return ptr->operator[](ind);
  }

  BOUT_HOST_DEVICE inline const T& operator[](size_type ind) const {
#ifndef __CUDA_ARCH__
    ASSERT3(0 <= ind && ind < size());
#endif
    return ptr->operator[](ind);
  }

  /*!
   * Exchange contents with another Array of the same type.
   * Sizes of the arrays may differ.
   */
  friend void swap(Array<T, Backing, CreationPolicy>& first, Array<T, Backing, CreationPolicy>& second) noexcept {
    using std::swap;
    swap(first.ptr, second.ptr);
  }

private:

  //Type defs to help keep things brief -- which backing do we use
  using dataBlock = Backing;
  using dataPtrType = typename CreationPolicy::template dataPtrType<dataBlock>;
  using CreationPolicy::create; /// Function to create instances of dataPtrType

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
  static storeType& store(bool cleanup=false) {
#ifdef _OPENMP
    static arenaType arena(omp_get_max_threads());
#else
    static arenaType arena(1);
#endif
    if (!cleanup) {
#ifdef _OPENMP
      return arena[omp_get_thread_num()];
#else
      return arena[0];
#endif
    }

    // Clean by deleting all data -- possible that just stores.clear() is
    // sufficient rather than looping over each entry.
    BOUT_OMP(single)
    {
      for (auto &stores : arena) {
        for (auto &p : stores) {
          auto &v = p.second;
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
      // create function inherited from CreationPolicy
      p = this->template create<dataBlock>(len);
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
template <typename T, typename Backing, typename CreationPolicy>
Array<T, Backing> copy(const Array<T, Backing, CreationPolicy>& other) {
  Array<T, Backing, CreationPolicy> a(other);
  a.ensureUnique();
  return a;
}

/// Compare arrays, which may have different backing and creation policies
template <typename T, typename B1, typename B2, typename C1, typename C2>
bool operator==(const Array<T, B1, C1>& lhs, const Array<T, B2, C2>& rhs) {
  return std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

#endif // __ARRAY_H__

