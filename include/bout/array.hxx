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
#include <vector>
#include <memory>
#include <numeric>
#include <mutex>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(BOUT_USE_CUDA) && defined(__CUDACC__)
#define BOUT_HOST_DEVICE __host__ __device__
#define BOUT_HOST __host__
#define BOUT_DEVICE __device__
#else
#define BOUT_HOST_DEVICE
#define BOUT_HOST
#define BOUT_DEVICE
#endif


#ifdef BOUT_HAS_UMPIRE
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

template <typename T>
struct ArrayData {
   ArrayData(int size) : len(size), clientUseCount(1), owner(true) { 
#ifdef BOUT_HAS_UMPIRE 
     auto& rm = umpire::ResourceManager::getInstance();
#ifdef BOUT_USE_CUDA
     //auto allocator = rm.getAllocator("UM");
     auto allocator = rm.getAllocator(umpire::resource::Pinned);
#else
     auto allocator = rm.getAllocator("HOST");
#endif
     data = static_cast<T*>(allocator.allocate(size * sizeof(T)));
#else
     data = new T[len]; 
#endif
  }

  
BOUT_HOST_DEVICE   ~ArrayData() { 
// __CUDA_ARCH__ is only defined device side
#ifndef __CUDA_ARCH__
   //  printf("UM dealloc len %d owner %d count %d\n",len,owner,clientUseCount);
#ifdef BOUT_HAS_UMPIRE
     if(data != nullptr && owner && (clientUseCount == 1)) { 
    //  printf("UM dealloc %d %p\n",len,data);
      auto& rm = umpire::ResourceManager::getInstance();
      rm.deallocate(data);
      data = nullptr;
      len = 0;
      clientUseCount = 0;
      owner = false;
     }
#else
     if(data != nullptr && owner && (clientUseCount == 1)) {
      delete[] data;
      data = nullptr;
      clientUseCount = 0;
      owner = false;
     }
#endif
#else
   //printf("device ~ArrayData()\n");     
#endif 
  }

BOUT_HOST_DEVICE ArrayData(const ArrayData& other) {
   len = other.len;
   data = other.data;
   owner = false;
}

BOUT_HOST_DEVICE inline iterator<T> begin() const { return data; }
BOUT_HOST_DEVICE inline iterator<T> end() const { return data + len; }
BOUT_HOST_DEVICE inline int size() const { return len; }
BOUT_HOST_DEVICE inline void dec_count() {clientUseCount--;}
BOUT_HOST_DEVICE inline void inc_count() {clientUseCount++;}
BOUT_HOST_DEVICE inline int use_count() { return clientUseCount;}

BOUT_HOST_DEVICE inline void operator=(ArrayData<T>& in) { 
#ifndef __CUDA_ARCH__
   std::copy(std::begin(in), std::end(in), begin()); 
#else
   int threadId = threadIdx.x + blockDim.x * threadIdx.y +
                 (blockDim.x * blockDim.y) * threadIdx.z;

   if(threadId < len) {
      data[threadId] = in[threadId];
   }
#endif
}

// this is the operator likely to be called on the device
BOUT_HOST_DEVICE inline void operator=(const ArrayData<T>& in) const { 
#ifndef __CUDA_ARCH__
   std::copy(std::begin(in), std::end(in), begin()); 
#else
   int threadId = threadIdx.x + blockDim.x * threadIdx.y +
                 (blockDim.x * blockDim.y) * threadIdx.z;

   if(threadId < len) {
      data[threadId] = in[threadId];
   }
#endif
}

BOUT_HOST_DEVICE inline T& operator[](int ind) { return data[ind]; };
BOUT_HOST_DEVICE inline  const T operator[](int ind) const { return data[ind]; };
private:
   T* data=nullptr;
   int len=0;
   int clientUseCount=0;
   bool owner = false;
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

template<typename T, typename Backing = ArrayData<T>>
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
  Array(size_type len) {
    ptr = get(len);

 //   printf("Array create %d @ %p\n",len,ptr->begin());  
  }
  
  /*!
   * Destructor. Releases the underlying dataBlock
   */
BOUT_HOST_DEVICE inline   ~Array() {
#ifndef __CUDA_ARCH__
    if(ptr) {
   //   printf("~Array count %d data %p\n",ptr->use_count(),ptr->begin());  
      release(ptr);
    }
#else
    //printf("~Array Device\n");
#endif
  }
  
  /*!
   * Copy constructor
   */
BOUT_HOST_DEVICE Array(const Array &other) {
    ptr = other.ptr;
#ifndef __CUDA_ARCH__
    ptr->inc_count();
#endif
  }

  /*!
   * Assignment operator
   * After this both Arrays share the same dataBlock
   *
   * Uses copy-and-swap idiom
   */

BOUT_HOST  inline Array& operator=(Array other) noexcept {
    swap(*this, other);
    return *this;
  }

BOUT_DEVICE  inline void operator=(const Array other) const {
    ptr->operator=(*other.ptr);
  }

  /*!
   * Move constructor
   */
  Array(Array&& other) noexcept {
    swap(*this, other);
  }

  /*!
   * Reallocate the array with size = \p new_size
   *
   * Note that this invalidates the existing data!
   */
  void reallocate(size_type new_size) {
   // printf("Array realloc attempt %d\n",new_size);
    int old_size = 0;
    if(ptr) {
      old_size = ptr->size();
      release(ptr);
    }
    ptr = get(new_size);
   // printf("Array realloc from %d to %d @ %p\n",old_size,new_size,ptr->begin());  
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
    static bool value = false;
    return value;
  }
  /*!
   * Release data. After this the Array is empty and any data access
   * will be invalid
   */
  void clear() noexcept {
    release(ptr);
  }


  static void cleanup() {
     // maybe do some umpire pool cleanup
  } 

  /*!
   * Returns true if the Array is empty
   */
BOUT_HOST_DEVICE inline   bool empty() const noexcept {
    return ptr == nullptr;
  }

  /*!
   * Return size of the array. Zero if the array is empty.
   */
BOUT_HOST_DEVICE inline   size_type size() const noexcept {
    if(!ptr)
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
   // printf("Array unique()\n");
    if(ptr) {
      return ptr->use_count() == 1;
    }
    return true;
  }

  /*!
   * Ensures that this Array does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() {
   // printf("Array ensureUnique\n");
    if(!ptr || unique())
      return;

    // Get a new (unique) block of data
    dataPtrType p = get(size());

    //Make copy of the underlying data
    p->operator=((*ptr));

    //Update the local pointer and release old
    //Can't just do ptr=p as need to try to add to store.
    release(ptr);
    ptr = std::move(p);
  }

  //////////////////////////////////////////////////////////
  // Iterators

BOUT_HOST_DEVICE inline iterator<T> begin() noexcept { 
   if(ptr) {
      return (ptr->begin());
   }
   else {
      return(nullptr);
   }
}

BOUT_HOST_DEVICE inline iterator<T> end() noexcept { 
   if(ptr) {
      return (ptr->end());
   }
   else {
      return(nullptr);
   }
}

BOUT_HOST_DEVICE inline const_iterator<T> begin() const noexcept { 
   if(ptr) {
      return (ptr->begin());
   }
   else {
      return(nullptr);
   }
}

BOUT_HOST_DEVICE inline const_iterator<T> end() const noexcept { 
   if(ptr) {
      return (ptr->end());
   }
   else {
      return(nullptr);
   }
}

  //////////////////////////////////////////////////////////
  // Element access

  /*!
   * Access a data element. This will fail if the Array is empty (so ptr is null),
   * or if ind is out of bounds. For efficiency no checking is performed,
   * so the user should perform checks.
   */
BOUT_HOST_DEVICE inline   T& operator[](size_type ind) {
    return ptr->operator[](ind);
  }

BOUT_HOST_DEVICE inline   const T& operator[](size_type ind) const {
    return ptr->operator[](ind);
  }

  /*!
   * Exchange contents with another Array of the same type.
   * Sizes of the arrays may differ.
   */
  friend void swap(Array<T> &first, Array<T> &second) noexcept {
    using std::swap;
    swap(first.ptr, second.ptr);
  }

private:

  //Type defs to help keep things brief -- which backing do we use
  using dataBlock = Backing;
  using dataPtrType = dataBlock*;
  /*!
   * Pointer to the data container object owned by this Array. 
   * May be null
   */
  dataPtrType ptr = nullptr;


  /*!
   * Returns a pointer to a dataBlock object of size \p len with no
   * references. This is either from the store, or newly allocated
   *
   * Expects \p len >= 0
   */
  dataPtrType get(size_type len) {
#ifdef BOUT_HAS_UMPIRE 
     auto& rm = umpire::ResourceManager::getInstance();
#ifdef BOUT_USE_CUDA
     auto allocator = rm.getAllocator("UM");
#else
     auto allocator = rm.getAllocator("HOST");
#endif
     dataPtrType pp = static_cast<dataPtrType>(allocator.allocate(sizeof(dataBlock)));
     dataPtrType p = new (pp) dataBlock(len);
 //    printf("UM  Array %d %p %p\n",len,p,pp);
#else
     dataPtrType p = new dataBlock(len);
#endif
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

  void release(dataPtrType& d) {
    if (!d) {
      return;
    }

   // printf("release size %d count %d data %p dataBlock %p\n",d->size(),d->use_count(),d->begin(),d);
    // Reduce reference count, and if zero return to store
    if (d->use_count() == 1) {
#ifdef BOUT_HAS_UMPIRE
      auto& rm = umpire::ResourceManager::getInstance();
      rm.deallocate(d->begin());
      rm.deallocate(d);
      umpire::Allocator alloc = rm.getAllocator("UM");
    //  printf("Umpire UM Allocations %d\n",alloc.getAllocationCount());
#else
      delete d;
#endif
      d = nullptr;
    } else {
      d->dec_count();
    }
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
#endif // __ARRAY_H__

