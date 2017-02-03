/*
 * Handle arrays of data
 * 
 * Provides an interface to create, iterate over and release
 * arrays of templated tyoes.
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

#include <map>
#include <vector>

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
 * the underlying memory (ArrayData) is put into
 * a map, rather than being freed. 
 * If the same size arrays are used repeatedly then this
 * avoids the need to use new and delete.
 * 
 */
template<typename T>
class Array {
public:
  typedef T data_type;
  
  /*!
   * Create an empty array
   * 
   * Array a();
   * a.empty(); // True
   * 
   */
  Array() : ptr(nullptr) {}
  
  /*!
   * Create an array of given length
   */
  Array(int len) {
    ptr = get(len);
    ptr->refs++;
  }
  
  /*!
   * Destructor. Releases the underlying ArrayData
   */
  ~Array() {
    release(ptr);
  }
  
  /*!
   * Copy constructor
   */
  Array(const Array &other) {
    ptr = other.ptr;
    ptr->refs++;
  }

  /*!
   * Assignment operator
   * After this both Arrays share the same ArrayData
   */
  Array& operator=(const Array &other) {
    ArrayData* const old = ptr;

    // Add reference
    ptr = other.ptr;
    ptr->refs++;

    // Release the old data
    release(old);
    
    return *this;
  }
  
#if __cplusplus >= 201103L
  /*!
   * Move constructor
   */
  Array(Array&& other) {
    ptr = other.ptr;
    other.ptr = nullptr;
  }

  /*! 
   * Move assignment
   */
  Array& operator=(Array &&other) {
    ArrayData* const old = ptr;

    ptr = other.ptr;
    other.ptr = nullptr;

    release(old);
    
    return *this;
  }
#endif

  /*!
   * Release data. After this the Array is empty and any data access
   * will be invalid
   */
  void clear() {
    release(ptr);
    ptr = nullptr;
  }

  /*!
   * Delete all data from the store
   */
  static void cleanup() {
    for(auto &p : store) {
      auto &v = p.second;
      for(ArrayData* a : v) {
        delete a;
      }
      v.clear();
    }
    // Don't use the store anymore
    use_store = false;
  }
  
  /*!
   * Returns true if the Array is empty
   */
  bool empty() const {
    return !ptr;
  }

  /*!
   * Return size of the array. Zero if the array is empty.
   */
  int size() const {
    if(!ptr)
      return 0;

    return ptr->len;
  }
  
  /*!
   * Returns true if the data is unique to this Array.
   * 
   */
  bool unique() const {
    return ptr->refs == 1;
  }

  /*!
   * Ensures that this Array does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() {
    if(!ptr || unique())
      return;

    // Get a new (unique) block of data
    ArrayData* p = get(size());

    // Copy values
    for(iterator it = begin(), ip = p->begin(); it != end(); ++it, ++ip)
      *ip = *it;

    ArrayData *old = ptr;
    ptr = p;
    ptr->refs++;
    
    release(old);
  }
  
  //////////////////////////////////////////////////////////
  // Iterators

  typedef T* iterator;

  iterator begin() {
    return (ptr) ? ptr->data : nullptr;
  }

  iterator end() {
    return (ptr) ? ptr->data + ptr->len : nullptr;
  }

  // Const iterators
  typedef const T* const_iterator;
  
  const_iterator begin() const {
    return (ptr) ? ptr->data : nullptr;
  }

  const_iterator end() const {
    return (ptr) ? ptr->data + ptr->len : nullptr;
  }

  //////////////////////////////////////////////////////////
  // Element access

  /*!
   * Access a data element. This will fail if the Array is empty (so ptr is null),
   * or if ind is out of bounds. For efficiency no checking is performed,
   * so the user should perform checks.
   */
  T& operator[](int ind) {
    return ptr->data[ind];
  }
  const T& operator[](int ind) const {
    return ptr->data[ind];
  }
private:

  /*!
   * ArrayData holds the actual data, and reference count
   * Handles the allocation and deletion of data
   */
  struct ArrayData {
    int refs;   ///< Number of references to this data
    int len;    ///< Size of the array
    T *data;    ///< Array of data
    
    ArrayData(int size) : refs(0), len(size) {
      data = new T[len];
    }
    ~ArrayData() {
      delete[] data;
    }
    iterator begin() {
      return data;
    }
    iterator end() {
      return data + len;
    }
  };

  /*!
   * Pointer to the ArrayData object owned by this Array. 
   * May be null
   */
  ArrayData* ptr;

  /*!
   * This maps from array size (int) to vectors of pointers to ArrayData objects
   * For each data type T, an instance should be declared once (and once only)
   */
  static std::map< int, std::vector<ArrayData* > > store;
  static bool use_store; ///< Should the store be used?
  
  /*!
   * Returns a pointer to an ArrayData object with no
   * references. This is either from the store, or newly allocated
   */
  ArrayData* get(int len) {
    std::vector<ArrayData* >& st = store[len];
    if(st.empty()) {
      return new ArrayData(len);
    }
    ArrayData *p = st.back();
    st.pop_back();
    return p;
  }

  /*!
   * Release an ArrayData object, reducing its reference count by one. 
   * If no more references, then put back into the store.
   */
  void release(ArrayData *d) {
    if(!d)
      return;
    
    // Reduce reference count, and if zero return to store
    if(!--d->refs) {
      if(use_store) {
        // Put back into store
        store[d->len].push_back(d);
      }else {
        delete d;
      }
    }
  }
  
};

/*!
 * Create a copy of an Array, which does not share data
 */ 
template<typename T>
Array<T>& copy(const Array<T> &other) {
  Array<T> a(other);
  a.ensureUnique();
  return a;
}

#endif // __ARRAY_H__

