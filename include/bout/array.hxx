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
 * This behaviour can be disabled by calling the static function useStore:
 *
 * Array<dcomplex>::useStore(false); // Disables memory store
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
#pragma omp atomic
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
#pragma omp atomic
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
#pragma omp atomic
    ptr->refs++;

    // Release the old data
    release(old);
    
    return *this;
  }

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

  /*!
   * Holds a static variable which controls whether
   * memory blocks (ArrayData) are put into a store
   * or new/deleted each time. 
   *
   * The variable is initialised to true on first use,
   * but can be set to false by passing "false" as input.
   * Once set to false it can't be changed back to true.
   */
  static bool useStore( bool keep_using = true ) {
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
  void clear() {
    release(ptr);
    ptr = nullptr;
  }

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

  /*!
   * Exchange contents with another Array of the same type.
   * Sizes of the arrays may differ.
   *
   * This is called by the template function swap(Array&, Array&)
   */
  void swap(Array<T> &other) {
    ArrayData* tmp_ptr = ptr;
    ptr = other.ptr;
    other.ptr = tmp_ptr;
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
   *
   * By putting the static store inside a function it is initialised on first use,
   * and doesn't need to be separately declared for each type T
   *
   * Inputs
   * ------
   *
   * @param[in] cleanup   If set to true, deletes all ArrayData and clears the store
   */
  static std::map< int, std::vector<ArrayData* > > & store(bool cleanup=false) {
    static std::map< int, std::vector<ArrayData* > > store = {};
    
    if (!cleanup) {
      return store;
    }
    
    // Clean by deleting all data
    for (auto &p : store) {
      auto &v = p.second;
      for (ArrayData* a : v) {
        delete a;
      }
      v.clear();
    }
    store.clear();
    
    return store;
  }
  
  /*!
   * Returns a pointer to an ArrayData object with no
   * references. This is either from the store, or newly allocated
   */
  ArrayData* get(int len) {
    ArrayData *p;
#pragma omp critical (store)
    {
      std::vector<ArrayData* >& st = store()[len];
      if (!st.empty()) {
        p = st.back();
        st.pop_back();
      } else {
        p = new ArrayData(len);
      }
    }
    return p;
  }
  
  /*!
   * Release an ArrayData object, reducing its reference count by one. 
   * If no more references, then put back into the store.
   */
  void release(ArrayData *d) {
    if (!d)
      return;
    
    // Reduce reference count, and if zero return to store
#pragma omp critical (store)
    {
      if (!--d->refs) {
        if (useStore()) {
          // Put back into store
          store()[d->len].push_back(d);
        } else {
          delete d;
        }
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

/*!
 * Exchange contents of two Arrays of the same type.
 * Sizes of the arrays may differ.
 */
template<typename T>
void swap(Array<T> &a, Array<T> &b) {
  a.swap(b);
}

#endif // __ARRAY_H__

