#include <complex>
#include "bout_types.hxx"
#include "bout/array.hxx"


template<typename T, typename Backing>
bool Array<T, Backing>::useStore( bool keep_using) noexcept{
  static bool value = true;
  if (keep_using) {
    return value; 
  }
  // Change to false
  value = false;
  return value;
}

template<typename T, typename Backing>
void Array<T, Backing>::cleanup() {
  // Clean the store, deleting data
  store(true);
  // Don't use the store anymore. This is so that array releases
  // after cleanup() get deleted rather than put into the store
  useStore(false);
}

template<typename T, typename Backing>
void Array<T, Backing>::ensureUnique() {
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

template<typename T, typename Backing>
void Array<T, Backing>::release(dataPtrType& d) noexcept {
  if (!d)
    return;

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

template<typename T, typename Backing>
typename Array<T, Backing>::storeType& Array<T, Backing>::store(bool cleanup) {
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

  //Store should now be empty but we need to return something,
  //so return an empty storeType from the arena.
  return arena[0];
}

template<typename T, typename Backing>
typename Array<T, Backing>::dataPtrType Array<T, Backing>::get(size_type len) {
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


template<typename T, typename Backing>
Array<T, Backing>::~Array<T, Backing>() noexcept{
  release(ptr);
}

template class Array<int, ArrayData<int>>;
template class Array<std::complex<BoutReal>, ArrayData<std::complex<BoutReal>>>;
template class Array<BoutReal, ArrayData<BoutReal>>;


