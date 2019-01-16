#ifndef __ARRAYND_HXX__
#define __ARRAYND_HXX__

#include "output.hxx"
#include "bout/array.hxx"

template <typename T, int ndim>
class ArrayND {
public:
  using data_type = T;
  using slice_type = ArrayND<T, ndim - 1>;
  using size_type = int;
  using shape_type = decltype(
      std::tuple_cat(std::tuple<size_type>{}, typename slice_type::shape_type{}));

  constexpr static int ndims = ndim;

  ArrayND() : n1(0){};

  template <typename... allSizes>
  ArrayND(size_type n1, allSizes... sizes) : n1(n1) {
    static_assert(sizeof...(sizes) == ndim - 1,
                  "Incorrect number of dimension sizes passed as arguments to ArrayND.");

    data = Array<slice_type>(n1);
    for (auto& i : data) {
      i = slice_type(sizes...);
    }
  }
  ArrayND(const ArrayND& other) : n1(other.n1), data(other.data) {
    // Prevent copy on write for ArrayND
    data.ensureUnique();
  }

  ArrayND& operator=(const ArrayND& other) {
    n1 = other.n1;
    data = other.data;
    // Prevent copy on write for ArrayND
    data.ensureUnique();
    return *this;
  }
  template <typename... allSizes>
  inline T& operator()(size_type i1, allSizes... sizes) {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1](sizes...);
  }
  inline slice_type& operator[](size_type i1) {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1];
  }
  inline slice_type& operator()(size_type i1) {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1];
  }
  template <typename... allSizes>
  inline const T& operator()(size_type i1, allSizes... sizes) const {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1](sizes...);
  }
  inline const slice_type& operator[](size_type i1) const {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1];
  }
  inline const slice_type& operator()(size_type i1) const {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1];
  }

  ArrayND& operator=(const T& val) {
    for (auto& i : data) {
      i = val;
    };
    return *this;
  };

  slice_type* begin() { return std::begin(data); };
  const slice_type* begin() const { return std::begin(data); };
  slice_type* end() { return std::end(data); };
  const slice_type* end() const { return std::end(data); };

  auto shape() const -> shape_type {
    return std::tuple_cat(std::make_tuple(n1), data[0].shape());
  }

  size_type size() const { return n1 * data[0].size(); }

  bool empty() const { return size() == 0; }

  /*!
   * Ensures that this ArrayND does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() { data.ensureUnique(); }

private:
  size_type n1;
  Array<slice_type> data;
};

template <typename T>
class ArrayND<T, 1> {
public:
  using data_type = T;
  using size_type = int;
  using shape_type = std::tuple<size_type>;
  constexpr static int ndim = 1;

  constexpr static int ndims = ndim;

  ArrayND() : n1(0){};
  ArrayND(size_type n1) : n1(n1) { data = Array<T>(n1); }
  // Should only end up calling this if we pass too many dimension sizes
  template <typename... allSizes>
  ArrayND(size_type n1, allSizes... sizes) : n1(n1) {
    static_assert(
        sizeof...(sizes) == ndim - 1,
        "Incorrect number of dimension sizes passed as arguments to ArrayND<T, 1>.");
  }

  ArrayND(const ArrayND& other) : n1(other.n1), data(other.data) {
    // Prevent copy on write for ArrayND
    data.ensureUnique();
  }

  ArrayND& operator=(const ArrayND& other) {
    n1 = other.n1;
    data = other.data;
    // Prevent copy on write for ArrayND
    data.ensureUnique();
    return *this;
  }

  inline T& operator()(size_type i1) {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1];
  }
  inline T& operator[](size_type i1) {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1];
  }
  inline const T& operator()(size_type i1) const {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1];
  }
  inline const T& operator[](size_type i1) const {
    ASSERT2(0 <= i1 && i1 < n1);
    return data[i1];
  }

  ArrayND& operator=(const T& val) {
    for (auto& i : data) {
      i = val;
    };
    return *this;
  };

  T* begin() { return std::begin(data); };
  const T* begin() const { return std::begin(data); };
  T* end() { return std::end(data); };
  const T* end() const { return std::end(data); };

  std::tuple<size_type> shape() const { return std::make_tuple(n1); }

  size_type size() const { return n1; }

  bool empty() const { return size() == 0; }

  /*!
   * Ensures that this ArrayND does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() { data.ensureUnique(); }

private:
  size_type n1;
  Array<T> data;
};

template <typename T>
using Array1D = ArrayND<T, 1>;

template <typename T>
using Array2D = ArrayND<T, 2>;

template <typename T>
using Array3D = ArrayND<T, 3>;

#endif
