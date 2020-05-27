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

  ArrayND() : len(0){};

  template <typename... allSizes>
  ArrayND(size_type len, allSizes... sizes) : len(len) {
    static_assert(sizeof...(sizes) == ndim - 1,
                  "Incorrect number of dimension sizes passed as arguments to ArrayND.");

    data = Array<slice_type>(len);
    for (auto& i : data) {
      i = slice_type(sizes...);
    }
  }
  ArrayND(const ArrayND& other) : len(other.len), data(other.data) {
  }

  ArrayND& operator=(const ArrayND& other) {
    len = other.len;
    data = other.data;
    return *this;
  }

  template <typename... allSizes>
  inline T& operator()(size_type i1, allSizes... sizes) {
    ASSERT2(0 <= i1 && i1 < len);
    return data[i1](sizes...);
  }
  inline slice_type& operator[](size_type i1) {
    ASSERT2(0 <= i1 && i1 < len);
    return data[i1];
  }
  inline slice_type& operator()(size_type i1) {
    ASSERT2(0 <= i1 && i1 < len);
    return data[i1];
  }

  template <typename... allSizes>
  inline const T& operator()(size_type i1, allSizes... sizes) const {
    ASSERT2(0 <= i1 && i1 < len);
    return data[i1](sizes...);
  }
  inline const slice_type& operator[](size_type i1) const {
    ASSERT2(0 <= i1 && i1 < len);
    return data[i1];
  }
  inline const slice_type& operator()(size_type i1) const {
    ASSERT2(0 <= i1 && i1 < len);
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

  // Note we assume a non-jagged array, i.e. we assume
  // data[0].shape() == data[j].shape() for 0 <= j < len
  constexpr auto shape() const -> shape_type {
    return std::tuple_cat(std::make_tuple(len), data[0].shape());
  }

  // Note we assume a non-jagged array, i.e. we assume
  // data[0].size() == data[j].size() for 0 <= j < len
  size_type size() const { return len * data[0].size(); }

  bool empty() const { return size() == 0; }

  bool unique() const {
    bool res = data.unique();
    for (const auto& i : data) {
      res = res && i.unique();
      // Try to return early if we can.
      if (!res)
        return res;
    };
    return res;
  }

  void ensureUnique() {
    data.ensureUnique();
    for (auto& i : data) {
      i.ensureUnique();
    };
  }

  void pack(ArrayND<T, 1>& flat) {
    ASSERT1(size() == flat.size());
    size_type currentPosition = 0;
    currentPosition = setFromFlatArray(flat, currentPosition);
    ASSERT1(currentPosition == size());
  };

  void pack(Array<T>& flat) {
    ASSERT1(size() == flat.size());
    size_type currentPosition = 0;
    currentPosition = setFromFlatArray(flat, currentPosition);
    ASSERT1(currentPosition == size());
  };

  size_type setFromFlatArray(Array<T>& flat, size_type currentPos = 0) {
    ArrayND<T, 1> tmp(flat.size());
    tmp.data = flat;
    return this->setFromFlatArray(tmp, currentPos);
  }

  size_type setFromFlatArray(ArrayND<T, 1>& flat, size_type currentPos = 0) {
    ASSERT1(flat.size() >= size());

    for (auto& i : data) {
      currentPos = i.setFromFlatArray(flat, currentPos);
    }

    return currentPos;
  };

  ArrayND<T, 1> flatten() {
    ArrayND<T, 1> result(size());
    size_type start = 0;
    for (auto& i : data) {
      result.setSpan(start, i.size(), i.flatten());
      start += i.size();
    }

    return result;
  };

  Array<slice_type> data;

private:
  size_type len;
};

template <typename T>
class ArrayND<T, 1> {
public:
  using data_type = T;
  using size_type = int;
  using shape_type = std::tuple<size_type>;
  constexpr static int ndim = 1;

  constexpr static int ndims = ndim;

  ArrayND() : len(0){};
  ArrayND(size_type len) : len(len) { data = Array<T>(len); }
  // Should only end up calling this if we pass too many dimension sizes
  template <typename... allSizes>
  ArrayND(size_type len, allSizes... sizes) : len(len) {
    static_assert(
        sizeof...(sizes) == ndim - 1,
        "Incorrect number of dimension sizes passed as arguments to ArrayND<T, 1>.");
  }

  ArrayND(const ArrayND& other) : len(other.len), data(other.data) {
  }

  ArrayND& operator=(const ArrayND& other) {
    len = other.len;
    data = other.data;
    return *this;
  }

  inline T& operator()(size_type i1) {
    ASSERT2(0 <= i1 && i1 < len);
    return data[i1];
  }
  inline T& operator[](size_type i1) {
    ASSERT2(0 <= i1 && i1 < len);
    return data[i1];
  }
  inline const T& operator()(size_type i1) const {
    ASSERT2(0 <= i1 && i1 < len);
    return data[i1];
  }
  inline const T& operator[](size_type i1) const {
    ASSERT2(0 <= i1 && i1 < len);
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

  constexpr std::tuple<size_type> shape() const { return std::make_tuple(len); }

  size_type size() const { return len; }

  bool empty() const { return size() == 0; }

  bool unique() const { return data.unique(); }

  /*!
   * Ensures that this ArrayND does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() { data.ensureUnique(); }

  size_type setFromFlatArray(Array<T>& flat, size_type currentPos = 0) {
    ArrayND<T, 1> tmp(flat.size());
    tmp.data = flat;
    return this->setFromFlatArray(tmp, currentPos);
  }

  size_type setFromFlatArray(ArrayND<T, 1>& flat, size_type currentPos = 0) {
    ASSERT1(flat.size() >= size());
    setSpan(0, size(), flat.getSpan(currentPos, size()));
    return currentPos + size();
  };

  ArrayND<T, 1> flatten() { return *this; }

  void setSpan(size_type start, size_type count, ArrayND<T, 1> source) {
    this->setSpan(start, count, source.data);
  }

  void setSpan(size_type start, size_type count, Array<T> source) {
    ASSERT1(start >= 0);
    ASSERT1(start < len);
    ASSERT1(start + count - 1 <= len);
    ASSERT1(source.size() >= count);
    for (size_type i = 0; i < count; i++) {
      data[start + i] = source[i];
    }
  };

  Array<T> getSpan(size_type start, size_type count) {
    ASSERT1(start >= 0);
    ASSERT1(start < len);
    ASSERT1(start + count - 1 <= len);
    Array<T> result{count};
    for (size_type i = 0; i < count; i++) {
      result[i] = data[start + i];
    }
    return result;
  };

  Array<T> data;

private:
  size_type len;
};

template <typename T>
using Array1D = ArrayND<T, 1>;

template <typename T>
using Array2D = ArrayND<T, 2>;

template <typename T>
using Array3D = ArrayND<T, 3>;

#endif
