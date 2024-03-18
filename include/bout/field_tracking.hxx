#pragma once

template <class F, class D>
class FieldTracking : public F {
public:
  using F::F;
  FieldTracking(const F& f) : F(f){};

  template <class... Args>
  FieldTracking& operator=(const Args&... args) {
    F::operator=(args...);
    static_cast<D*>(this)->onChange();
  }

  template <class... Args>
  FieldTracking& operator*=(const Args&... args) {
    F::operator/=(args...);
    static_cast<D*>(this)->onChange();
    return *this;
  }

  template <class... Args>
  FieldTracking& operator/=(const Args&... args) {
    F::operator/=(args...);
    static_cast<D*>(this)->onChange();
    return *this;
  }

  template <class... Args>
  FieldTracking& operator+=(const Args&... args) {
    F::operator/=(args...);
    static_cast<D*>(this)->onChange();
    return *this;
  }

  template <class... Args>
  FieldTracking& operator-=(const Args&... args) {
    F::operator/=(args...);
    static_cast<D*>(this)->onChange();
    return *this;
  }

  operator const F&() const { return static_cast<F*>(this); }
};

template <class F, class D>
F operator*(FieldTracking<F, D> f1, F f2) {
  return *static_cast<F*>(&f1) * f2;
}

template <class F, class D>
F operator*(FieldTracking<F, D> f1, FieldTracking<F, D> f2) {
  return *static_cast<F*>(&f1) * *static_cast<F*>(&f2);
}

template <class F, class D>
F operator*(F f1, FieldTracking<F, D> f2) {
  return f1 * *static_cast<F*>(&f2);
}

template <class F, class D>
F operator/(FieldTracking<F, D> f1, F f2) {
  return *static_cast<F*>(&f1) / f2;
}

template <class F, class D>
F operator/(FieldTracking<F, D> f1, FieldTracking<F, D> f2) {
  return *static_cast<F*>(&f1) / *static_cast<F*>(&f2);
}

template <class F, class D>
F operator/(F f1, FieldTracking<F, D> f2) {
  return f1 / *static_cast<F*>(&f2);
}

template <class F, class D>
F operator+(FieldTracking<F, D> f1, F f2) {
  return *static_cast<F*>(&f1) + f2;
}

template <class F, class D>
F operator+(FieldTracking<F, D> f1, FieldTracking<F, D> f2) {
  return *static_cast<F*>(&f1) + *static_cast<F*>(&f2);
}

template <class F, class D>
F operator+(F f1, FieldTracking<F, D> f2) {
  return f1 + *static_cast<F*>(&f2);
}

template <class F, class D>
F operator-(FieldTracking<F, D> f1, F f2) {
  return *static_cast<F*>(&f1) - f2;
}

template <class F, class D>
F operator-(FieldTracking<F, D> f1, FieldTracking<F, D> f2) {
  return *static_cast<F*>(&f1) - *static_cast<F*>(&f2);
}

template <class F, class D>
F operator-(F f1, FieldTracking<F, D> f2) {
  return f1 - *static_cast<F*>(&f2);
}
