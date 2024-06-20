#pragma once

#include "bout/bout_types.hxx"
#include <memory>
#include <string>

namespace bout {

class Scalar {
  BoutReal data{};
  std::unique_ptr<Scalar> deriv{nullptr};

public:
  Scalar(BoutReal value) : data(value) {}

  static int size() { return 1; }

  Scalar& timeDeriv() {
    if (deriv == nullptr) {
      deriv = std::make_unique<Scalar>(0.0);
    }
    return *deriv;
  }

  Scalar& operator=(BoutReal value) {
    data = value;
    return *this;
  }
  operator BoutReal() const { return data; }

  std::string name;
};

inline Scalar& ddt(Scalar& scalar) { return scalar.timeDeriv(); }

} // namespace bout
