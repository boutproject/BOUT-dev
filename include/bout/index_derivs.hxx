/*!************************************************************************
 * \file index_derivs.hxx
 *
 * Definition of available derivative methods and registration within store
 *
 **************************************************************************
 * Copyright 2018
 *    D.Dickinson, P.Hill, B.Dudson
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef __INDEX_DERIVS_HXX__
#define __INDEX_DERIVS_HXX__

#include <functional>
#include <iostream>

#include <bout/assert.hxx>
#include <bout/constants.hxx>
#include <bout/deriv_store.hxx>
#include <bout/index_derivs_interface.hxx>
#include <bout/region.hxx>
#include <bout/scorepwrapper.hxx>
#include <bout/template_combinations.hxx>

#include <bout_types.hxx>
#include <fft.hxx>
#include <interpolation.hxx>
#include <msg_stack.hxx>
#include <stencils.hxx>
#include <unused.hxx>

class Field3D;
class Field2D;

const BoutReal WENO_SMALL = 1.0e-8; // Small number for WENO schemes

struct metaData {
  // Would rather use a std::string here but this causes the
  // metaData struct to be non-trivially destrucible which
  // can prevent using temporary instances of this. Instead
  // we'll use char* for now.
  // const std::string key;
  const char* key;
  const int nGuards;
  const DERIV derivType; // Can be used to identify the type of the derivative
};

/// Provide an easy way to report a Region's statistics
inline std::ostream& operator<<(std::ostream& out, const metaData& meta) {
  out << "key : " << meta.key;
  out << ", ";
  out << "nGuards : " << meta.nGuards;
  out << ", ";
  out << "type : " << toString(meta.derivType);
  return out;
}

/// Here we define a helper class that provides a means to use a supplied
/// stencil using functor to calculate a derivative over the entire field.
/// Note we currently have a different interface for some of the derivative types
/// to avoid needing different classes to represent the different operations
/// The use of a functor here makes it possible to wrap up metaData into the
/// type as well.
template <typename FF>
class DerivativeType {
public:
  template <DIRECTION direction, STAGGER stagger, int nGuards, typename T>
  void standard(const T& var, T& result, const std::string& region) const {
    AUTO_TRACE();
    ASSERT2(meta.derivType == DERIV::Standard || meta.derivType == DERIV::StandardSecond
            || meta.derivType == DERIV::StandardFourth)
    ASSERT2(var.getMesh()->getNguard(direction) >= nGuards);

    BOUT_FOR(i, var.getRegion(region)) {
      result[i] = apply(populateStencil<direction, stagger, nGuards>(var, i));
    }
    return;
  }

  template <DIRECTION direction, STAGGER stagger, int nGuards, typename T>
  void upwindOrFlux(const T& vel, const T& var, T& result, const std::string& region) const {
    AUTO_TRACE();
    ASSERT2(meta.derivType == DERIV::Upwind || meta.derivType == DERIV::Flux)
    ASSERT2(var.getMesh()->getNguard(direction) >= nGuards);

    if (meta.derivType == DERIV::Flux || stagger != STAGGER::None) {
      BOUT_FOR(i, var.getRegion(region)) {
        result[i] = apply(populateStencil<direction, stagger, nGuards>(vel, i),
                          populateStencil<direction, STAGGER::None, nGuards>(var, i));
      }
    } else {
      BOUT_FOR(i, var.getRegion(region)) {
        result[i] =
            apply(vel[i], populateStencil<direction, STAGGER::None, nGuards>(var, i));
      }
    }
    return;
  }

  static constexpr FF func{};
  static constexpr metaData meta{FF::meta};

  BoutReal apply(const stencil& f) const { return func(f); }
  BoutReal apply(BoutReal v, const stencil& f) const { return func(v, f); }
  BoutReal apply(const stencil& v, const stencil& f) const { return func(v, f); }
};

// Redundant definitions because C++
// Not necessary in C++17
template <class FF>
constexpr FF DerivativeType<FF>::func;
template <class FF>
constexpr metaData DerivativeType<FF>::meta;

/////////////////////////////////////////////////////////////////////////////////
/// Following code is for dealing with registering a method/methods for all
/// template combinations, in conjunction with the template_combinations code.
/////////////////////////////////////////////////////////////////////////////////

struct registerMethod {
  template <typename Direction, typename Stagger, typename FieldTypeContainer,
            typename Method>
  void operator()(Direction, Stagger, FieldTypeContainer, Method) {
    AUTO_TRACE();
    using namespace std::placeholders;

    // Now we want to get the actual field type out of the TypeContainer
    // used to pass this around
    using FieldType = typename FieldTypeContainer::type;

    Method method{};

    // Note whilst this should be known at compile time using this directly in the
    // template parameters below causes problems for old versions of gcc/libstdc++
    // (tested with 4.8.3) so we currently use a hacky workaround. Once we drop
    // support for these versions the branching in the case statement below can be
    // removed and we can use nGuard directly in the template statement.
    const int nGuards = method.meta.nGuards;

    auto& derivativeRegister = DerivativeStore<FieldType>::getInstance();

    switch (method.meta.derivType) {
    case (DERIV::Standard):
    case (DERIV::StandardSecond):
    case (DERIV::StandardFourth): {
      if (nGuards == 1) {
        const auto theFunc = std::bind(
            // Method to store in function
            &Method::template standard<Direction::value, Stagger::value, 1, FieldType>,
            // Arguments -- first is hidden this of type-bound, others are placeholders
            // for input field, output field, region
            method, _1, _2, _3);
        derivativeRegister.registerDerivative(theFunc, Direction{}, Stagger{}, method);
      } else {
        const auto theFunc = std::bind(
            // Method to store in function
            &Method::template standard<Direction::value, Stagger::value, 2, FieldType>,
            // Arguments -- first is hidden this of type-bound, others are placeholders
            // for input field, output field, region
            method, _1, _2, _3);
        derivativeRegister.registerDerivative(theFunc, Direction{}, Stagger{}, method);
      }
      break;
    }
    case (DERIV::Upwind):
    case (DERIV::Flux): {
      if (nGuards == 1) {
        const auto theFunc = std::bind(
            // Method to store in function
            &Method::template upwindOrFlux<Direction::value, Stagger::value, 1,
                                           FieldType>,
            // Arguments -- first is hidden this of type-bound, others are placeholders
            // for input field, output field, region
            method, _1, _2, _3, _4);
        derivativeRegister.registerDerivative(theFunc, Direction{}, Stagger{}, method);
      } else {
        const auto theFunc = std::bind(
            // Method to store in function
            &Method::template upwindOrFlux<Direction::value, Stagger::value, 2,
                                           FieldType>,
            // Arguments -- first is hidden this of type-bound, others are placeholders
            // for input field, output field, region
            method, _1, _2, _3, _4);
        derivativeRegister.registerDerivative(theFunc, Direction{}, Stagger{}, method);
      }
      break;
    }
    default:
      throw BoutException("Unhandled derivative method in registerMethod.");
    };
  }
};

#define DEFINE_STANDARD_DERIV_CORE(name, key, nGuards, type)                        \
  struct name {                                                                     \
    BoutReal operator()(const stencil& f) const;                                    \
    BoutReal operator()(BoutReal UNUSED(vc), const stencil& UNUSED(f)) const {      \
      return BoutNaN;                                                               \
    };                                                                              \
    BoutReal operator()(const stencil& UNUSED(v), const stencil& UNUSED(f)) const { \
      return BoutNaN;                                                               \
    };                                                                              \
    static constexpr metaData meta = {key, nGuards, type};                          \
  };
#define DEFINE_STANDARD_DERIV(name, key, nGuards, type) \
  DEFINE_STANDARD_DERIV_CORE(name, key, nGuards, type)  \
  BoutReal name::operator()(const stencil& f) const

#define DEFINE_UPWIND_DERIV_CORE(name, key, nGuards, type)                          \
  struct name {                                                                     \
    BoutReal operator()(const stencil& UNUSED(f)) const { return BoutNaN; };        \
    BoutReal operator()(BoutReal vc, const stencil& f) const;                       \
    BoutReal operator()(const stencil& UNUSED(v), const stencil& UNUSED(f)) const { \
      return BoutNaN;                                                               \
    };                                                                              \
    static constexpr metaData meta = {key, nGuards, type};                          \
  };
#define DEFINE_UPWIND_DERIV(name, key, nGuards, type) \
  DEFINE_UPWIND_DERIV_CORE(name, key, nGuards, type)  \
  BoutReal name::operator()(BoutReal vc, const stencil& f) const

#define DEFINE_FLUX_DERIV_CORE(name, key, nGuards, type)                       \
  struct name {                                                                \
    BoutReal operator()(const stencil& UNUSED(f)) const { return BoutNaN; };   \
    BoutReal operator()(BoutReal UNUSED(vc), const stencil& UNUSED(f)) const { \
      return BoutNaN;                                                          \
    };                                                                         \
    BoutReal operator()(const stencil& v, const stencil& f) const;             \
    static constexpr metaData meta = {key, nGuards, type};                     \
  };
#define DEFINE_FLUX_DERIV(name, key, nGuards, type) \
  DEFINE_FLUX_DERIV_CORE(name, key, nGuards, type)  \
  BoutReal name::operator()(const stencil& v, const stencil& f) const

#define DEFINE_STANDARD_DERIV_STAGGERED(name, key, nGuards, type) \
  DEFINE_STANDARD_DERIV(name, key, nGuards, type)
#define DEFINE_UPWIND_DERIV_STAGGERED(name, key, nGuards, type) \
  DEFINE_FLUX_DERIV(name, key, nGuards, type)
#define DEFINE_FLUX_DERIV_STAGGERED(name, key, nGuards, type) \
  DEFINE_FLUX_DERIV(name, key, nGuards, type)

/// Some helper defines for now that allow us to wrap up enums
/// and the specific methods.
#define WRAP_ENUM(family, value) enumWrapper<family, family::value>

#define REGISTER_DERIVATIVE(name)                                                      \
  namespace {                                                                          \
  produceCombinations<Set<WRAP_ENUM(DIRECTION, X), WRAP_ENUM(DIRECTION, Y),            \
                          WRAP_ENUM(DIRECTION, YOrthogonal), WRAP_ENUM(DIRECTION, Z)>, \
                      Set<WRAP_ENUM(STAGGER, None)>,                                   \
                      Set<TypeContainer<Field3D>, TypeContainer<Field2D>>,             \
                      Set<DerivativeType<name>>>                                       \
      reg##name(registerMethod{});                                                     \
  }
#define REGISTER_STAGGERED_DERIVATIVE(name)                                            \
  namespace {                                                                          \
  produceCombinations<Set<WRAP_ENUM(DIRECTION, X), WRAP_ENUM(DIRECTION, Y),            \
                          WRAP_ENUM(DIRECTION, YOrthogonal), WRAP_ENUM(DIRECTION, Z)>, \
                      Set<WRAP_ENUM(STAGGER, C2L), WRAP_ENUM(STAGGER, L2C)>,           \
                      Set<TypeContainer<Field3D>, TypeContainer<Field2D>>,             \
                      Set<DerivativeType<name>>>                                       \
      reg##name(registerMethod{});                                                     \
  }

#define REGISTER_STANDARD_DERIVATIVE(name, key, nGuards, type) \
  DEFINE_STANDARD_DERIV_CORE(name, key, nGuards, type)         \
  REGISTER_DERIVATIVE(name)                                    \
  BoutReal name::operator()(const stencil& f) const

#define REGISTER_UPWIND_DERIVATIVE(name, key, nGuards, type) \
  DEFINE_UPWIND_DERIV_CORE(name, key, nGuards, type)         \
  REGISTER_DERIVATIVE(name)                                  \
  BoutReal name::operator()(BoutReal vc, const stencil& f) const

#define REGISTER_FLUX_DERIVATIVE(name, key, nGuards, type) \
  DEFINE_FLUX_DERIV_CORE(name, key, nGuards, type)         \
  REGISTER_DERIVATIVE(name)                                \
  BoutReal name::operator()(const stencil& v, const stencil& f) const

#define REGISTER_STANDARD_STAGGERED_DERIVATIVE(name, key, nGuards, type) \
  DEFINE_STANDARD_DERIV_CORE(name, key, nGuards, type)                   \
  REGISTER_STAGGERED_DERIVATIVE(name)                                    \
  BoutReal name::operator()(const stencil& f) const

#define REGISTER_STANDARD_DERIVATIVE_STAGGERED(name, key, nGuards, type) \
  REGISTER_STANDARD_STAGGERED_DERIVATIVE(name, key, nGuards, type)

#define REGISTER_UPWIND_STAGGERED_DERIVATIVE(name, key, nGuards, type) \
  /*Note staggered upwind looks like flux*/                            \
  DEFINE_FLUX_DERIV_CORE(name, key, nGuards, type)                     \
  REGISTER_STAGGERED_DERIVATIVE(name)                                  \
  BoutReal name::operator()(const stencil& v, const stencil& f) const

#define REGISTER_UPWIND_DERIVATIVE_STAGGERED(name, key, nGuards, type) \
  REGISTER_UPWIND_STAGGERED_DERIVATIVE(name, key, nGuards, type)

#define REGISTER_FLUX_STAGGERED_DERIVATIVE(name, key, nGuards, type) \
  DEFINE_FLUX_DERIV_CORE(name, key, nGuards, type)                   \
  REGISTER_STAGGERED_DERIVATIVE(name)                                \
  BoutReal name::operator()(const stencil& v, const stencil& f) const

#define REGISTER_FLUX_DERIVATIVE_STAGGERED(name, key, nGuards, type) \
  REGISTER_FLUX_STAGGERED_DERIVATIVE(name, key, nGuards, type)

#endif
