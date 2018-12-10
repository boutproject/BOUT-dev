///
/// Variant utilities
///
/// All in namespace bout::utils
///    variant
///    visit
///    holds_alternative
///    get
///
///    variantEqualTo
///    variantStaticCastOrThrow
///    variantToString
///
/// Internal implementation in bout::utils::details

#pragma once

#ifndef __VARIANT_HXX__
#define __VARIANT_HXX__

#include <variant>

#include "utils.hxx"

namespace bout {
namespace utils {

/// Import variant, visit into bout::utils namespace
using std::variant;
using std::visit;
using std::holds_alternative;
using std::get;

////////////////////////////////////////////////////////////
// Variant comparison

namespace details {

/// Compare two values.
///   Different types -> false
template <typename T, typename U>
struct CompareTypes {
  bool operator()(const T& UNUSED(v), const U& UNUSED(t)) { return false; }
};

/// Compare two values
///   Same type -> use `==` operator to compare
template <typename T>
struct CompareTypes<T, T> {
  bool operator()(const T& v, const T& t) { return v == t; }
};

/// A visitor for std::variant which compares
/// the value stored in the variant with a given value using CompareTypes
template <typename T>
struct IsEqual {
  const T& t;
  IsEqual(const T& t) : t(t) {}

  template <typename U>
  bool operator()(const U& u) {
    return CompareTypes<T, U>()(t, u);
  }
};
} // namespace details

/// Return true only if the given variant \p v
/// has the same type and value as \p t
///
/// Note: Handles the case that \p t is not of a type
/// which \v can hold.
template <typename Variant, typename T>
bool variantEqualTo(const Variant& v, const T& t) {
  return visit(details::IsEqual(t), v);
}

////////////////////////////////////////////////////////////
// Variant casting

namespace details {

/// Helper class for casting between two types
/// The general case is that casting can't be done
///
/// Note: This must be a class/struct because partial
/// specialisation is not allowed for functions
template <typename Target, typename Source, bool>
struct _StaticCastOrThrow {
  Target operator()(Source&& UNUSED(source)) { throw std::bad_cast{}; }
};

/// Specialised case (bool = true)
/// This version is used when casting can be done
template <typename Target, typename Source>
struct _StaticCastOrThrow<Target, Source, true> {
  Target operator()(Source&& source) {
    return static_cast<Target>(std::forward<Source>(source));
  }
};

/// Functor to perform static casting with std::visit
/// If the Target cannot be constructed from the Source
/// then an exception (std::bad_cast) will be thrown at run time.
///
/// Note: This needs to be at runtime because the particular
/// type which a variant is holding is only known at runtime.
template <typename Target>
struct StaticCastOrThrow {
  template <typename Source>
  Target operator()(Source&& source) const {
    return _StaticCastOrThrow<Target, Source,
                              std::is_constructible<Target, Source>::value>()(
        std::forward<Source>(source));
  }
};
} // namespace details

/// Cast a variant to a given type using static_cast
/// If this can't be done then a std::bad_cast exception is thrown
///
/// Note: \p T can be a type which variant \v cannot hold
/// in which case std::bad_cast will be thrown at runtime
template <typename Variant, typename T>
T variantStaticCastOrThrow(const Variant &v) {
  std::visit( details::StaticCastOrThrow<T>(), v );
}

namespace details {
  
struct ToString {
  template <typename T>
  std::string operator()(T&& val) {
    return toString(std::forward<T>(val));
  }
};
  
} // namespace details

template <typename Variant>
std::string variantToString(const Variant& v) { return visit(details::ToString(), v); }

} // namespace utils
} // namespace bout

#endif //__VARIANT_HXX__
