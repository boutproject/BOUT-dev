#ifndef BOUT_TRAITS_H
#define BOUT_TRAITS_H

#include <type_traits>

class Field;
class Field2D;
class Field3D;
class FieldPerp;
class Options;

namespace bout {
namespace utils {

namespace details {
/// Helper class for fold expressions pre-C++17
///
/// Taken from "C++ Templates: The Complete Guide, Second Edition"
///  Addison-Wesley, 2017
///  ISBN-13:  978-0-321-71412-1
///  ISBN-10:      0-321-71412-1
/// Copyright © 2017 by Addison-Wesley, David Vandevoorde, Nicolai
/// M. Josuttis, and Douglas Gregor.
constexpr bool and_all() { return true; }
template <class T>
constexpr bool and_all(T cond) {
  return cond;
}
template <class T, class... Ts>
constexpr bool and_all(T cond, Ts... conds) {
  return cond and and_all(conds...);
}
} // namespace details

/// If `T` is derived from `Field`, provides the member constant
/// `value` equal to `true`. Otherwise `value is `false`.
///
/// The following is C++14, but simplifies the use of `is_field`:
///
///     template <class T>
///     constexpr bool is_field_v = is_field<T>::value;
///
/// Examples
/// --------
///
///     template <class T>
///     void print_field(const T& field) {
///       static_assert(bout::utils::is_field<T>::value,
///           "print_field only works with Field2Ds, Field3Ds or FieldPerps")
///       // implementation
///     }
template <class T>
using is_Field = std::is_base_of<Field, T>;

/// If `T` is derived from `Field2D`, provides the member constant
/// `value` equal to `true`. Otherwise `value is `false`.
template <class T>
using is_Field2D = std::is_base_of<Field2D, T>;

/// If `T` is derived from `Field3D`, provides the member constant
/// `value` equal to `true`. Otherwise `value is `false`.
template <class T>
using is_Field3D = std::is_base_of<Field3D, T>;

/// If `T` is derived from `FieldPerp`, provides the member constant
/// `value` equal to `true`. Otherwise `value is `false`.
template <class T>
using is_FieldPerp = std::is_base_of<FieldPerp, T>;

/// If `T` is derived from `Options`, provides the member constant
/// `value` equal to `true`. Otherwise `value is `false`.
template <class T>
using is_Options = std::is_base_of<Options, T>;

/// Enable a function if all the Ts are subclasses of `Field`, and
/// returns the common type: i.e. `Field3D` if at least one argument
/// is `Field3D`, otherwise `Field2D` if they are all `Field2D`
///
/// This is most useful in two particular cases:
/// 1. when there are multiple overloads for a function but some only
///    make sense for fields (as opposed to `BoutReal`, say) or
///    vice-versa, and some overloads should not be used for fields
/// 2. when a function takes multiple fields and the return type is
///    also a field and must be "big enough"
///
/// In other cases, such as a function without overloads that only
/// works for fields, consider using `static_assert` with `is_Field`
/// to give a nice compile-time error
///
/// Examples
/// --------
///
/// Consider the following template function:
///
///     template <class T, class U, class V,
///          class ResultType = typename bout::utils::EnableIfField<T, U, V>>
///     auto where(const T& test, const U& gt0, const V& le0) -> ResultType {
///       // function body
///     }
///
/// This function only "appears" if `T`, `U` and `V` are all
/// subclasses of `Field`. `ResultType` is the common type of `T`, `U`
/// and `V`. If `T` and `U` are both `Field2D`, `ResultType` is
/// `Field2D` if `V` is `Field2D`, and `Field3D` if `V` is `Field3D`.
template <class... Ts>
using EnableIfField =
    typename std::enable_if<details::and_all(is_Field<Ts>::value ...),
                            typename std::common_type<Ts...>::type>::type;

/// Enable a function if all the Ts are subclasses of `Field2D`, and
/// returns the common type
template <class... Ts>
using EnableIfField2D =
    typename std::enable_if<details::and_all(is_Field2D<Ts>::value ...),
                            typename std::common_type<Ts...>::type>::type;

/// Enable a function if all the Ts are subclasses of `Field3D`, and
/// returns the common type
template <class... Ts>
using EnableIfField3D =
    typename std::enable_if<details::and_all(is_Field3D<Ts>::value ...),
                            typename std::common_type<Ts...>::type>::type;

/// Enable a function if all the Ts are subclasses of `FieldPerp`, and
/// returns the common type
template <class... Ts>
using EnableIfFieldPerp =
    typename std::enable_if<details::and_all(is_FieldPerp<Ts>::value ...),
                            typename std::common_type<Ts...>::type>::type;

/// Enable a function if T is a subclass of Options
template <class T>
using EnableIfOptions = std::enable_if_t<std::is_base_of<Options, T>::value>;
} // namespace utils
} // namespace bout

#endif // BOUT_TRAITS_H
