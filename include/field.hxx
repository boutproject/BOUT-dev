/*!
 * \file field.hxx
 * \brief field base class definition for differencing methods
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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
 */

class Field;

#ifndef FIELD_H
#define FIELD_H

#include "bout/build_config.hxx"

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

#include "field_data.hxx"

#include "bout/region.hxx"
#include "bout_types.hxx"
#include "boutcomm.hxx"
#include "boutexception.hxx"
#include <globals.hxx>
#include "msg_stack.hxx"
#include "bout/region.hxx"
#include "stencils.hxx"
#include "utils.hxx"
#include <bout/rvec.hxx>
#include "bout/traits.hxx"

#include "unused.hxx"

class Mesh;

/// Base class for scalar fields
class Field : public FieldData {
public:
  Field() = default;
  Field(const Field& other) = default;
  Field(Field&& other) = default;
  Field& operator=(const Field& other) = default;
  Field& operator=(Field&& other) = default;
  virtual ~Field() = default;

  Field(Mesh* localmesh, CELL_LOC location_in, DirectionTypes directions_in);

  /// Getters for DIRECTION types
  DirectionTypes getDirections() const {
    return directions;
  }
  YDirectionType getDirectionY() const {
    return directions.y;
  }
  ZDirectionType getDirectionZ() const {
    return directions.z;
  }

  /// Setters for *DirectionType
  virtual Field& setDirections(DirectionTypes directions_in) {
    directions = directions_in;
    return *this;
  }
  virtual Field& setDirectionY(YDirectionType y_type) {
    directions.y = y_type;
    return *this;
  }
  virtual Field& setDirectionZ(ZDirectionType z_type) {
    directions.z = z_type;
    return *this;
  }

  std::string name;

#if CHECK > 0
  // Routines to test guard/boundary cells set
  
  virtual bool bndryValid() {
    if(!bndry_xin)
      throw BoutException("Inner X guard cells not set\n");
    if(!bndry_xout)
      throw BoutException("Outer X guard cells not set\n");
    if(!bndry_yup)
      throw BoutException("Upper y guard cells not set\n");
    if(!bndry_ydown)
      throw BoutException("Lower y guard cells not set\n");
    return true;
  }
  
  /// Status of the 4 boundaries
  bool bndry_xin{true}, bndry_xout{true}, bndry_yup{true}, bndry_ydown{true};
#endif

  /*!
   * Return the number of nx points
   */
  virtual int getNx() const;
  /*!
   * Return the number of ny points
   */
  virtual int getNy() const;
  /*!
   * Return the number of nz points
   */
  virtual int getNz() const;

  friend void swap(Field& first, Field& second) noexcept {
    using std::swap;
    swap(static_cast<FieldData&>(first), static_cast<FieldData&>(second));
    swap(first.name, second.name);
    swap(first.directions, second.directions);
  }

private:
  /// Labels for the type of coordinate system this field is defined over
  DirectionTypes directions{YDirectionType::Standard, ZDirectionType::Standard};
};

/// Check if Fields have compatible meta-data
inline bool areFieldsCompatible(const Field& field1, const Field& field2) {
  return
      field1.getCoordinates() == field2.getCoordinates() &&
      field1.getMesh() == field2.getMesh() &&
      field1.getLocation() == field2.getLocation() &&
      areDirectionsCompatible(field1.getDirections(), field2.getDirections());
}

#if CHECKLEVEL >= 1
#define ASSERT1_FIELDS_COMPATIBLE(field1, field2)			\
  if ((field1).getLocation() != (field2).getLocation()){		\
    throw BoutException("Error in {:s}:{:d}\nFields at different position:" \
			"`{:s}` at {:s}, `{:s}` at {:s}",__FILE__,__LINE__, \
			#field1, toString((field1).getLocation()),	\
			#field2, toString((field2).getLocation()));	\
  }									\
  if ((field1).getCoordinates() != (field2).getCoordinates()){		\
    throw BoutException("Error in {:s}:{:d}\nFields have different coordinates:" \
			"`{:s}` at {:p}, `{:s}` at {:p}",__FILE__,__LINE__, \
			#field1, static_cast<void*>((field1).getCoordinates()), \
			#field2, static_cast<void*>((field2).getCoordinates())); \
  }								\
  if ((field1).getMesh() != (field2).getMesh()){			\
    throw BoutException("Error in {:s}:{:d}\nFields are on different Meshes:" \
			"`{:s}` at {:p}, `{:s}` at {:p}",__FILE__,__LINE__, \
			#field1, static_cast<void*>((field1).getMesh()), \
			#field2, static_cast<void*>((field2).getMesh())); \
  }									\
  if (!areDirectionsCompatible((field1).getDirections(),		\
			       (field2).getDirections())){		\
    throw BoutException("Error in {:s}:{:d}\nFields at different directions:" \
			"`{:s}` at {:s}, `{:s}` at {:s}",__FILE__,__LINE__, \
			#field1, toString((field1).getDirections()),	\
			#field2, toString((field2).getDirections()));	\
  }

#else
#define ASSERT1_FIELDS_COMPATIBLE(field1, field2);
#endif

/// Return an empty shell field of some type derived from Field, with metadata
/// copied and a data array that is allocated but not initialised.
template<typename T>
inline T emptyFrom(const T& f) {
  static_assert(bout::utils::is_Field<T>::value, "emptyFrom only works on Fields");
  return T(f.getMesh(), f.getLocation(), {f.getDirectionY(), f.getDirectionZ()}).allocate();
}

/// Return a field of some type derived from Field, with metadata copied from
/// another field and a data array allocated and initialised to zero.
template<typename T>
inline T zeroFrom(const T& f) {
  static_assert(bout::utils::is_Field<T>::value, "zeroFrom only works on Fields");
  T result{emptyFrom(f)};
  result = 0.;
  return result;
}

/// Return a field of some type derived from Field, with metadata copied from
/// another field and a data array allocated and filled with the given value.
template<typename T>
inline T filledFrom(const T& f, BoutReal fill_value) {
  static_assert(bout::utils::is_Field<T>::value, "filledFrom only works on Fields");
  T result{emptyFrom(f)};
  result = fill_value;
  return result;
}

/// Return a field of some type derived from Field, with metadata copied from
/// another field and a data array allocated and filled using a callable e.g. lambda function
///
/// e.g.
///   Field3D result = filledFrom(field, [&](const auto& index) {
///                                          return ...;
///                                      });
/// 
/// An optional third argument is the region string
template <
    typename T, typename Function,
    typename = decltype(std::declval<Function&>()(std::declval<typename T::ind_type&>()))>
inline T filledFrom(const T& f, Function func, std::string region_string = "RGN_ALL") {
  static_assert(bout::utils::is_Field<T>::value, "filledFrom only works on Fields");
  T result{emptyFrom(f)};
  BOUT_FOR(i, result.getRegion(region_string)) {
    result[i] = func(i);
  }
  return result;
}

/// Unary + operator. This doesn't do anything
template<typename T, typename = bout::utils::EnableIfField<T>>
T operator+(const T& f) {return f;}

namespace bout {
/// Check if all values of a field \p var are finite.  Loops over all points including the
/// boundaries by default (can be changed using the \p rgn argument)
/// If any element is not finite, throws an exception that includes the position of the
/// first found.
///
/// Note that checkFinite runs the check irrespective of CHECK level. It is intended to be
/// used during initialization, where we always want to check inputs, even for optimized
/// builds.
template<typename T>
inline void checkFinite(const T& f, const std::string& name="field", const std::string& rgn="RGN_ALL") {
  AUTO_TRACE();

  if (!f.isAllocated()) {
    throw BoutException("{:s} is not allocated", name);
  }

  BOUT_FOR_SERIAL(i, f.getRegion(rgn)) {
    if (!std::isfinite(f[i])) {
      throw BoutException("{:s} is not finite at {:s}", name, toString(i));
    }
  }
}

/// Check if all values of a field \p var are positive.  Loops over all points including
/// the boundaries by default (can be changed using the \p rgn argument)
/// If any element is not finite, throws an exception that includes the position of the
/// first found.
///
/// Note that checkPositive runs the check irrespective of CHECK level. It is intended to
/// be used during initialization, where we always want to check inputs, even for
/// optimized builds.
template<typename T>
inline void checkPositive(const T& f, const std::string& name="field", const std::string& rgn="RGN_ALL") {
  AUTO_TRACE();

  if (!f.isAllocated()) {
    throw BoutException("{:s} is not allocated", name);
  }

  BOUT_FOR_SERIAL(i, f.getRegion(rgn)) {
    if (f[i] <= 0.) {
      throw BoutException("{:s} ({:s} {:s}) is {:e} (not positive) at {:s}", name,
                          toString(f.getLocation()), toString(f.getDirections()), f[i],
                          toString(i));
    }
  }
}
} // namespace bout

//////////////// NON-MEMBER FUNCTIONS //////////////////

/// Convert \p f to field-aligned space in \p region (default: whole domain)
template <typename T>
inline T toFieldAligned(const T& f, const std::string& region = "RGN_ALL") {
  static_assert(bout::utils::is_Field<T>::value, "toFieldAligned only works on Fields");
  return f.getCoordinates()->getParallelTransform().toFieldAligned(f, region);
}

template<typename T>
[[deprecated("Please use toFieldAligned(const T& f, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline T toFieldAligned(const T& f, REGION region) {
  return toFieldAligned(f, toString(region));
}

/// Convert \p f from field-aligned space in \p region (default: whole domain)
template <typename T>
inline T fromFieldAligned(const T& f, const std::string& region = "RGN_ALL") {
  static_assert(bout::utils::is_Field<T>::value, "fromFieldAligned only works on Fields");
  return f.getCoordinates()->getParallelTransform().fromFieldAligned(f, region);
}

template<typename T>
[[deprecated("Please use fromFieldAligned(const T& f, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline T fromFieldAligned(const T& f, REGION region) {
  return fromFieldAligned(f, toString(region));
}

/// Minimum of \p f, excluding the boundary/guard cells by default
/// (can be changed with \p rgn argument).
///
/// By default this is only on the local processor, but setting \p
/// allpe true does a collective Allreduce over all processors.
///
/// @param[in] f      Input field
/// @param[in] allpe  Minimum over all processors?
/// @param[in] rgn    The region to calculate the result over
template <typename T, typename = bout::utils::EnableIfField<T>>
inline BoutReal min(const T& f, bool allpe = false,
                    const std::string& rgn = "RGN_NOBNDRY") {
  AUTO_TRACE();

  checkData(f);

  const auto region = f.getRegion(rgn);
  BoutReal result = f[*region.cbegin()];

  BOUT_FOR_OMP(i, region, parallel for reduction(min:result)) {
    if (f[i] < result) {
      result = f[i];
    }
  }

  if (allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }

  return result;
}

template<typename T, typename = bout::utils::EnableIfField<T>>
[[deprecated("Please use Field3D min(const Field3D& f, bool allpe, "
    "const std::string& region = \"RGN_NOBNDRY\") instead")]]
inline BoutReal min(const T& f, bool allpe, REGION rgn) {
  return min(f, allpe, toString(rgn));
}

/// Returns true if all elements of \p f over \p region are equal. By
/// default only checks the local processor, use \p allpe to check
/// globally
///
/// @param[in] f       The field to check
/// @param[in] allpe   Check over all processors
/// @param[in] region  The region to check for uniformity over
template <typename T, typename = bout::utils::EnableIfField<T>>
inline bool isUniform(const T& f, bool allpe = false,
                      const std::string& region = "RGN_ALL") {
  bool result = true;
  auto element = f[*f.getRegion(region).begin()];
  // TODO: maybe parallise this loop, as the early return is unlikely
  BOUT_FOR_SERIAL(i, f.getRegion(region)) {
    if (f[i] != element) {
      result = false;
      break;
    }
  }
  if (allpe) {
    bool localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_C_BOOL, MPI_LOR, BoutComm::get());
  }
  return result;
}

/// Returns the value of the first element of \p f (in the region \p
/// region if given). If checks are enabled, then throws an exception
/// if \p f is not uniform over \p region. By default only checks the
/// local processor, use \p allpe to check globally
///
/// @param[in] f       The field to check
/// @param[in] allpe   Check over all processors
/// @param[in] region  The region to assume is uniform
template <typename T, typename = bout::utils::EnableIfField<T>>
inline BoutReal getUniform(const T& f, bool allpe = false,
                           const std::string& region = "RGN_ALL") {
#if CHECK > 1
  if (not isUniform(f, allpe, region)) {
    throw BoutException("Requested getUniform({}, {}, {}) but Field is not const", f.name,
                        allpe, region);
  }
#endif
  return f[*f.getRegion(region).begin()];
}

/// Maximum of \p r, excluding the boundary/guard cells by default
/// (can be changed with \p rgn argument).
///
/// By default this is only on the local processor, but setting \p
/// allpe to true does a collective Allreduce over all processors.
///
/// @param[in] f      Input field
/// @param[in] allpe  Maximum over all processors?
/// @param[in] rgn    The region to calculate the result over
template <typename T, typename = bout::utils::EnableIfField<T>>
inline BoutReal max(const T& f, bool allpe = false,
                    const std::string& rgn = "RGN_NOBNDRY") {
  AUTO_TRACE();

  checkData(f);

  const auto region = f.getRegion(rgn);
  BoutReal result = f[*region.cbegin()];

  BOUT_FOR_OMP(i, region, parallel for reduction(max:result)) {
    if (f[i] > result) {
      result = f[i];
    }
  }

  if (allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
  }

  return result;
}

template<typename T, typename = bout::utils::EnableIfField<T>>
[[deprecated("Please use Field3D max(const Field3D& f, bool allpe, "
    "const std::string& region = \"RGN_NOBNDRY\") instead")]]
inline BoutReal max(const T& f, bool allpe, REGION rgn) {
  return max(f, allpe, toString(rgn));
}

/// Mean of \p f, excluding the boundary/guard cells by default (can
/// be changed with \p rgn argument).
///
/// By default this is only on the local processor, but setting \p
/// allpe to true does a collective Allreduce over all processors.
///
/// @param[in] f      Input field
/// @param[in] allpe  Mean over all processors?
/// @param[in] rgn    The region to calculate the result over
template <typename T, typename = bout::utils::EnableIfField<T>>
inline BoutReal mean(const T& f, bool allpe = false,
                     const std::string& rgn = "RGN_NOBNDRY") {
  AUTO_TRACE();

  checkData(f);

  // Intitialise the cummulative sum and counter
  BoutReal result = 0.;
  int count = 0;

  BOUT_FOR_OMP(i, f.getRegion(rgn), parallel for reduction(+:result,count)) {
    result += f[i];
    count += 1;
  }

  if (allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get());
    int localcount = count;
    MPI_Allreduce(&localcount, &count, 1, MPI_INT, MPI_SUM, BoutComm::get());
  }

  return result / static_cast<BoutReal>(count);
}

template<typename T, typename = bout::utils::EnableIfField<T>>
[[deprecated("Please use Field3D mean(const Field3D& f, bool allpe, "
    "const std::string& region = \"RGN_NOBNDRY\") instead")]]
inline BoutReal mean(const T& f, bool allpe, REGION rgn) {
  return mean(f, allpe, toString(rgn));
}

/// Exponent: pow(lhs, lhs) is \p lhs raised to the power of \p rhs
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument)
/// If CHECK >= 3 then the result will be checked for non-finite numbers
template<typename T, typename = bout::utils::EnableIfField<T>>
T pow(const T& lhs, const T& rhs, const std::string& rgn = "RGN_ALL") {
  AUTO_TRACE();

  ASSERT1(areFieldsCompatible(lhs, rhs));

  T result{emptyFrom(lhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs[i]); }

  checkData(result);
  return result;
}
template<typename T, typename = bout::utils::EnableIfField<T>>
[[deprecated("Please use pow(const T& lhs, const T& rhs"
    "const std::string& region = \"RGN_ALL\") instead")]]
inline T pow(const T& lhs, const T& rhs, REGION rgn) {
  return pow(lhs, rhs, toString(rgn));
}

template<typename T, typename = bout::utils::EnableIfField<T>>
T pow(const T &lhs, BoutReal rhs, const std::string& rgn = "RGN_ALL") {
  AUTO_TRACE();

  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  T result{emptyFrom(lhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs); }

  checkData(result);
  return result;
}
template<typename T, typename = bout::utils::EnableIfField<T>>
[[deprecated("Please use pow(const T& lhs, BoutReal rhs"
    "const std::string& region = \"RGN_ALL\") instead")]]
inline T pow(const T& lhs, BoutReal rhs, REGION rgn) {
  return pow(lhs, rhs, toString(rgn));
}

template<typename T, typename = bout::utils::EnableIfField<T>>
T pow(BoutReal lhs, const T &rhs, const std::string& rgn = "RGN_ALL") {
  AUTO_TRACE();

  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  // Define and allocate the output result
  T result{emptyFrom(rhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs, rhs[i]); }

  checkData(result);
  return result;
}
template<typename T, typename = bout::utils::EnableIfField<T>>
[[deprecated("Please use pow(BoutReal lhs, const T& rhs"
    "const std::string& region = \"RGN_ALL\") instead")]]
inline T pow(BoutReal lhs, const T& rhs, REGION rgn) {
  return pow(lhs, rhs, toString(rgn));
}


/*!
 * This macro takes a function \p func, which is
 * assumed to operate on a single BoutReal and return
 * a single BoutReal, and wraps it up into a function
 * of a Field called \p name.
 *
 * @param name  The name of the function to define
 * @param func  The function to apply to each value
 *
 * If CHECK >= 1, checks if the Field is allocated
 *
 * Loops over the entire domain, applies function,
 * and uses checkData() to, if CHECK >= 3, check
 * result for non-finite numbers
 *
 */
#ifdef FIELD_FUNC
#error This macro has already been defined
#else
#define FIELD_FUNC(name, func)                                                       \
  template<typename T, typename = bout::utils::EnableIfField<T>>                     \
  inline T name(const T &f, const std::string& rgn = "RGN_ALL") {                    \
    AUTO_TRACE();                                                                    \
    /* Check if the input is allocated */                                            \
    checkData(f);                                                                    \
    /* Define and allocate the output result */                                      \
    T result{emptyFrom(f)};                                                          \
    BOUT_FOR(d, result.getRegion(rgn)) { result[d] = func(f[d]); }                   \
    checkData(result);                                                               \
    return result;                                                                   \
  }                                                                                  \
  template<typename T, typename = bout::utils::EnableIfField<T>>                     \
  [[deprecated("Please use func(const T& f, "                                   \
      "const std::string& region = \"RGN_ALL\") instead")]]                          \
  inline T name(const T& f, REGION region) {                                         \
    return name(f, toString(region));                                                \
  }
#endif

/// Square root of \p f over region \p rgn
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(sqrt, ::sqrt)

/// Absolute value (modulus, |f|) of \p f over region \p rgn
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(abs, ::fabs)

/// Exponential: \f$\exp(f)\f$ is e to the power of \p f, over region
/// \p rgn
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(exp, ::exp)

/// Natural logarithm of \p f over region \p rgn, inverse of
/// exponential
///
///     \f$\ln(\exp(f)) = f\f$
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the rgn argument)
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(log, ::log)

/// Sine trigonometric function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(sin, ::sin)

/// Cosine trigonometric function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(cos, ::cos)

/// Tangent trigonometric function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(tan, ::tan)

/// Hyperbolic sine trigonometric function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(sinh, ::sinh)

/// Hyperbolic cosine trigonometric function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(cosh, ::cosh)

/// Hyperbolic tangent trigonometric function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
FIELD_FUNC(tanh, ::tanh)

/// Check if all values of a field \p var are finite.
/// Loops over all points including the boundaries by
/// default (can be changed using the \p rgn argument
template<typename T, typename = bout::utils::EnableIfField<T>>
inline bool finite(const T &f, const std::string& rgn = "RGN_ALL") {
  AUTO_TRACE();

  if (!f.isAllocated()) {
    return false;
  }

  BOUT_FOR_SERIAL(i, f.getRegion(rgn)) {
    if (!std::isfinite(f[i])) {
      return false;
    }
  }

  return true;
}
template<typename T, typename = bout::utils::EnableIfField<T>>
[[deprecated("Please use bool finite(const Field3D& f, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline bool finite(const T& f, REGION rgn) {
  return finite(f, toString(rgn));
}

/// Makes a copy of a field \p f, ensuring that the underlying data is
/// not shared.
template<typename T, typename = bout::utils::EnableIfField<T>>
T copy(const T &f) {
  T result = f;
  result.allocate();
  return result;
}

/// Apply a floor value \p f to a field \p var. Any value lower than
/// the floor is set to the floor.
///
/// @param[in] var  Variable to apply floor to
/// @param[in] f    The floor value
/// @param[in] rgn  The region to calculate the result over
template<typename T, typename = bout::utils::EnableIfField<T>>
inline T floor(const T& var, BoutReal f, const std::string& rgn = "RGN_ALL") {
  checkData(var);
  T result = copy(var);

  BOUT_FOR(d, var.getRegion(rgn)) {
    if (result[d] < f) {
      result[d] = f;
    }
  }

  return result;
}
template<typename T, typename = bout::utils::EnableIfField<T>>
[[deprecated("Please use floor(const T& var, BoutReal f, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline T floor(const T& var, BoutReal f, REGION rgn) {
  return floor(var, f, toString(rgn));
}

#undef FIELD_FUNC

#endif /* FIELD_H */
