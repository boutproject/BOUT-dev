/*!
 * \file field.hxx
 * \brief field base class definition for differencing methods
 *
 **************************************************************************
 * Copyright 2010 - 2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#include <cmath>
#include <cstdio>
#include <optional>
#include <string>
#include <type_traits>

#include "bout/bout_types.hxx"
#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/build_config.hxx"
#include "bout/field_data.hxx"
#include "bout/region.hxx"
#include "bout/traits.hxx"
#include "bout/utils.hxx"
#include <bout/globals.hxx>
#include <bout/rvec.hxx>

#include "bout/fieldops.hxx"

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
  DirectionTypes getDirections() const { return directions; }
  YDirectionType getDirectionY() const { return directions.y; }
  ZDirectionType getDirectionZ() const { return directions.z; }

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

  bool isFci() const;

#if CHECK > 0
  // Routines to test guard/boundary cells set

  virtual bool bndryValid() {
    if (!bndry_xin) {
      throw BoutException("Inner X guard cells not set\n");
    }
    if (!bndry_xout) {
      throw BoutException("Outer X guard cells not set\n");
    }
    if (!bndry_yup) {
      throw BoutException("Upper y guard cells not set\n");
    }
    if (!bndry_ydown) {
      throw BoutException("Lower y guard cells not set\n");
    }
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

  /// Get the total number of points
  virtual int size() const = 0;

  friend void swap(Field& first, Field& second) noexcept {
    using std::swap;
    swap(static_cast<FieldData&>(first), static_cast<FieldData&>(second));
    swap(first.name, second.name);
    swap(first.directions, second.directions);
  }

  /// Dummy functions to increase portability
  virtual void setRegion([[maybe_unused]] size_t regionID) {}
  virtual void setRegion([[maybe_unused]] std::optional<size_t> regionID) {}
  virtual void setRegion([[maybe_unused]] const std::string& region_name) {}
  virtual void resetRegion() {}
  virtual std::optional<size_t> getRegionID() const { return {}; }
  virtual bool hasParallelSlices() const { return true; }
  virtual void calcParallelSlices() {}
  virtual void splitParallelSlices() {}
  virtual void clearParallelSlices() {}
  virtual size_t numberParallelSlices() const { return 0; }

private:
  /// Labels for the type of coordinate system this field is defined over
  DirectionTypes directions{YDirectionType::Standard, ZDirectionType::Standard};
};

/// Check if Fields have compatible meta-data
inline bool areFieldsCompatible(const Field& field1, const Field& field2) {
  return field1.getCoordinates() == field2.getCoordinates()
         && field1.getMesh() == field2.getMesh()
         && field1.getLocation() == field2.getLocation()
         && areDirectionsCompatible(field1.getDirections(), field2.getDirections());
}

#if CHECKLEVEL >= 1
#define ASSERT1_FIELDS_COMPATIBLE(field1, field2)                                        \
  if ((field1).getLocation() != (field2).getLocation()) {                                \
    throw BoutException("Error in {:s}:{:d}\nFields at different position:"              \
                        "`{:s}` at {:s}, `{:s}` at {:s}",                                \
                        __FILE__, __LINE__, #field1, toString((field1).getLocation()),   \
                        #field2, toString((field2).getLocation()));                      \
  }                                                                                      \
  if ((field1).getCoordinates() != (field2).getCoordinates()) {                          \
    throw BoutException("Error in {:s}:{:d}\nFields have different coordinates:"         \
                        "`{:s}` at {:p}, `{:s}` at {:p}",                                \
                        __FILE__, __LINE__, #field1,                                     \
                        static_cast<void*>((field1).getCoordinates()), #field2,          \
                        static_cast<void*>((field2).getCoordinates()));                  \
  }                                                                                      \
  if ((field1).getMesh() != (field2).getMesh()) {                                        \
    throw BoutException("Error in {:s}:{:d}\nFields are on different Meshes:"            \
                        "`{:s}` at {:p}, `{:s}` at {:p}",                                \
                        __FILE__, __LINE__, #field1,                                     \
                        static_cast<void*>((field1).getMesh()), #field2,                 \
                        static_cast<void*>((field2).getMesh()));                         \
  }                                                                                      \
  if (!areDirectionsCompatible((field1).getDirections(), (field2).getDirections())) {    \
    throw BoutException("Error in {:s}:{:d}\nFields at different directions:"            \
                        "`{:s}` at {:s}, `{:s}` at {:s}",                                \
                        __FILE__, __LINE__, #field1, toString((field1).getDirections()), \
                        #field2, toString((field2).getDirections()));                    \
  }

#define ASSERT1_EXPR_COMPATIBLE(expr1, expr2)                                          \
  if ((expr1).getLocation() != (expr2).getLocation()) {                                \
    throw BoutException("Error in {:s}:{:d}\nFields at different position:"            \
                        "`{:s}` at {:s}, `{:s}` at {:s}",                              \
                        __FILE__, __LINE__, #expr1, toString((expr1).getLocation()),   \
                        #expr2, toString((expr2).getLocation()));                      \
  }                                                                                    \
  if ((expr1).getMesh() != (expr2).getMesh()) {                                        \
    throw BoutException("Error in {:s}:{:d}\nFields are on different Meshes:"          \
                        "`{:s}` at {:p}, `{:s}` at {:p}",                              \
                        __FILE__, __LINE__, #expr1,                                    \
                        static_cast<void*>((expr1).getMesh()), #expr2,                 \
                        static_cast<void*>((expr2).getMesh()));                        \
  }                                                                                    \
  if (!areDirectionsCompatible((expr1).getDirections(), (expr2).getDirections())) {    \
    throw BoutException("Error in {:s}:{:d}\nFields at different directions:"          \
                        "`{:s}` at {:s}, `{:s}` at {:s}",                              \
                        __FILE__, __LINE__, #expr1, toString((expr1).getDirections()), \
                        #expr2, toString((expr2).getDirections()));                    \
  }

#else
#define ASSERT1_FIELDS_COMPATIBLE(field1, field2) ;
#define ASSERT1_EXPR_COMPATIBLE(expr1, expr2) ;
#endif

/// Return an empty shell field of some type derived from Field, with metadata
/// copied and a data array that is allocated but not initialised.
template <typename T>
inline T emptyFrom(const T& f) {
  static_assert(bout::utils::is_Field_v<T>, "emptyFrom only works on Fields");
  return T(f.getMesh(), f.getLocation(),
           DirectionTypes{f.getDirectionY(), f.getDirectionZ()}, f.getRegionID())
      .allocate();
}

/// Return a field of some type derived from Field, with metadata copied from
/// another field and a data array allocated and initialised to zero.
template <typename T>
inline T zeroFrom(const T& f) {
  static_assert(bout::utils::is_Field_v<T>, "zeroFrom only works on Fields");
  T result{emptyFrom(f)};
  result = 0.;
  return result;
}

/// Return a field of some type derived from Field, with metadata copied from
/// another field and a data array allocated and filled with the given value.
template <typename T>
inline T filledFrom(const T& f, BoutReal fill_value) {
  static_assert(bout::utils::is_Field_v<T>, "filledFrom only works on Fields");
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
  static_assert(bout::utils::is_Field_v<T>, "filledFrom only works on Fields");
  T result{emptyFrom(f)};
  BOUT_FOR(i, result.getRegion(region_string)) { result[i] = func(i); }
  return result;
}

/// Unary + operator. This doesn't do anything
template <typename T, typename = bout::utils::EnableIfField<T>>
T operator+(const T& f) {
  return f;
}

namespace bout {
/// Check if all values of a field \p var are finite.  Loops over all points including the
/// boundaries by default (can be changed using the \p rgn argument)
/// If any element is not finite, throws an exception that includes the position of the
/// first found.
///
/// Note that checkFinite runs the check irrespective of CHECK level. It is intended to be
/// used during initialization, where we always want to check inputs, even for optimized
/// builds.
template <typename T>
inline void checkFinite(const T& f, const std::string& name = "field",
                        const std::string& rgn = "RGN_ALL") {

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
template <typename T>
inline void checkPositive(const T& f, const std::string& name = "field",
                          const std::string& rgn = "RGN_ALL") {

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
  static_assert(bout::utils::is_Field_v<T>, "toFieldAligned only works on Fields");
  ASSERT3(f.getCoordinates() != nullptr);
  return f.getCoordinates()->getParallelTransform().toFieldAligned(f, region);
}

/// Convert \p f from field-aligned space in \p region (default: whole domain)
template <typename T>
inline T fromFieldAligned(const T& f, const std::string& region = "RGN_ALL") {
  static_assert(bout::utils::is_Field_v<T>, "fromFieldAligned only works on Fields");
  ASSERT3(f.getCoordinates() != nullptr);
  return f.getCoordinates()->getParallelTransform().fromFieldAligned(f, region);
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

template <typename ResT, typename L, typename R, typename Func>
inline BoutReal min(const BinaryExpr<ResT, L, R, Func>& f, bool allpe = false,
                    const std::string& rgn = "RGN_NOBNDRY") {
  const auto& region = f.getMesh()->template getRegion<ResT>(rgn);
  const auto reduction_view =
      makeReductionView(static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),
                        region.getLinearIndices());
  BoutReal result =
      bout::reduce::Min::finalize(reduceExpr<bout::reduce::Min>(reduction_view));

  if (allpe) {
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }

  return result;
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
    // by default we only check for exact equality, as that should cover most cases
    if (f[i] != element and (not almost_equal(f[i], element, 10))) {
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
inline BoutReal getUniform(const T& f, [[maybe_unused]] bool allpe = false,
                           const std::string& region = "RGN_ALL") {
#if CHECK > 1
  if (not isUniform(f, allpe, region)) {
    const BoutReal f1 = min(f, allpe, region);
    const BoutReal f2 = max(f, allpe, region);
    throw BoutException("Requested getUniform({}, {}, {}) but Field is not const "
                        "([{:.15f}...{:.15f}] Δ={:e} {:e}ε)",
                        f.name, allpe, region, f1, f2, f2 - f1,
                        (f2 - f1) / (f1 + f2) / std::numeric_limits<BoutReal>::epsilon());
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

template <typename ResT, typename L, typename R, typename Func>
inline BoutReal max(const BinaryExpr<ResT, L, R, Func>& f, bool allpe = false,
                    const std::string& rgn = "RGN_NOBNDRY") {
  const auto& region = f.getMesh()->template getRegion<ResT>(rgn);
  const auto reduction_view =
      makeReductionView(static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),
                        region.getLinearIndices());
  BoutReal result =
      bout::reduce::Max::finalize(reduceExpr<bout::reduce::Max>(reduction_view));

  if (allpe) {
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
  }

  return result;
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

template <typename ResT, typename L, typename R, typename Func>
inline BoutReal mean(const BinaryExpr<ResT, L, R, Func>& f, bool allpe = false,
                     const std::string& rgn = "RGN_NOBNDRY") {
  const auto& region = f.getMesh()->template getRegion<ResT>(rgn);
  const auto reduction_view =
      makeReductionView(static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),
                        region.getLinearIndices());
  auto state = reduceExpr<bout::reduce::Mean>(reduction_view);

  if (allpe) {
    BoutReal localsum = state.sum;
    int localcount = state.count;
    MPI_Allreduce(&localsum, &state.sum, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get());
    MPI_Allreduce(&localcount, &state.count, 1, MPI_INT, MPI_SUM, BoutComm::get());
  }

  return bout::reduce::Mean::finalize(state);
}

namespace bout::op {
struct Pow {
  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& L,
                                                        const RView& R) const {
    return ::pow(L(idx), R(idx));
  }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(BoutReal a, BoutReal b) const {
    return ::pow(a, b);
  }
};
}; // namespace bout::op

namespace bout::detail {
template <typename T>
std::optional<int> getPerpYIndex(const T& value);

template <typename ResT, typename L, typename R, typename Func>
std::optional<int> getPerpYIndex(const BinaryExpr<ResT, L, R, Func>& expr);

template <typename ResT>
std::optional<size_t> getPowRegionID(const Mesh* mesh, const std::string& region_name) {
  if constexpr (std::is_same_v<ResT, Field3D>) {
    return bout::detail::getField3DRegionID(mesh, region_name);
  } else {
    return std::nullopt;
  }
}

template <typename ResT, typename L, typename R, typename LView, typename RView,
          typename IndType>
auto makePowExpr(const LView& lhs_view, const RView& rhs_view, Mesh* mesh,
                 CELL_LOC location, DirectionTypes directions,
                 std::optional<size_t> regionID, const Region<IndType>& region,
                 std::optional<int> yindex = std::nullopt) {
  return BinaryExpr<ResT, L, R, bout::op::Pow>{lhs_view, rhs_view, bout::op::Pow{},
                                               mesh,     location, directions,
                                               regionID, region,   yindex};
}
} // namespace bout::detail

/// Exponent: pow(lhs, lhs) is \p lhs raised to the power of \p rhs
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument)
/// If CHECK >= 3 then the result will be checked for non-finite numbers
template <typename L, typename R>
std::enable_if_t<is_expr_field2d_v<L> && is_expr_field2d_v<R>,
                 BinaryExpr<Field2D, L, R, bout::op::Pow>>
pow(const L& lhs, const R& rhs) {
  ASSERT1_EXPR_COMPATIBLE(lhs, rhs);
  return bout::detail::makePowExpr<Field2D, L, R>(
      static_cast<typename L::View>(lhs), static_cast<typename R::View>(rhs),
      lhs.getMesh(), lhs.getLocation(), lhs.getDirections(), std::nullopt,
      lhs.getMesh()->getRegion2D("RGN_ALL"));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field2d_v<L> && is_expr_field2d_v<R>,
                 BinaryExpr<Field2D, L, R, bout::op::Pow>>
pow(const L& lhs, const R& rhs, const std::string& rgn) {
  ASSERT1_EXPR_COMPATIBLE(lhs, rhs);
  return bout::detail::makePowExpr<Field2D, L, R>(
      static_cast<typename L::View>(lhs), static_cast<typename R::View>(rhs),
      lhs.getMesh(), lhs.getLocation(), lhs.getDirections(), std::nullopt,
      lhs.getMesh()->getRegion2D(rgn));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field3d_v<L> && is_expr_field3d_v<R>,
                 BinaryExpr<Field3D, L, R, bout::op::Pow>>
pow(const L& lhs, const R& rhs) {
  ASSERT1_EXPR_COMPATIBLE(lhs, rhs);
  auto regionID = lhs.getMesh()->getCommonRegion(lhs.getRegionID(), rhs.getRegionID());
  return bout::detail::makePowExpr<Field3D, L, R>(
      static_cast<typename L::View>(lhs), static_cast<typename R::View>(rhs),
      lhs.getMesh(), lhs.getLocation(), lhs.getDirections(), regionID,
      (regionID.has_value() ? lhs.getMesh()->getRegion(regionID.value())
                            : lhs.getMesh()->getRegion("RGN_ALL")),
      bout::detail::getPerpYIndex(lhs));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field3d_v<L> && is_expr_field3d_v<R>,
                 BinaryExpr<Field3D, L, R, bout::op::Pow>>
pow(const L& lhs, const R& rhs, const std::string& rgn) {
  ASSERT1_EXPR_COMPATIBLE(lhs, rhs);
  return bout::detail::makePowExpr<Field3D, L, R>(
      static_cast<typename L::View>(lhs), static_cast<typename R::View>(rhs),
      lhs.getMesh(), lhs.getLocation(), lhs.getDirections(),
      bout::detail::getPowRegionID<Field3D>(lhs.getMesh(), rgn),
      lhs.getMesh()->getRegion(rgn), bout::detail::getPerpYIndex(lhs));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field3d_v<L> && is_expr_field2d_v<R>,
                 BinaryExpr<Field3D, L, R, bout::op::Pow>>
pow(const L& lhs, const R& rhs) {
  ASSERT1_EXPR_COMPATIBLE(lhs, rhs);
  int mesh_nz = lhs.getMesh()->LocalNz;
  return bout::detail::makePowExpr<Field3D, L, R>(
      static_cast<typename L::View>(lhs),
      static_cast<typename R::View>(rhs).setScale(1, mesh_nz), lhs.getMesh(),
      lhs.getLocation(), lhs.getDirections(), lhs.getRegionID(),
      lhs.getMesh()->getRegion("RGN_ALL"), bout::detail::getPerpYIndex(lhs));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field3d_v<L> && is_expr_field2d_v<R>,
                 BinaryExpr<Field3D, L, R, bout::op::Pow>>
pow(const L& lhs, const R& rhs, const std::string& rgn) {
  ASSERT1_EXPR_COMPATIBLE(lhs, rhs);
  int mesh_nz = lhs.getMesh()->LocalNz;
  return bout::detail::makePowExpr<Field3D, L, R>(
      static_cast<typename L::View>(lhs),
      static_cast<typename R::View>(rhs).setScale(1, mesh_nz), lhs.getMesh(),
      lhs.getLocation(), lhs.getDirections(),
      bout::detail::getPowRegionID<Field3D>(lhs.getMesh(), rgn),
      lhs.getMesh()->getRegion(rgn), bout::detail::getPerpYIndex(lhs));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field2d_v<L> && is_expr_field3d_v<R>,
                 BinaryExpr<Field3D, L, R, bout::op::Pow>>
pow(const L& lhs, const R& rhs) {
  ASSERT1_EXPR_COMPATIBLE(lhs, rhs);
  int mesh_nz = rhs.getMesh()->LocalNz;
  return bout::detail::makePowExpr<Field3D, L, R>(
      static_cast<typename L::View>(lhs).setScale(1, mesh_nz),
      static_cast<typename R::View>(rhs), rhs.getMesh(), rhs.getLocation(),
      rhs.getDirections(), rhs.getRegionID(), rhs.getMesh()->getRegion("RGN_ALL"),
      bout::detail::getPerpYIndex(rhs));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field2d_v<L> && is_expr_field3d_v<R>,
                 BinaryExpr<Field3D, L, R, bout::op::Pow>>
pow(const L& lhs, const R& rhs, const std::string& rgn) {
  ASSERT1_EXPR_COMPATIBLE(lhs, rhs);
  int mesh_nz = rhs.getMesh()->LocalNz;
  return bout::detail::makePowExpr<Field3D, L, R>(
      static_cast<typename L::View>(lhs).setScale(1, mesh_nz),
      static_cast<typename R::View>(rhs), rhs.getMesh(), rhs.getLocation(),
      rhs.getDirections(), bout::detail::getPowRegionID<Field3D>(rhs.getMesh(), rgn),
      rhs.getMesh()->getRegion(rgn), bout::detail::getPerpYIndex(rhs));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field2d_v<L> && is_expr_constant_v<R>,
                 BinaryExpr<Field2D, L, Constant<R>, bout::op::Pow>>
pow(const L& lhs, R rhs) {
  return bout::detail::makePowExpr<Field2D, L, Constant<R>>(
      static_cast<typename L::View>(lhs), static_cast<typename Constant<R>::View>(rhs),
      lhs.getMesh(), lhs.getLocation(), lhs.getDirections(), std::nullopt,
      lhs.getMesh()->getRegion2D("RGN_ALL"));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field2d_v<L> && is_expr_constant_v<R>,
                 BinaryExpr<Field2D, L, Constant<R>, bout::op::Pow>>
pow(const L& lhs, R rhs, const std::string& rgn) {
  return bout::detail::makePowExpr<Field2D, L, Constant<R>>(
      static_cast<typename L::View>(lhs), static_cast<typename Constant<R>::View>(rhs),
      lhs.getMesh(), lhs.getLocation(), lhs.getDirections(), std::nullopt,
      lhs.getMesh()->getRegion2D(rgn));
}

template <typename L, typename R>
std::enable_if_t<is_expr_constant_v<L> && is_expr_field2d_v<R>,
                 BinaryExpr<Field2D, Constant<L>, R, bout::op::Pow>>
pow(L lhs, const R& rhs) {
  return bout::detail::makePowExpr<Field2D, Constant<L>, R>(
      static_cast<typename Constant<L>::View>(lhs), static_cast<typename R::View>(rhs),
      rhs.getMesh(), rhs.getLocation(), rhs.getDirections(), std::nullopt,
      rhs.getMesh()->getRegion2D("RGN_ALL"));
}

template <typename L, typename R>
std::enable_if_t<is_expr_constant_v<L> && is_expr_field2d_v<R>,
                 BinaryExpr<Field2D, Constant<L>, R, bout::op::Pow>>
pow(L lhs, const R& rhs, const std::string& rgn) {
  return bout::detail::makePowExpr<Field2D, Constant<L>, R>(
      static_cast<typename Constant<L>::View>(lhs), static_cast<typename R::View>(rhs),
      rhs.getMesh(), rhs.getLocation(), rhs.getDirections(), std::nullopt,
      rhs.getMesh()->getRegion2D(rgn));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field3d_v<L> && is_expr_constant_v<R>,
                 BinaryExpr<Field3D, L, Constant<R>, bout::op::Pow>>
pow(const L& lhs, R rhs) {
  return bout::detail::makePowExpr<Field3D, L, Constant<R>>(
      static_cast<typename L::View>(lhs), static_cast<typename Constant<R>::View>(rhs),
      lhs.getMesh(), lhs.getLocation(), lhs.getDirections(), lhs.getRegionID(),
      lhs.getMesh()->getRegion("RGN_ALL"), bout::detail::getPerpYIndex(lhs));
}

template <typename L, typename R>
std::enable_if_t<is_expr_field3d_v<L> && is_expr_constant_v<R>,
                 BinaryExpr<Field3D, L, Constant<R>, bout::op::Pow>>
pow(const L& lhs, R rhs, const std::string& rgn) {
  return bout::detail::makePowExpr<Field3D, L, Constant<R>>(
      static_cast<typename L::View>(lhs), static_cast<typename Constant<R>::View>(rhs),
      lhs.getMesh(), lhs.getLocation(), lhs.getDirections(),
      bout::detail::getPowRegionID<Field3D>(lhs.getMesh(), rgn),
      lhs.getMesh()->getRegion(rgn), bout::detail::getPerpYIndex(lhs));
}

template <typename L, typename R>
std::enable_if_t<is_expr_constant_v<L> && is_expr_field3d_v<R>,
                 BinaryExpr<Field3D, Constant<L>, R, bout::op::Pow>>
pow(L lhs, const R& rhs) {
  return bout::detail::makePowExpr<Field3D, Constant<L>, R>(
      static_cast<typename Constant<L>::View>(lhs), static_cast<typename R::View>(rhs),
      rhs.getMesh(), rhs.getLocation(), rhs.getDirections(), rhs.getRegionID(),
      rhs.getMesh()->getRegion("RGN_ALL"), bout::detail::getPerpYIndex(rhs));
}

template <typename L, typename R>
std::enable_if_t<is_expr_constant_v<L> && is_expr_field3d_v<R>,
                 BinaryExpr<Field3D, Constant<L>, R, bout::op::Pow>>
pow(L lhs, const R& rhs, const std::string& rgn) {
  return bout::detail::makePowExpr<Field3D, Constant<L>, R>(
      static_cast<typename Constant<L>::View>(lhs), static_cast<typename R::View>(rhs),
      rhs.getMesh(), rhs.getLocation(), rhs.getDirections(),
      bout::detail::getPowRegionID<Field3D>(rhs.getMesh(), rgn),
      rhs.getMesh()->getRegion(rgn), bout::detail::getPerpYIndex(rhs));
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
class Field3DParallel;
class FieldPerp;

namespace bout::detail {
template <typename T>
using UnaryFieldResult_t =
    std::conditional_t<std::is_same_v<std::decay_t<T>, ::Field3DParallel>, ::Field3D,
                       std::decay_t<T>>;

template <typename T>
std::optional<size_t> getUnaryRegionID(const Mesh* mesh, const std::string& region_name) {
  if constexpr (std::is_same_v<UnaryFieldResult_t<T>, ::Field3D>) {
    return bout::detail::getField3DRegionID(mesh, region_name);
  } else {
    return std::nullopt;
  }
}

template <typename T>
std::optional<int> getPerpYIndex(const T& value) {
  if constexpr (std::is_same_v<std::decay_t<T>, ::FieldPerp>) {
    return value.getIndex();
  } else {
    return std::nullopt;
  }
}

template <typename ResT, typename L, typename R, typename Func>
std::optional<int> getPerpYIndex(const BinaryExpr<ResT, L, R, Func>& expr) {
  if constexpr (std::is_same_v<ResT, ::FieldPerp>) {
    return expr.getIndex();
  } else {
    return std::nullopt;
  }
}
} // namespace bout::detail

#ifdef FIELD_FUNC
#error This macro has already been defined
#else
#define FIELD_FUNC(name, func)                                                          \
  namespace bout::op {                                                                  \
  struct name {                                                                         \
    template <typename LView, typename RView>                                           \
    BOUT_HOST_DEVICE BoutReal operator()(int idx, const LView& L, const RView&) const { \
      return func(L(idx));                                                              \
    }                                                                                   \
  };                                                                                    \
  };                                                                                    \
  template <typename T, typename = bout::utils::EnableIfField<T>>                       \
  inline auto name(const T& f, const std::string& rgn = "RGN_ALL") {                    \
    using ResT = bout::detail::UnaryFieldResult_t<T>;                                   \
    return BinaryExpr<ResT, T, T, bout::op::name>{                                      \
        static_cast<typename T::View>(f),                                               \
        static_cast<typename T::View>(f),                                               \
        bout::op::name{},                                                               \
        f.getMesh(),                                                                    \
        f.getLocation(),                                                                \
        f.getDirections(),                                                              \
        bout::detail::getUnaryRegionID<T>(f.getMesh(), rgn),                            \
        f.getMesh()->template getRegion<ResT>(rgn),                                     \
        bout::detail::getPerpYIndex(f)};                                                \
  }                                                                                     \
  template <typename ResT, typename L, typename R, typename Func>                       \
  inline auto name(const BinaryExpr<ResT, L, R, Func>& f) {                             \
    using UnaryResT = bout::detail::UnaryFieldResult_t<ResT>;                           \
    return BinaryExpr<UnaryResT, BinaryExpr<ResT, L, R, Func>,                          \
                      BinaryExpr<ResT, L, R, Func>, bout::op::name>{                    \
        static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),                    \
        static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),                    \
        bout::op::name{},                                                               \
        f.getMesh(),                                                                    \
        f.getLocation(),                                                                \
        f.getDirections(),                                                              \
        f.getRegionID(),                                                                \
        f.indices,                                                                      \
        bout::detail::getPerpYIndex(f)};                                                \
  }                                                                                     \
  template <typename ResT, typename L, typename R, typename Func>                       \
  inline auto name(const BinaryExpr<ResT, L, R, Func>& f, const std::string& rgn) {     \
    using UnaryResT = bout::detail::UnaryFieldResult_t<ResT>;                           \
    return BinaryExpr<UnaryResT, BinaryExpr<ResT, L, R, Func>,                          \
                      BinaryExpr<ResT, L, R, Func>, bout::op::name>{                    \
        static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),                    \
        static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),                    \
        bout::op::name{},                                                               \
        f.getMesh(),                                                                    \
        f.getLocation(),                                                                \
        f.getDirections(),                                                              \
        bout::detail::getUnaryRegionID<UnaryResT>(f.getMesh(), rgn),                    \
        f.getMesh()->template getRegion<UnaryResT>(rgn),                                \
        bout::detail::getPerpYIndex(f)};                                                \
  }
#endif

namespace bout::op {
struct Square {
  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BoutReal operator()(int idx, const LView& L, const RView&) const {
    const BoutReal value = L(idx);
    return ::SQ(value);
  }
};

struct Floor {
  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& L,
                                                        const RView& R) const {
    const BoutReal value = L(idx);
    const BoutReal floor_value = R(idx);
    return value < floor_value ? floor_value : value;
  }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(BoutReal value,
                                                        BoutReal floor_value) const {
    return value < floor_value ? floor_value : value;
  }
};
}; // namespace bout::op

template <typename T, typename = bout::utils::EnableIfField<T>>
inline auto SQ(const T& f, const std::string& rgn = "RGN_ALL") {
  using ResT = bout::detail::UnaryFieldResult_t<T>;
  return BinaryExpr<ResT, T, T, bout::op::Square>{
      static_cast<typename T::View>(f),
      static_cast<typename T::View>(f),
      bout::op::Square{},
      f.getMesh(),
      f.getLocation(),
      f.getDirections(),
      bout::detail::getUnaryRegionID<T>(f.getMesh(), rgn),
      f.getMesh()->template getRegion<ResT>(rgn),
      bout::detail::getPerpYIndex(f)};
}

template <typename ResT, typename L, typename R, typename Func>
inline auto SQ(const BinaryExpr<ResT, L, R, Func>& f) {
  using UnaryResT = bout::detail::UnaryFieldResult_t<ResT>;
  return BinaryExpr<UnaryResT, BinaryExpr<ResT, L, R, Func>, BinaryExpr<ResT, L, R, Func>,
                    bout::op::Square>{
      static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),
      static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),
      bout::op::Square{},
      f.getMesh(),
      f.getLocation(),
      f.getDirections(),
      f.getRegionID(),
      f.indices,
      bout::detail::getPerpYIndex(f)};
}

template <typename ResT, typename L, typename R, typename Func>
inline auto SQ(const BinaryExpr<ResT, L, R, Func>& f, const std::string& rgn) {
  using UnaryResT = bout::detail::UnaryFieldResult_t<ResT>;
  return BinaryExpr<UnaryResT, BinaryExpr<ResT, L, R, Func>, BinaryExpr<ResT, L, R, Func>,
                    bout::op::Square>{
      static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),
      static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(f),
      bout::op::Square{},
      f.getMesh(),
      f.getLocation(),
      f.getDirections(),
      bout::detail::getUnaryRegionID<UnaryResT>(f.getMesh(), rgn),
      f.getMesh()->template getRegion<UnaryResT>(rgn),
      bout::detail::getPerpYIndex(f)};
}

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
/// \f[\ln(\exp(f)) = f\f]
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
template <typename T, typename = bout::utils::EnableIfField<T>>
inline bool finite(const T& f, const std::string& rgn = "RGN_ALL") {

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

/// Makes a copy of a field \p f, ensuring that the underlying data is
/// not shared.
template <typename T, typename = bout::utils::EnableIfField<T>>
T copy(const T& f) {
  T result = f;
  result.allocate();
  return result;
}

class Field3DParallel;

/// Apply a floor value \p f to a field \p var. Any value lower than
/// the floor is set to the floor.
///
/// @param[in] var  Variable to apply floor to
/// @param[in] f    The floor value
/// @param[in] rgn  The region to calculate the result over
template <typename T, typename = bout::utils::EnableIfField<T>>
inline auto floor(const T& var, BoutReal f, const std::string& rgn = "RGN_ALL") {
  using ResT = bout::detail::UnaryFieldResult_t<T>;
  return BinaryExpr<ResT, T, Constant<BoutReal>, bout::op::Floor>{
      static_cast<typename T::View>(var),
      static_cast<typename Constant<BoutReal>::View>(f),
      bout::op::Floor{},
      var.getMesh(),
      var.getLocation(),
      var.getDirections(),
      bout::detail::getUnaryRegionID<T>(var.getMesh(), rgn),
      var.getMesh()->template getRegion<ResT>(rgn),
      bout::detail::getPerpYIndex(var)};
}

template <typename ResT, typename L, typename R, typename Func>
inline auto floor(const BinaryExpr<ResT, L, R, Func>& var, BoutReal f) {
  using UnaryResT = bout::detail::UnaryFieldResult_t<ResT>;
  return BinaryExpr<UnaryResT, BinaryExpr<ResT, L, R, Func>, Constant<BoutReal>,
                    bout::op::Floor>{
      static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(var),
      static_cast<typename Constant<BoutReal>::View>(f),
      bout::op::Floor{},
      var.getMesh(),
      var.getLocation(),
      var.getDirections(),
      var.getRegionID(),
      var.indices,
      bout::detail::getPerpYIndex(var)};
}

template <typename ResT, typename L, typename R, typename Func>
inline auto floor(const BinaryExpr<ResT, L, R, Func>& var, BoutReal f,
                  const std::string& rgn) {
  using UnaryResT = bout::detail::UnaryFieldResult_t<ResT>;
  return BinaryExpr<UnaryResT, BinaryExpr<ResT, L, R, Func>, Constant<BoutReal>,
                    bout::op::Floor>{
      static_cast<typename BinaryExpr<ResT, L, R, Func>::View>(var),
      static_cast<typename Constant<BoutReal>::View>(f),
      bout::op::Floor{},
      var.getMesh(),
      var.getLocation(),
      var.getDirections(),
      bout::detail::getUnaryRegionID<UnaryResT>(var.getMesh(), rgn),
      var.getMesh()->template getRegion<UnaryResT>(rgn),
      bout::detail::getPerpYIndex(var)};
}

#undef FIELD_FUNC

#endif /* FIELD_H */
