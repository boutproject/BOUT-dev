/*!
 * \file facefield3d.cxx
 *
 * \brief Implementation of FaceField3D_t class for face-centered fields
 *
 * \author BOUT++ Team
 *
 **************************************************************************
 * Copyright 2024 BOUT++ Contributors
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

#include "bout/facefield3d.hxx"
#include "bout/mesh.hxx"
#include "bout/globals.hxx"

template <typename T>
FaceField3D_t<T>::FaceField3D_t(Mesh* localmesh) 
    : fx_(localmesh, CELL_XLOW), 
      fy_(localmesh, CELL_YLOW), 
      fz_(localmesh, CELL_ZLOW) {
  if (localmesh == nullptr) {
    // Use global mesh if not specified
    localmesh = bout::globals::mesh;
  }
}

template <typename T>
FaceField3D_t<T>::FaceField3D_t(const FaceField3D_t& other)
    : fx_(other.fx_), fy_(other.fy_), fz_(other.fz_) {
#if BOUT_USE_TRACK
  name = other.name;
#endif
}

template <typename T>
FaceField3D_t<T>::FaceField3D_t(FaceField3D_t&& other) noexcept
    : fx_(std::move(other.fx_)), 
      fy_(std::move(other.fy_)), 
      fz_(std::move(other.fz_)) {
#if BOUT_USE_TRACK
  name = std::move(other.name);
#endif
  deriv = other.deriv;
  other.deriv = nullptr;
}

template <typename T>
FaceField3D_t<T>& FaceField3D_t<T>::operator=(T val) {
  fx_ = val;
  fy_ = val;
  fz_ = val;
  return *this;
}

template <typename T>
FaceField3D_t<T>& FaceField3D_t<T>::operator=(const FaceField3D_t& rhs) {
  if (this != &rhs) {
    fx_ = rhs.fx_;
    fy_ = rhs.fy_;
    fz_ = rhs.fz_;
#if BOUT_USE_TRACK
    name = rhs.name;
#endif
  }
  return *this;
}

template <typename T>
FaceField3D_t<T>& FaceField3D_t<T>::operator=(FaceField3D_t&& rhs) noexcept {
  if (this != &rhs) {
    fx_ = std::move(rhs.fx_);
    fy_ = std::move(rhs.fy_);
    fz_ = std::move(rhs.fz_);
#if BOUT_USE_TRACK
    name = std::move(rhs.name);
#endif
    deriv = rhs.deriv;
    rhs.deriv = nullptr;
  }
  return *this;
}

template <typename T>
FaceField3D_t<T>& FaceField3D_t<T>::operator+=(const FaceField3D_t& rhs) {
  fx_ += rhs.fx_;
  fy_ += rhs.fy_;
  fz_ += rhs.fz_;
  return *this;
}

template <typename T>
FaceField3D_t<T>& FaceField3D_t<T>::operator-=(const FaceField3D_t& rhs) {
  fx_ -= rhs.fx_;
  fy_ -= rhs.fy_;
  fz_ -= rhs.fz_;
  return *this;
}

template <typename T>
FaceField3D_t<T>& FaceField3D_t<T>::operator*=(T rhs) {
  fx_ *= rhs;
  fy_ *= rhs;
  fz_ *= rhs;
  return *this;
}

template <typename T>
FaceField3D_t<T>& FaceField3D_t<T>::operator*=(const Field3D& rhs) {
  fx_ *= rhs;
  fy_ *= rhs;
  fz_ *= rhs;
  return *this;
}

template <typename T>
FaceField3D_t<T>& FaceField3D_t<T>::operator/=(T rhs) {
  fx_ /= rhs;
  fy_ /= rhs;
  fz_ /= rhs;
  return *this;
}

template <typename T>
FaceField3D_t<T>& FaceField3D_t<T>::operator/=(const Field3D& rhs) {
  fx_ /= rhs;
  fy_ /= rhs;
  fz_ /= rhs;
  return *this;
}

template <typename T>
FaceField3D_t<T> FaceField3D_t<T>::operator-() const {
  FaceField3D_t<T> result(getMesh());
  result.fx_ = -fx_;
  result.fy_ = -fy_;
  result.fz_ = -fz_;
  return result;
}

template <typename T>
void FaceField3D_t<T>::applyBoundary(bool init) {
  fx_.applyBoundary(init);
  fy_.applyBoundary(init);
  fz_.applyBoundary(init);
}

template <typename T>
void FaceField3D_t<T>::applyBoundary(const std::string& condition) {
  fx_.applyBoundary(condition);
  fy_.applyBoundary(condition);
  fz_.applyBoundary(condition);
}

template <typename T>
void FaceField3D_t<T>::applyBoundary(const std::string& region, 
                                     const std::string& condition) {
  fx_.applyBoundary(region, condition);
  fy_.applyBoundary(region, condition);
  fz_.applyBoundary(region, condition);
}

template <typename T>
void FaceField3D_t<T>::applyParallelBoundary(bool init) {
  fx_.applyParallelBoundary(init);
  fy_.applyParallelBoundary(init);
  fz_.applyParallelBoundary(init);
}

template <typename T>
void FaceField3D_t<T>::applyParallelBoundary(const std::string& condition) {
  fx_.applyParallelBoundary(condition);
  fy_.applyParallelBoundary(condition);
  fz_.applyParallelBoundary(condition);
}

template <typename T>
void FaceField3D_t<T>::applyParallelBoundary(const std::string& region,
                                             const std::string& condition) {
  fx_.applyParallelBoundary(region, condition);
  fy_.applyParallelBoundary(region, condition);
  fz_.applyParallelBoundary(region, condition);
}

template <typename T>
FieldData* FaceField3D_t<T>::timeDeriv() {
  if (deriv == nullptr) {
    // Create a new FaceField3D for the time derivative
    deriv = new FaceField3D_t<T>(getMesh());
  }
  return deriv;
}

// Binary operators implementation
template <typename T>
FaceField3D_t<T> operator+(const FaceField3D_t<T>& lhs, const FaceField3D_t<T>& rhs) {
  FaceField3D_t<T> result(lhs);
  result += rhs;
  return result;
}

template <typename T>
FaceField3D_t<T> operator-(const FaceField3D_t<T>& lhs, const FaceField3D_t<T>& rhs) {
  FaceField3D_t<T> result(lhs);
  result -= rhs;
  return result;
}

template <typename T>
FaceField3D_t<T> operator*(const FaceField3D_t<T>& lhs, T rhs) {
  FaceField3D_t<T> result(lhs);
  result *= rhs;
  return result;
}

template <typename T>
FaceField3D_t<T> operator*(T lhs, const FaceField3D_t<T>& rhs) {
  return rhs * lhs;  // Commutative
}

template <typename T>
FaceField3D_t<T> operator*(const FaceField3D_t<T>& lhs, const Field3D& rhs) {
  FaceField3D_t<T> result(lhs);
  result *= rhs;
  return result;
}

template <typename T>
FaceField3D_t<T> operator*(const Field3D& lhs, const FaceField3D_t<T>& rhs) {
  return rhs * lhs;  // Commutative
}

template <typename T>
FaceField3D_t<T> operator/(const FaceField3D_t<T>& lhs, T rhs) {
  FaceField3D_t<T> result(lhs);
  result /= rhs;
  return result;
}

template <typename T>
FaceField3D_t<T> operator/(const FaceField3D_t<T>& lhs, const Field3D& rhs) {
  FaceField3D_t<T> result(lhs);
  result /= rhs;
  return result;
}

// Explicit instantiation for BoutReal
template class FaceField3D_t<BoutReal>;

// Explicit instantiation of binary operators
template FaceField3D operator+(const FaceField3D& lhs, const FaceField3D& rhs);
template FaceField3D operator-(const FaceField3D& lhs, const FaceField3D& rhs);
template FaceField3D operator*(const FaceField3D& lhs, BoutReal rhs);
template FaceField3D operator*(BoutReal lhs, const FaceField3D& rhs);
template FaceField3D operator*(const FaceField3D& lhs, const Field3D& rhs);
template FaceField3D operator*(const Field3D& lhs, const FaceField3D& rhs);
template FaceField3D operator/(const FaceField3D& lhs, BoutReal rhs);
template FaceField3D operator/(const FaceField3D& lhs, const Field3D& rhs);