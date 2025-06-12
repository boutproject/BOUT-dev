/*!
 * \file facefield3d.hxx
 *
 * \brief Class for face-centered 3D fields. Stores field data on cell faces.
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

#pragma once
#ifndef BOUT_FACEFIELD3D_H
#define BOUT_FACEFIELD3D_H

#include "bout/field3d.hxx"
#include "bout/field_data.hxx"

/*!
 * \brief Represents a 3D field with components on cell faces
 * 
 * This class stores field data on the faces of computational cells
 * in all three directions (X, Y, Z). Each component is stored as
 * a Field3D with appropriate staggered location (CELL_XLOW, CELL_YLOW, CELL_ZLOW).
 * 
 * This is primarily used for flux calculations in finite-volume methods
 * where fluxes are naturally defined on cell faces.
 * 
 * Example usage:
 * \code{cpp}
 * FaceField3D flux(mesh);
 * flux = 1.0;  // Set all face values to 1.0
 * 
 * // Access individual components
 * Field3D& fx = flux.x();
 * Field3D& fy = flux.y();
 * Field3D& fz = flux.z();
 * \endcode
 */
template <typename T = BoutReal>
class FaceField3D_t : public FieldData {
public:
  /*!
   * \brief Constructor
   * 
   * Creates a FaceField3D with three Field3D components,
   * each allocated with the appropriate staggered location:
   * - x component: CELL_XLOW
   * - y component: CELL_YLOW
   * - z component: CELL_ZLOW
   * 
   * \param mesh The mesh to use for field allocation (optional)
   */
  explicit FaceField3D_t(Mesh* mesh = nullptr);

  /*!
   * \brief Copy constructor
   */
  FaceField3D_t(const FaceField3D_t& other);

  /*!
   * \brief Move constructor
   */
  FaceField3D_t(FaceField3D_t&& other) noexcept;

  /*!
   * \brief Destructor
   */
  ~FaceField3D_t() override = default;

  // Component access methods
  /*!
   * \brief Access the X-face component
   * \return Reference to the Field3D storing X-face data
   */
  Field3D& x() { return fx_; }
  const Field3D& x() const { return fx_; }

  /*!
   * \brief Access the Y-face component
   * \return Reference to the Field3D storing Y-face data
   */
  Field3D& y() { return fy_; }
  const Field3D& y() const { return fy_; }

  /*!
   * \brief Access the Z-face component
   * \return Reference to the Field3D storing Z-face data
   */
  Field3D& z() { return fz_; }
  const Field3D& z() const { return fz_; }

  // Assignment operators
  /*!
   * \brief Assign a scalar value to all face components
   */
  FaceField3D_t& operator=(T val);

  /*!
   * \brief Copy assignment
   */
  FaceField3D_t& operator=(const FaceField3D_t& rhs);

  /*!
   * \brief Move assignment
   */
  FaceField3D_t& operator=(FaceField3D_t&& rhs) noexcept;

  // Arithmetic operators
  /*!
   * \brief Add another FaceField3D component-wise
   */
  FaceField3D_t& operator+=(const FaceField3D_t& rhs);

  /*!
   * \brief Subtract another FaceField3D component-wise
   */
  FaceField3D_t& operator-=(const FaceField3D_t& rhs);

  /*!
   * \brief Multiply by a scalar
   */
  FaceField3D_t& operator*=(T rhs);

  /*!
   * \brief Multiply by a Field3D (component-wise)
   */
  FaceField3D_t& operator*=(const Field3D& rhs);

  /*!
   * \brief Divide by a scalar
   */
  FaceField3D_t& operator/=(T rhs);

  /*!
   * \brief Divide by a Field3D (component-wise)
   */
  FaceField3D_t& operator/=(const Field3D& rhs);

  // Unary operators
  /*!
   * \brief Unary negation
   */
  FaceField3D_t operator-() const;

  // Utility methods
  /*!
   * \brief Get the mesh associated with this field
   */
  Mesh* getMesh() const { return fx_.getMesh(); }

  /*!
   * \brief Apply boundary conditions to all components
   */
  void applyBoundary(bool init = false) override;

  /*!
   * \brief Apply boundary conditions to all components
   */
  void applyBoundary(const std::string& condition);

  /*!
   * \brief Apply boundary conditions to all components
   */
  void applyBoundary(const std::string& region, const std::string& condition);

  /*!
   * \brief Apply parallel boundary conditions to all components
   */
  void applyParallelBoundary(bool init = false) override;

  /*!
   * \brief Apply parallel boundary conditions to all components
   */
  void applyParallelBoundary(const std::string& condition);

  /*!
   * \brief Apply parallel boundary conditions to all components
   */
  void applyParallelBoundary(const std::string& region, const std::string& condition);

  // Virtual methods from FieldData
  bool isReal() const override { return true; }
  bool is3D() const override { return true; }
  int byteSize() const override { return 3 * sizeof(T); }
  int BoutRealSize() const override { return (is_real<T>::value) ? 3 : 0; }
  
#if BOUT_USE_TRACK
  std::string name;
#endif

  // MetaData interface
  void setLocation(CELL_LOC loc) override {
    // For FaceField3D, location is inherently face-based
    // This is a no-op as each component has its own location
  }
  
  CELL_LOC getLocation() const override {
    // Return a default - actual locations are per-component
    return CELL_LOC::deflt;
  }

  // Make FieldData::setBoundary visible
  using FieldData::setBoundary;
  
  void setBoundary(const std::string& name) override {
    fx_.setBoundary(name);
    fy_.setBoundary(name);
    fz_.setBoundary(name);
  }

  FieldData* timeDeriv() override;

  void setDirectionY(YDirectionType y) override {
    fx_.setDirectionY(y);
    fy_.setDirectionY(y);
    fz_.setDirectionY(y);
  }

private:
  Field3D fx_; ///< X-face component (located at CELL_XLOW)
  Field3D fy_; ///< Y-face component (located at CELL_YLOW)
  Field3D fz_; ///< Z-face component (located at CELL_ZLOW)

  // Track time derivative if allocated
  FaceField3D_t<T>* deriv{nullptr};
};

// Type alias for the common case
using FaceField3D = FaceField3D_t<BoutReal>;

// Binary operators
template <typename T>
FaceField3D_t<T> operator+(const FaceField3D_t<T>& lhs, const FaceField3D_t<T>& rhs);

template <typename T>
FaceField3D_t<T> operator-(const FaceField3D_t<T>& lhs, const FaceField3D_t<T>& rhs);

template <typename T>
FaceField3D_t<T> operator*(const FaceField3D_t<T>& lhs, T rhs);

template <typename T>
FaceField3D_t<T> operator*(T lhs, const FaceField3D_t<T>& rhs);

template <typename T>
FaceField3D_t<T> operator*(const FaceField3D_t<T>& lhs, const Field3D& rhs);

template <typename T>
FaceField3D_t<T> operator*(const Field3D& lhs, const FaceField3D_t<T>& rhs);

template <typename T>
FaceField3D_t<T> operator/(const FaceField3D_t<T>& lhs, T rhs);

template <typename T>
FaceField3D_t<T> operator/(const FaceField3D_t<T>& lhs, const Field3D& rhs);

#endif // BOUT_FACEFIELD3D_H