/**************************************************************************
 * Functions to interpolate between cell locations (e.g. lower Y and centred)
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
 **************************************************************************/

#ifndef __INTERP_H__
#define __INTERP_H__

#include "bout_types.hxx"
#include "field3d.hxx"
#include "mask.hxx"
#include "stencils.hxx"
#include "utils.hxx"

/// Perform interpolation between centre -> shifted or vice-versa
/*!
  Interpolate using 4th-order staggered formula

  @param[in] s  Input stencil. mm -> -3/2, m -> -1/2, p -> +1/2, pp -> +3/2
*/
inline BoutReal interp(const stencil& s) {
  return (9. * (s.m + s.p) - s.mm - s.pp) / 16.;
}

/// Interpolate to a give cell location
/*!
  Interpolate between different cell locations

  NOTE: This requires communication if the result is required in guard cells
  NOTE: Since corner guard cells cannot be communicated, it never makes sense
  to calculate interpolation in guard cells. If guard cell values are required,
  we must communicate (unless interpolating in z). Since mesh->communicate()
  communicates both x- and y-guard cells by default, there is no difference
  between RGN_ALL, RGN_NOX and RGN_NOY.

  @param[in]   var  Input variable
  @param[in]   loc  Location of output values
  @param[in]   region  Region where output will be calculated
*/
template <typename T>
const T interp_to(const T& var, CELL_LOC loc, REGION region = RGN_ALL) {
  AUTO_TRACE();
  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "interp_to must be templated with one of Field2D or Field3D.");

  Mesh* fieldmesh = var.getMesh();
  T result(fieldmesh);

  if ((loc != CELL_CENTRE && loc != CELL_DEFAULT) && (fieldmesh->StaggerGrids == false)) {
    throw BoutException("Asked to interpolate, but StaggerGrids is disabled!");
  }

  if (fieldmesh->StaggerGrids && (var.getLocation() != loc)) {

    // Staggered grids enabled, and need to perform interpolation
    TRACE("Interpolating %s -> %s", CELL_LOC_STRING(var.getLocation()).c_str(),
          CELL_LOC_STRING(loc).c_str());

    if (region != RGN_NOBNDRY) {
      // result is requested in some boundary region(s)
      result = var; // NOTE: This is just for boundaries. FIX!
    }
    // NOTE: invalidateGuards() is called in Field3D::alloctate() if the data
    // block is not already allocated, so will be called here if
    // region==RGN_NOBNDRY
    result.allocate();
    result.setLocation(loc);

    // Cell location of the input field
    const CELL_LOC location = var.getLocation();

    if ((location == CELL_CENTRE) || (loc == CELL_CENTRE)) {
      // Going between centred and shifted

      // Get the non-centre location for interpolation direction
      const CELL_LOC dir = (loc == CELL_CENTRE) ? location : loc;

      switch (dir) {
      case CELL_XLOW: {
        // At least 2 boundary cells needed for interpolation in x-direction
        ASSERT0(fieldmesh->xstart >= 2);

        if ((location == CELL_CENTRE) && (loc == CELL_XLOW)) { // C2L
          BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
            // Producing a stencil centred around a lower X value
            result[i] = interp(populateStencil<DIRECTION::X, STAGGER::C2L, 2>(var, i));
          }
        } else if (location == CELL_XLOW) { // L2C
          BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
            // Stencil centred around a cell centre
            result[i] = interp(populateStencil<DIRECTION::X, STAGGER::L2C, 2>(var, i));
          }
        }

        break;
      }
      case CELL_YLOW: {
        // At least 2 boundary cells needed for interpolation in y-direction
        ASSERT0(fieldmesh->ystart >= 2);

        if (var.hasYupYdown() && ((&var.yup() != &var) || (&var.ydown() != &var))) {
          // Field "var" has distinct yup and ydown fields which
          // will be used to calculate a derivative along
          // the magnetic field
          throw BoutException(
              "At the moment, fields with yup/ydown cannot use interp_to.\n"
              "If we implement a 3-point stencil for interpolate or double-up\n"
              "/double-down fields, then we can use this case.");

          if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) { // C2L
            BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
              // Producing a stencil centred around a lower X value
              result[i] = interp(
                  populateStencil<DIRECTION::YOrthogonal, STAGGER::C2L, 2>(var, i));
            }
          } else if (location == CELL_YLOW) { // L2C
            BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
              // Stencil centred around a cell centre
              result[i] = interp(
                  populateStencil<DIRECTION::YOrthogonal, STAGGER::L2C, 2>(var, i));
            }
          }
        } else {
          // var has no yup/ydown fields, so we need to shift into field-aligned
          // coordinates

          const T var_fa = fieldmesh->toFieldAligned(var);
          if (region != RGN_NOBNDRY) {
            // repeat the hack above for boundary points
            // this avoids a duplicate toFieldAligned call if we had called
            // result = toFieldAligned(result)
            // to get the boundary cells
            //
            // result is requested in some boundary region(s)
            result = var_fa; // NOTE: This is just for boundaries. FIX!
            result.allocate();
            result.setLocation(loc);
          }

          if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) { // C2L
            BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
              // Producing a stencil centred around a lower X value
              result[i] = interp(
                  populateStencil<DIRECTION::YAligned, STAGGER::C2L, 2>(var_fa, i));
            }
          } else if (location == CELL_YLOW) { // L2C
            BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
              // Stencil centred around a cell centre
              result[i] = interp(
                  populateStencil<DIRECTION::YAligned, STAGGER::L2C, 2>(var_fa, i));
            }
          }

          result = fieldmesh->fromFieldAligned(result);
        }
        break;
      }
      case CELL_ZLOW: {

        if ((location == CELL_CENTRE) && (loc == CELL_ZLOW)) { // C2L
          BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
            // Producing a stencil centred around a lower X value
            result[i] = interp(populateStencil<DIRECTION::Z, STAGGER::C2L, 2>(var, i));
          }
        } else if (location == CELL_ZLOW) { // L2C
          BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
            // Stencil centred around a cell centre
            result[i] = interp(populateStencil<DIRECTION::Z, STAGGER::L2C, 2>(var, i));
          }
        }
        break;
      }
      default: {
        // This should never happen
        throw BoutException("Unsupported direction of interpolation\n"
                            " - don't know how to interpolate to %s",
                            CELL_LOC_STRING(loc).c_str());
      }
      };

      if ((dir != CELL_ZLOW) && (region != RGN_NOBNDRY)) {
        fieldmesh->communicate(result);
      }

    } else {
      // Shifted -> shifted
      // For now, shift to centre then to final location loc
      // We probably should not rely on this, but it might work if one of the
      // shifts is in the z-direction where guard cells aren't needed.
      result = interp_to(interp_to(var, CELL_CENTRE), loc, region);
    }
    return result;
  } else {
    // Nothing to do - just return unchanged
    // Copying into result to return as returning var may increase the number of
    // references to the var data whilst returning result doesn't
    result = var;
    return result;
  }
}

/// Print out the cell location (for debugging)
void printLocation(const Field3D &var);

const char *strLocation(CELL_LOC loc);

/// Interpolate a field onto a perturbed set of points
const Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                          const Field3D &delta_z);

const Field3D interpolate(const Field2D &f, const Field3D &delta_x,
                          const Field3D &delta_z);
const Field3D interpolate(const Field2D &f, const Field3D &delta_x);

////////////////////////////////////////

class Interpolation {
protected:
  Mesh* localmesh{nullptr};

  // 3D vector of points to skip (true -> skip this point)
  BoutMask skip_mask;

public:
  Interpolation(int y_offset = 0, Mesh* localmeshIn = nullptr)
      : localmesh(localmeshIn == nullptr ? mesh : localmeshIn),
        skip_mask(*localmesh, false), y_offset(y_offset) {}
  Interpolation(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : Interpolation(y_offset, mesh) {
    skip_mask = mask;
  }
  virtual ~Interpolation() {}

  virtual void calcWeights(const Field3D &delta_x, const Field3D &delta_z) = 0;
  virtual void calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                           const BoutMask &mask) = 0;

  virtual Field3D interpolate(const Field3D &f) const = 0;
  virtual Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                              const Field3D &delta_z) = 0;
  virtual Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                              const Field3D &delta_z, const BoutMask &mask) = 0;

  void setMask(const BoutMask &mask) { skip_mask = mask; }

  // Interpolate using the field at (x,y+y_offset,z), rather than (x,y,z)
  int y_offset;
  void setYOffset(int offset) { y_offset = offset; }
};

class HermiteSpline : public Interpolation {
protected:
  /// This is protected rather than private so that it can be
  /// extended and used by HermiteSplineMonotonic
  
  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  // Basis functions for cubic Hermite spline interpolation
  //    see http://en.wikipedia.org/wiki/Cubic_Hermite_spline
  // The h00 and h01 basis functions are applied to the function itself
  // and the h10 and h11 basis functions are applied to its derivative
  // along the interpolation direction.

  Field3D h00_x;
  Field3D h01_x;
  Field3D h10_x;
  Field3D h11_x;
  Field3D h00_z;
  Field3D h01_z;
  Field3D h10_z;
  Field3D h11_z;

public:
  HermiteSpline(Mesh *mesh = nullptr) : HermiteSpline(0, mesh) {}
  HermiteSpline(int y_offset = 0, Mesh *mesh = nullptr);
  HermiteSpline(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : HermiteSpline(y_offset, mesh) {
    skip_mask = mask;
  }

  /// Callback function for InterpolationFactory
  static Interpolation *CreateHermiteSpline(Mesh *mesh) {
    return new HermiteSpline(mesh);
  }

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z) override;
  void calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                   const BoutMask &mask) override;

  // Use precalculated weights
  Field3D interpolate(const Field3D &f) const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                      const Field3D &delta_z) override;
  Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z,
                      const BoutMask &mask) override;
};


/// Monotonic Hermite spline interpolator
///
/// Similar to HermiteSpline, so uses most of the same code.
/// Forces the interpolated result to be in the range of the
/// neighbouring cell values. This prevents unphysical overshoots,
/// but also degrades accuracy near maxima and minima.
/// Perhaps should only impose near boundaries, since that is where
/// problems most obviously occur.
class MonotonicHermiteSpline : public HermiteSpline {
public:
  MonotonicHermiteSpline(Mesh *mesh = nullptr) : HermiteSpline(0, mesh) {}
  MonotonicHermiteSpline(int y_offset = 0, Mesh *mesh = nullptr)
      : HermiteSpline(y_offset, mesh) {}
  MonotonicHermiteSpline(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : HermiteSpline(mask, y_offset, mesh) {}

  /// Callback function for InterpolationFactory
  static Interpolation *CreateMonotonicHermiteSpline(Mesh *mesh) {
    return new MonotonicHermiteSpline(mesh);
  }
  
  /// Interpolate using precalculated weights.
  /// This function is called by the other interpolate functions
  /// in the base class HermiteSpline.
  Field3D interpolate(const Field3D &f) const override;
};

class Lagrange4pt : public Interpolation {
  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  Field3D t_x, t_z;

public:
  Lagrange4pt(Mesh *mesh = nullptr) : Lagrange4pt(0, mesh) {}
  Lagrange4pt(int y_offset = 0, Mesh *mesh = nullptr);
  Lagrange4pt(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : Lagrange4pt(y_offset, mesh) {
    skip_mask = mask;
  }

  /// Callback function for InterpolationFactory
  static Interpolation *CreateLagrange4pt(Mesh *mesh) { return new Lagrange4pt(mesh); }

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z) override;
  void calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                   const BoutMask &mask) override;

  // Use precalculated weights
  Field3D interpolate(const Field3D &f) const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                      const Field3D &delta_z) override;
  Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z,
                      const BoutMask &mask) override;
  BoutReal lagrange_4pt(BoutReal v2m, BoutReal vm, BoutReal vp, BoutReal v2p,
                        BoutReal offset) const;
  BoutReal lagrange_4pt(const BoutReal v[], BoutReal offset) const;
};

class Bilinear : public Interpolation {
  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  Field3D w0, w1, w2, w3;

public:
  Bilinear(Mesh *mesh = nullptr) : Bilinear(0, mesh) {}
  Bilinear(int y_offset = 0, Mesh *mesh = nullptr);
  Bilinear(const BoutMask &mask, int y_offset = 0, Mesh *mesh = nullptr)
      : Bilinear(y_offset, mesh) {
    skip_mask = mask;
  }

  /// Callback function for InterpolationFactory
  static Interpolation *CreateBilinear(Mesh *mesh) { return new Bilinear(mesh); }

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z) override;
  void calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                   const BoutMask &mask) override;

  // Use precalculated weights
  Field3D interpolate(const Field3D &f) const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                      const Field3D &delta_z) override;
  Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z,
                      const BoutMask &mask) override;
};

#endif // __INTERP_H__
