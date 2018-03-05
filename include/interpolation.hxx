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
#include "utils.hxx"

/// Interpolate to a give cell location
const Field3D interp_to(const Field3D &var, CELL_LOC loc, REGION region=RGN_ALL);
const Field2D interp_to(const Field2D &var, CELL_LOC loc, REGION region=RGN_ALL);

/// Print out the cell location (for debugging)
void printLocation(const Field3D &var);

const char* strLocation(CELL_LOC loc);


/// Interpolate a field onto a perturbed set of points
const Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z);

const Field3D interpolate(const Field2D &f, const Field3D &delta_x, const Field3D &delta_z);
const Field3D interpolate(const Field2D &f, const Field3D &delta_x);

////////////////////////////////////////

class Interpolation {
protected:
  // 3D vector of points to skip (true -> skip this point)
  BoutMask skip_mask;

  Mesh *localmesh;

public:
  Interpolation(int y_offset = 0, Mesh *mesh = nullptr)
      : localmesh(mesh), y_offset(y_offset) {
    if (mesh == nullptr) {
      localmesh = mesh;
    }
  }
  Interpolation(BoutMask mask, int y_offset = 0, Mesh *mesh = nullptr)
      : Interpolation(y_offset, mesh) {
    skip_mask = mask;
  }
  virtual ~Interpolation() {}

  virtual void calcWeights(const Field3D &delta_x, const Field3D &delta_z) = 0;
  virtual void calcWeights(const Field3D &delta_x, const Field3D &delta_z, BoutMask mask) = 0;

  virtual Field3D interpolate(const Field3D& f) const = 0;
  virtual Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z) = 0;
  virtual Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z, BoutMask mask) = 0;

  void setMask(BoutMask mask) { skip_mask = mask; }

  // Interpolate using the field at (x,y+y_offset,z), rather than (x,y,z)
  int y_offset;
  void setYOffset(int offset) { y_offset = offset; }
};

class HermiteSpline : public Interpolation {
  Tensor<int> i_corner;      // x-index of bottom-left grid point
  Tensor<int> k_corner;      // z-index of bottom-left grid point

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
  HermiteSpline(BoutMask mask, int y_offset=0, Mesh *mesh = nullptr) : HermiteSpline(y_offset, mesh) {
    skip_mask = mask;}

  /// Callback function for InterpolationFactory
  static Interpolation* CreateHermiteSpline(Mesh *mesh) {
    return new HermiteSpline(mesh);
  }

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z);
  void calcWeights(const Field3D &delta_x, const Field3D &delta_z, BoutMask mask);

  // Use precalculated weights
  Field3D interpolate(const Field3D& f) const;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z);
  Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z, BoutMask mask);

};

class Lagrange4pt : public Interpolation {
  Tensor<int> i_corner;      // x-index of bottom-left grid point
  Tensor<int> k_corner;      // z-index of bottom-left grid point

  Field3D t_x, t_z;

public:
  Lagrange4pt(Mesh *mesh = nullptr) : Lagrange4pt(0, mesh) {}
  Lagrange4pt(int y_offset = 0, Mesh *mesh = nullptr);
  Lagrange4pt(BoutMask mask, int y_offset=0, Mesh *mesh = nullptr) : Lagrange4pt(y_offset, mesh) {
    skip_mask = mask;}

  /// Callback function for InterpolationFactory
  static Interpolation* CreateLagrange4pt(Mesh *mesh) {
    return new Lagrange4pt(mesh);
  }

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z);
  void calcWeights(const Field3D &delta_x, const Field3D &delta_z, BoutMask mask);

  // Use precalculated weights
  Field3D interpolate(const Field3D& f) const;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z);
  Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z, BoutMask mask);
  BoutReal lagrange_4pt(BoutReal v2m,BoutReal vm,BoutReal vp,BoutReal v2p,BoutReal offset) const;
  BoutReal lagrange_4pt(const BoutReal v[],BoutReal offset) const;
};

class Bilinear : public Interpolation {
  Tensor<int> i_corner;      // x-index of bottom-left grid point
  Tensor<int> k_corner;      // z-index of bottom-left grid point

  Field3D w0, w1, w2, w3;

public:
  Bilinear(Mesh *mesh = nullptr) : Bilinear(0, mesh) {}
  Bilinear(int y_offset = 0, Mesh *mesh = nullptr);
  Bilinear(BoutMask mask, int y_offset=0, Mesh *mesh = nullptr) : Bilinear(y_offset, mesh) {
    skip_mask = mask;}

  /// Callback function for InterpolationFactory
  static Interpolation* CreateBilinear(Mesh *mesh) {
    return new Bilinear(mesh);
  }

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z);
  void calcWeights(const Field3D &delta_x, const Field3D &delta_z, BoutMask mask);

  // Use precalculated weights
  Field3D interpolate(const Field3D& f) const;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z);
  Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z, BoutMask mask);
};

#endif // __INTERP_H__
