/**************************************************************************
 * Flux-coordinate Independent parallel derivatives
 *
 **************************************************************************
 * Copyright 2014 B.D.Dudson, P. Hill
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

#ifndef __FCITRANSFORM_H__
#define __FCITRANSFORM_H__

#include <bout/paralleltransform.hxx>
#include <bout/deprecated.hxx>
#include <interpolation.hxx>
#include <mask.hxx>
#include <parallel_boundary_region.hxx>
#include <unused.hxx>

/*!
 * Field line map - contains the coefficients for interpolation
 */
class FCIMap {
  /// Interpolation object
  Interpolation *interp;

  /// Private constructor - must be initialised with mesh
  FCIMap();
public:
  /// dir MUST be either +1 or -1
  FCIMap(Mesh& mesh, int dir, bool zperiodic);
  DEPRECATED(FCIMap(Mesh &mesh, int dir, bool UNUSED(yperiodic), bool zperiodic))
      : FCIMap(mesh, dir, zperiodic) {}

  int dir;                     /**< Direction of map */

  BoutMask boundary_mask;      /**< boundary mask - has the field line left the domain */
  Field3D y_prime;             /**< distance to intersection with boundary */

  BoundaryRegionPar* boundary; /**< boundary region */

  const Field3D interpolate(Field3D &f) const { return interp->interpolate(f); }
};

/*!
 * Flux Coordinate Independent method for parallel derivatives
 */
class FCITransform : public ParallelTransform {
public:
  DEPRECATED(FCITransform(Mesh &mesh, bool UNUSED(yperiodic), bool zperiodic))
      : FCITransform(mesh, zperiodic) {
        if (mesh.ystart > 1)
          // FCITransform can only use myg=1 because it only loads grid
          // information for one point forward or back along the magnetic
          // field, so it cannot set Field3D::yup2_field or
          // Field3D::ydown2_field
          throw BoutException("FCI method must use only one y-guard cell: set option myg=1");
      }
  FCITransform(Mesh &mesh, bool zperiodic = true)
      : mesh(mesh), forward_map(mesh, +1, zperiodic), backward_map(mesh, -1, zperiodic),
        zperiodic(zperiodic) {
          if (mesh.ystart > 1)
            // FCITransform can only use myg=1 because it only loads grid
            // information for one point forward or back along the magnetic
            // field, so it cannot set Field3D::yup2_field or
            // Field3D::ydown2_field
            throw BoutException("FCI method must use only one y-guard cell: set option myg=1");
        }

  void calcYUpDown(Field3D &f) override;

  const Field3D toFieldAligned(const Field3D &UNUSED(f), const REGION UNUSED(region)) override {
    throw BoutException("FCI method cannot transform into field aligned grid");
  }

  const Field3D fromFieldAligned(const Field3D &UNUSED(f), const REGION UNUSED(region)) override {
    throw BoutException("FCI method cannot transform into field aligned grid");
  }

  bool canToFromFieldAligned() override{
    return false;
  }
private:
  FCITransform();

  Mesh& mesh;

  FCIMap forward_map;           /**< FCI map for field lines in +ve y */
  FCIMap backward_map;          /**< FCI map for field lines in -ve y */

  bool zperiodic;               /**< Is the z-direction periodic? */
};

#endif // __FCITRANSFORM_H__
