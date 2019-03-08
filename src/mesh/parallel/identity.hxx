/**************************************************************************
 * Parallel derivatives without any transform
 *
 **************************************************************************
 * Copyright 2015, 2016 B.D.Dudson, D. Dickinson
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

#ifndef __IDENTITYTRANSFORM_H__
#define __IDENTITYTRANSFORM_H__

#include "bout/paralleltransform.hxx"
#include "field3d.hxx"
#include "unused.hxx"

/*!
 * This class implements the simplest form of ParallelTransform
 * where the domain is a logically rectangular domain, and
 * yup() and ydown() refer to the same field.
 */
class ParallelTransformIdentity : public ParallelTransform {
public:
  ParallelTransformIdentity(Mesh& mesh_in) : ParallelTransform(mesh_in) {
    // check the coordinate system used for the grid data source
    checkInputGrid();
  }

  /*!
   * Merges the yup and ydown() fields of f, so that
   * f.yup() = f.ydown() = f
   */
  void calcParallelSlices(Field3D& f) override;

  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  const Field3D toFieldAligned(const Field3D& f, const REGION UNUSED(region)) override {
    return f;
  }

  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  const Field3D fromFieldAligned(const Field3D& f, const REGION UNUSED(region)) override {
    return f;
  }

  bool canToFromFieldAligned() override { return true; }

protected:
  void checkInputGrid() override;
};

#endif // __IDENTITYTRANSFORM_H__
