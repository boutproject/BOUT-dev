/**************************************************************************
 * Shifted metric parallel derivatives, implemented with interpolation
 * instead of FFT
 *
 **************************************************************************
 * Copyright 2018 B.D.Dudson, P. Hill, J. Omotani, J. Parker
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

#ifndef __SHIFTEDINTERP_H__
#define __SHIFTEDINTERP_H__

#include <bout/paralleltransform.hxx>
#include <interpolation.hxx>

/*!
 * Shifted metric method
 * Each Y location is shifted in Z with respect to its neighbours
 * so that the grid is orthogonal in X-Z, but requires interpolation
 * to calculate the values of points along field-lines.
 *
 * In this implementation the interpolation is done using Interpolation objects
 */
class ShiftedMetricInterp : public ParallelTransform {
public:
  ShiftedMetricInterp() = delete;
  ShiftedMetricInterp(Mesh &mesh);
  ~ShiftedMetricInterp() {
    delete interp_yup;
    delete interp_ydown;
    delete interp_to_aligned;
    delete interp_from_aligned;
  }
  
  /*!
   * Calculates the yup() and ydown() fields of f
   * by interpolating f through a toroidal shift angle
   */ 
  void calcParallelSlices(Field3D &f) override;
  
  /*!
   * Uses interpolation of f through a toroidal shift angle to align the grid
   * points with the y coordinate (along magnetic field usually). 
   * 
   * Note that the returned field will no longer be orthogonal in X-Z, and the
   * metric tensor will need to be changed if X derivatives are used.
   */
  const Field3D toFieldAligned(const Field3D &f, REGION region = RGN_ALL) override;

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  const Field3D fromFieldAligned(const Field3D &f, REGION region = RGN_ALL) override;

  bool canToFromFieldAligned() override{
    return false;
  }

  std::vector<ParallelTransform::positionsAndWeights> getWeightsForYUpApproximation(int i, int j, int k) {
    return interp_yup->getWeightsForYApproximation(i,j,k,1);
  }
  std::vector<ParallelTransform::positionsAndWeights> getWeightsForYDownApproximation(int i, int j, int k) {
    return interp_ydown->getWeightsForYApproximation(i,j,k,-1);
  }

private:
  Mesh &localmesh; ///< The mesh this paralleltransform is part of

  /// This is the shift in toroidal angle (z) which takes a point from
  /// X-Z orthogonal to field-aligned along Y.
  Field2D zShift;

  /// Interpolation objects for yup and ydown transformations
  //Interpolation *interp_yup, *interp_ydown;

  /// Interpolation objects for yup and ydown transformations
  Interpolation *interp_yup, *interp_ydown;

  /// Interpolation objects for shifting to and from field-aligned coordinates
  Interpolation *interp_to_aligned, *interp_from_aligned;
};

#endif // __SHIFTEDINTERP_H__
