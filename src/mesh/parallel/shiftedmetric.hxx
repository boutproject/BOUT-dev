/**************************************************************************
 * Shifted-metric parallel derivatives
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

#ifndef __SHIFTEDMETRIC_H__
#define __SHIFTEDMETRIC_H__

#include "bout/paralleltransform.hxx"

#include "field2d.hxx"

/*!
 * Shifted metric method
 * Each Y location is shifted in Z with respect to its neighbours
 * so that the grid is orthogonal in X-Z, but requires interpolation
 * to calculate the values of points along field-lines.
 *
 * In this implementation the interpolation is done using FFTs in Z
 */
class ShiftedMetric : public ParallelTransform {
public:
  ShiftedMetric() = delete;
  ShiftedMetric(Mesh& mesh, CELL_LOC location, Field2D zShift, BoutReal zlength_in);

  /*!
   * Calculates the yup() and ydown() fields of f
   * by taking FFTs in Z and applying a phase shift.
   */
  void calcParallelSlices(Field3D& f) override;

  /*!
   * Uses FFTs and a phase shift to align the grid points
   * with the y coordinate (along magnetic field usually).
   *
   * Note that the returned field will no longer be orthogonal
   * in X-Z, and the metric tensor will need to be changed
   * if X derivatives are used.
   */
  const Field3D toFieldAligned(const Field3D& f, const REGION region = RGN_ALL) override;

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  const Field3D fromFieldAligned(const Field3D& f,
                                 const REGION region = RGN_ALL) override;

  bool canToFromFieldAligned() override { return true; }

protected:
  void checkInputGrid() override;

private:
  CELL_LOC location{CELL_CENTRE};

  /// This is the shift in toroidal angle (z) which takes a point from
  /// X-Z orthogonal to field-aligned along Y.
  Field2D zShift;

  /// Length of the z-domain in radians
  BoutReal zlength{0.};

  int nmodes;

  Tensor<dcomplex> toAlignedPhs;   ///< Cache of phase shifts for transforming from X-Z
                                   /// orthogonal coordinates to field-aligned coordinates
  Tensor<dcomplex> fromAlignedPhs; ///< Cache of phase shifts for transforming from
                                   /// field-aligned coordinates to X-Z orthogonal
  /// coordinates

  /// Helper POD for parallel slice phase shifts
  struct ParallelSlicePhase {
    Tensor<dcomplex> phase_shift;
    int y_offset;
  };

  /// Cache of phase shifts for the parallel slices. Slices are stored
  /// in the following order:
  ///     {+1, ..., +n, -1, ..., -n}
  /// slice[i] stores offset i+1
  /// slice[2*i + 1] stores offset -(i+1)
  /// where i goes from 0 to (n-1), with n the number of y guard cells
  std::vector<ParallelSlicePhase> parallel_slice_phases;

  /*!
   * Shift a 2D field in Z.
   * Since 2D fields are constant in Z, this has no effect
   */
  const Field2D shiftZ(const Field2D& f, const Field2D& UNUSED(zangle),
                       const REGION UNUSED(region) = RGN_NOX) const {
    return f;
  };

  /*!
   * Shift a 3D field \p f in Z by the given \p zangle
   *
   * @param[in] f  The field to shift
   * @param[in] zangle   Toroidal angle (z)
   *
   */
  const Field3D shiftZ(const Field3D& f, const Field2D& zangle,
                       const REGION region = RGN_NOX) const;

  /*!
   * Shift a 3D field \p f by the given phase \p phs in Z
   *
   * Calculates FFT in Z, multiplies by the complex phase
   * and inverse FFTS.
   *
   * @param[in] f  The field to shift
   * @param[in] phs  The phase to shift by
   * @param[in] y_direction_out  The value to set yDirectionType of the result to
   */
  const Field3D shiftZ(const Field3D& f, const Tensor<dcomplex>& phs,
                       const YDirectionType y_direction_out,
                       const REGION region = RGN_NOX) const;

  /*!
   * Shift a given 1D array, assumed to be in Z, by the given \p zangle
   *
   * @param[in] in  A 1D array of length \p len
   * @param[in] len  Length of the in and out arrays
   * @param[in] zangle  The angle (z coordinate) to shift by
   * @param[out] out  A 1D array of length \p len, already allocated
   */
  void shiftZ(const BoutReal* in, int len, BoutReal zangle, BoutReal* out) const;

  /*!
   * Shift a given 1D array, assumed to be in Z, by the given \p zangle
   *
   * @param[in] in  A 1D array of length mesh.LocalNz
   * @param[in] phs Phase shift, assumed to have length (mesh.LocalNz/2 + 1) i.e. the
   * number of modes
   * @param[out] out  A 1D array of length mesh.LocalNz, already allocated
   */
  void shiftZ(const BoutReal* in, const dcomplex* phs, BoutReal* out) const;

  /// Calculate and store the phases for to/from field aligned and for
  /// the parallel slices using zShift
  void cachePhases();

  /// Shift a 3D field \p f in Z to all the parallel slices in \p phases
  ///
  /// @param[in] f      The field to shift
  /// @param[in] phases The phase and offset information for each parallel slice
  /// @return The shifted parallel slices
  std::vector<Field3D> shiftZ(const Field3D& f,
                              const std::vector<ParallelSlicePhase>& phases) const;
};

#endif // __SHIFTEDMETRIC_H__
