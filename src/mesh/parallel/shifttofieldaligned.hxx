/**************************************************************************
 * Parallel derivatives use shifted-metric scheme implemented by shifting
 * to globally field-aligned coordinates
 *
 **************************************************************************
 * Copyright 2018 B.D.Dudson, P. Hill, J. Omotani
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

#ifndef __SHIFTTOFIELDALIGNED_H__
#define __SHIFTTOFIELDALIGNED_H__

#include <bout/paralleltransform.hxx>

class ShiftToFieldAligned : public ParallelTransform {
public:
  ShiftToFieldAligned(Mesh& mesh_in);
  /// Use an existing zShift
  ShiftToFieldAligned(Mesh& mesh_in, const Field2D& zShift) { init(mesh_in, zShift); }

  /*!
   * Do nothing here since we don't cache field-aligned fields
   */
  void calcYUpDown(Field3D& UNUSED(f)) override {}

  /*!
   * Uses FFTs and a phase shift to align the grid points
   * with the y coordinate (along magnetic field usually).
   *
   * Note that the returned field will no longer be orthogonal
   * in X-Z, and the metric tensor will need to be changed
   * if X derivatives are used.
   */
  const Field3D toFieldAligned(const Field3D& f, const REGION region = RGN_NOX) override {
    return implementations.at(f.getLocation()).toFieldAligned(f, region);
  }

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  const Field3D fromFieldAligned(const Field3D& f,
                                 const REGION region = RGN_NOX) override {
    return implementations.at(f.getLocation()).fromFieldAligned(f, region);
  }

  COORDINATE_SYSTEM getCoordinateSystem() const { return COORDINATE_SYSTEM::Orthogonal; }

  bool canToFromFieldAligned() override { return true; }

private:
  class Implementation {
  public:
    Implementation(Mesh& mesh_in, const CELL_LOC location_in, const Field2D& zShift_in);

    /*!
     * Uses FFTs and a phase shift to align the grid points
     * with the y coordinate (along magnetic field usually).
     *
     * Note that the returned field will no longer be orthogonal
     * in X-Z, and the metric tensor will need to be changed
     * if X derivatives are used.
     */
    const Field3D toFieldAligned(const Field3D& f, const REGION region = RGN_NOX);

    /*!
     * Converts a field back to X-Z orthogonal coordinates
     * from field aligned coordinates.
     */
    const Field3D fromFieldAligned(const Field3D& f, const REGION region = RGN_NOX);

  private:
    Mesh& mesh;
    const CELL_LOC location;
    Field2D zShift;

    /// Calculate and store the phases for to/from field aligned using zShift
    void cachePhases();

    int nmodes;

    Tensor<dcomplex> toAlignedPhs; ///< Cache of phase shifts for transforming from X-Z
                                   ///  orthogonal coordinates to field-aligned
                                   ///  coordinates.
    Tensor<dcomplex> fromAlignedPhs; ///< Cache of phase shifts for transforming from
                                     ///  field-aligned coordinates to X-Z orthogonal
                                     ///  coordinates.

    /*!
     * Shift a 3D field \p f by the given phase \p phs in Z
     *
     * Calculates FFT in Z, multiplies by the complex phase
     * and inverse FFTS.
     *
     * @param[in] f  The field to shift
     * @param[in] phs  The phase to shift by
     */
    const Field3D shiftZ(const Field3D& f, const Tensor<dcomplex>& phs,
                         const REGION region = RGN_NOX);

    /*!
     * Shift a given 1D array, assumed to be in Z, by the given \p zangle
     *
     * @param[in] in  A 1D array of length mesh.LocalNz
     * @param[in] phs Phase shift, assumed to have length (mesh.LocalNz/2 + 1)
     * i.e. the number of modes
     * @param[out] out  A 1D array of length mesh.LocalNz, already allocated
     */
    void shiftZ(const BoutReal* in, const dcomplex* phs, BoutReal* out);
  };

  std::map<CELL_LOC, Implementation> implementations;

  // populate the implementations map
  void init(Mesh& mesh_in, const Field2D& zShift);
};

#endif // __SHIFTTOFIELDALIGNED_H__
