#ifndef __SHIFTEDMETRIC_TRANSFORM_HXX__
#define __SHIFTEDMETRIC_TRANSFORM_HXX__

#include "bout/paralleltransform.hxx"
#include "dcomplex.hxx"
#include "field2d.hxx"
#include "unused.hxx"
#include <vector>

#include "bout_types.hxx"
#include "field3d.hxx"

class Mesh;

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
  ShiftedMetric(Mesh& mesh);

  /*!
   * Calculates the yup() and ydown() fields of f
   * by taking FFTs in Z and applying a phase shift.
   */
  void calcYUpDown(Field3D& f) override;

  /*!
   * Uses FFTs and a phase shift to align the grid points
   * with the y coordinate (along magnetic field usually).
   *
   * Note that the returned field will no longer be orthogonal
   * in X-Z, and the metric tensor will need to be changed
   * if X derivatives are used.
   */
  const Field3D toFieldAligned(const Field3D& f) override;

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  const Field3D fromFieldAligned(const Field3D& f) override;

  bool canToFromFieldAligned() override { return true; }

  /// A 3D array, implemented as nested vectors
  typedef std::vector<std::vector<std::vector<dcomplex>>> arr3Dvec;

private:
  Mesh& mesh; ///< The mesh this paralleltransform is part of

  /// This is the shift in toroidal angle (z) which takes a point from
  /// X-Z orthogonal to field-aligned along Y.
  Field2D zShift;
  std::vector<dcomplex>
      cmplx; ///< A temporary array, used for input/output to fft routines
  std::vector<dcomplex>
      cmplxLoc; ///< A temporary array, used for input/output to fft routines

  arr3Dvec toAlignedPhs;   ///< Cache of phase shifts for transforming from X-Z orthogonal
                           ///coordinates to field-aligned coordinates
  arr3Dvec fromAlignedPhs; ///< Cache of phase shifts for transforming from field-aligned
                           ///coordinates to X-Z orthogonal coordinates

  arr3Dvec yupPhs;   ///< Cache of phase shifts for calculating yup fields
  arr3Dvec ydownPhs; ///< Cache of phase shifts for calculating ydown fields

  /*!
   * Shift a 2D field in Z.
   * Since 2D fields are constant in Z, this has no effect
   */
  const Field2D shiftZ(const Field2D& f, const Field2D& UNUSED(zangle)) { return f; };

  /*!
   * Shift a 3D field \p f in Z by the given \p zangle
   *
   * @param[in] f  The field to shift
   * @param[in] zangle   Toroidal angle (z)
   *
   */
  const Field3D shiftZ(const Field3D& f, const Field2D& zangle);

  /*!
   * Shift a 3D field \p f by the given phase \p phs in Z
   *
   * Calculates FFT in Z, multiplies by the complex phase
   * and inverse FFTS.
   *
   * @param[in] f  The field to shift
   * @param[in] phs  The phase to shift by
   */
  const Field3D shiftZ(const Field3D& f, const arr3Dvec& phs);

  /*!
   * Shift a given 1D array, assumed to be in Z, by the given \p zangle
   *
   * @param[in] in  A 1D array of length \p len
   * @param[in] len  Length of the in and out arrays
   * @param[in] zangle  The angle (z coordinate) to shift by
   * @param[out] out  A 1D array of length \p len, already allocated
   */
  void shiftZ(const BoutReal* in, int len, BoutReal zangle, BoutReal* out);

  /*!
   * Shift a given 1D array, assumed to be in Z, by the given \p zangle
   *
   * @param[in] in  A 1D array of length mesh.LocalNz
   * @param[in] phs Phase shift, assumed to have length (mesh.LocalNz/2 + 1) i.e. the
   * number of modes
   * @param[out] out  A 1D array of length mesh.LocalNz, already allocated
   */
  void shiftZ(const BoutReal* in, const std::vector<dcomplex>& phs, BoutReal* out);
};

#endif // __SHIFTEDMETRIC_TRANSFORM_HXX__
