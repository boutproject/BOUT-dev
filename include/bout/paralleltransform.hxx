/*
 * Abstract class, representing transformation to calculate
 * values along Y
 */

#ifndef __PARALLELTRANSFORM_H__
#define __PARALLELTRANSFORM_H__

#include <field3d.hxx>
#include <flexible.hxx>
#include <boutexception.hxx>
#include <dcomplex.hxx>
#include <unused.hxx>

class Mesh;

#include <vector>

/*!
 * Calculates the values of a field along the magnetic
 * field (y direction)
 *
 * This is used inside Mesh to represent a variety of possible
 * numerical schemes, including Shifted Metric and FCI
 *
 */
class ParallelTransform {
public:
  virtual ~ParallelTransform() {}

  /// Given a 3D field, calculate and set the Y up down fields
  virtual void calcYUpDown(Field3D &f) = 0;
  
  /// Convert a 3D field into field-aligned coordinates
  /// so that the y index is along the magnetic field
  virtual const Field3D toFieldAligned(const Field3D &f) = 0;
  
  /// Convert back from field-aligned coordinates
  /// into standard form
  virtual const Field3D fromFieldAligned(const Field3D &f) = 0;

  virtual bool canToFromFieldAligned() = 0;
};


/*!
 * This class implements the simplest form of ParallelTransform
 * where the domain is a logically rectangular domain, and 
 * yup() and ydown() refer to the same field.
 */
class ParallelTransformIdentity : public ParallelTransform {
public:
  /*!
   * Merges the yup and ydown() fields of f, so that
   * f.yup() = f.ydown() = f
   */ 
  void calcYUpDown(Field3D &f) override {f.mergeYupYdown();}
  
  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */ 
  const Field3D toFieldAligned(const Field3D &f) override {
    return f;
  }
  
  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  const Field3D fromFieldAligned(const Field3D &f) override {
    return f;
  }

  bool canToFromFieldAligned() override{
    return true;
  }
};

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
  ShiftedMetric(Mesh &mesh);
  
  /*!
   * Calculates the yup() and ydown() fields of f
   * by taking FFTs in Z and applying a phase shift.
   */ 
  void calcYUpDown(Field3D &f) override;
  
  /*!
   * Uses FFTs and a phase shift to align the grid points
   * with the y coordinate (along magnetic field usually). 
   * 
   * Note that the returned field will no longer be orthogonal
   * in X-Z, and the metric tensor will need to be changed 
   * if X derivatives are used.
   */
  const Field3D toFieldAligned(const Field3D &f) override;

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  const Field3D fromFieldAligned(const Field3D &f) override;

  bool canToFromFieldAligned() override{
    return true;
  }

  /// A 3D array, implemented as nested vectors
  typedef std::vector<std::vector<std::vector<dcomplex>>> arr3Dvec;
private:
  ShiftedMetric();

  Mesh &mesh; ///< The mesh this paralleltransform is part of

  /// This is the shift in toroidal angle (z) which takes a point from
  /// X-Z orthogonal to field-aligned along Y.
  std::vector<dcomplex> cmplx; ///< A temporary array, used for input/output to fft routines
  std::vector<dcomplex> cmplxLoc; ///< A temporary array, used for input/output to fft routines

  arr3Dvec getToAlignedPhs(CELL_LOC location = CELL_CENTRE); ///< Get phase shifts, calculating if necessary;
  arr3Dvec getFromAlignedPhs(CELL_LOC location = CELL_CENTRE); ///< Get phase shifts, calculating if necessary;
  arr3Dvec getYupPhs1(CELL_LOC location = CELL_CENTRE); ///< Get phase shifts, calculating if necessary;
  arr3Dvec getYdownPhs1(CELL_LOC location = CELL_CENTRE); ///< Get phase shifts, calculating if necessary;
  arr3Dvec getYupPhs2(CELL_LOC location = CELL_CENTRE); ///< Get phase shifts, calculating if necessary;
  arr3Dvec getYdownPhs2(CELL_LOC location = CELL_CENTRE); ///< Get phase shifts, calculating if necessary;

  arr3Dvec toAlignedPhs_CENTRE; ///< Cache of phase shifts for transforming from X-Z orthogonal coordinates to field-aligned coordinates. Cell centre version.
  arr3Dvec fromAlignedPhs_CENTRE; ///< Cache of phase shifts for transforming from field-aligned coordinates to X-Z orthogonal coordinates. Cell centre version.
  arr3Dvec toAlignedPhs_XLOW; ///< Cache of phase shifts for transforming from X-Z orthogonal coordinates to field-aligned coordinates. Interpolated to CELL_XLOW.
  arr3Dvec fromAlignedPhs_XLOW; ///< Cache of phase shifts for transforming from field-aligned coordinates to X-Z orthogonal coordinates. Interpolated to CELL_XLOW.
  arr3Dvec toAlignedPhs_YLOW; ///< Cache of phase shifts for transforming from X-Z orthogonal coordinates to field-aligned coordinates. Interpolated to CELL_YLOW.
  arr3Dvec fromAlignedPhs_YLOW; ///< Cache of phase shifts for transforming from field-aligned coordinates to X-Z orthogonal coordinates. Interpolated to CELL_YLOW.

  arr3Dvec yupPhs1_CENTRE; ///< Cache of phase shifts for calculating yup1 fields. Cell centre version.
  arr3Dvec ydownPhs1_CENTRE; ///< Cache of phase shifts for calculating ydown1 fields. Cell centre version.
  arr3Dvec yupPhs2_CENTRE; ///< Cache of phase shifts for calculating yup2 fields. Cell centre version.
  arr3Dvec ydownPhs2_CENTRE; ///< Cache of phase shifts for calculating ydown2 fields. Cell centre version.
  arr3Dvec yupPhs1_XLOW; ///< Cache of phase shifts for calculating yup1 fields. Interpolated to CELL_XLOW.
  arr3Dvec ydownPhs1_XLOW; ///< Cache of phase shifts for calculating ydown1 fields. Interpolated to CELL_XLOW.
  arr3Dvec yupPhs2_XLOW; ///< Cache of phase shifts for calculating yup2 fields. Interpolated to CELL_XLOW.
  arr3Dvec ydownPhs2_XLOW; ///< Cache of phase shifts for calculating ydown2 fields. Interpolated to CELL_XLOW.
  arr3Dvec yupPhs1_YLOW; ///< Cache of phase shifts for calculating yup1 fields. Interpolated to CELL_YLOW.
  arr3Dvec ydownPhs1_YLOW; ///< Cache of phase shifts for calculating ydown1 fields. Interpolated to CELL_YLOW.
  arr3Dvec yupPhs2_YLOW; ///< Cache of phase shifts for calculating yup2 fields. Interpolated to CELL_YLOW.
  arr3Dvec ydownPhs2_YLOW; ///< Cache of phase shifts for calculating ydown2 fields. Interpolated to CELL_YLOW.
  Flexible<Field2D> zShift;

  /*!
   * Shift a 2D field in Z. 
   * Since 2D fields are constant in Z, this has no effect
   */
  const Field2D shiftZ(const Field2D &f, const Field2D &UNUSED(zangle)){return f;};

  /*!
   * Shift a 3D field \p f in Z by the given \p zangle
   *
   * @param[in] f  The field to shift
   * @param[in] zangle   Toroidal angle (z)
   *
   */ 
  const Field3D shiftZ(const Field3D &f, const Field2D &zangle);

  /*!
   * Shift a 3D field \p f by the given phase \p phs in Z
   *
   * Calculates FFT in Z, multiplies by the complex phase
   * and inverse FFTS.
   *
   * @param[in] f  The field to shift
   * @param[in] phs  The phase to shift by
   */
  const Field3D shiftZ(const Field3D &f, const arr3Dvec &phs);

  /*!
   * Shift a given 1D array, assumed to be in Z, by the given \p zangle
   *
   * @param[in] in  A 1D array of length \p len
   * @param[in] len  Length of the in and out arrays
   * @param[in] zangle  The angle (z coordinate) to shift by
   * @param[out] out  A 1D array of length \p len, already allocated
   */
  void shiftZ(const BoutReal *in, int len, BoutReal zangle,  BoutReal *out);

  /*!
   * Shift a given 1D array, assumed to be in Z, by the given \p zangle
   *
   * @param[in] in  A 1D array of length mesh.LocalNz
   * @param[in] phs Phase shift, assumed to have length (mesh.LocalNz/2 + 1) i.e. the number of modes
   * @param[out] out  A 1D array of length mesh.LocalNz, already allocated
   */
  void shiftZ(const BoutReal *in, const std::vector<dcomplex> &phs, BoutReal *out);
};


#endif // __PARALLELTRANSFORM_H__
