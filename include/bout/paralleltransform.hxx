/*
 * Abstract class, representing transformation to calculate
 * values along Y
 */

#ifndef __PARALLELTRANSFORM_H__
#define __PARALLELTRANSFORM_H__

#include <boutexception.hxx>
#include <dcomplex.hxx>
#include <field3d.hxx>
#include <unused.hxx>
#include <utils.hxx>

class Mesh;

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
  /// Constructor
  ParallelTransform(Mesh& mesh_in = *mesh) : thismesh(mesh_in) {}
  virtual ~ParallelTransform() {}

  /// Given a 3D field, calculate and set the Y up down fields
  virtual void calcYUpDown(Field3D &f) = 0;

  /// Calculate Yup and Ydown fields by integrating over mapped points
  /// This should be used for parallel divergence operators
  virtual void integrateYUpDown(Field3D &f) {
    return calcYUpDown(f);
  }
  
  /// Convert a 3D field into field-aligned coordinates
  /// so that the y index is along the magnetic field
  virtual const Field3D toFieldAligned(const Field3D &f, const REGION region = RGN_NOX) = 0;
  
  /// Convert back from field-aligned coordinates
  /// into standard form
  virtual const Field3D fromFieldAligned(const Field3D &f, const REGION region = RGN_NOX) = 0;

  /*!
   * Return the default coordinate system for Field3Ds when using this ParallelTransform
   */
  virtual CoordinateSystem getCoordinateSystem3D() const = 0;

  /*!
   * Return the default coordinate system for Field2Ds when using this ParallelTransform
   */
  virtual CoordinateSystem getCoordinateSystem2D() const = 0;

  virtual bool canToFromFieldAligned() = 0;

protected:
  Mesh& thismesh; ///< Mesh that this ParallelTransform belongs to
};


/*!
 * This class implements the simplest form of ParallelTransform
 * where the domain is a logically rectangular domain, and 
 * yup() and ydown() refer to the same field.
 */
class ParallelTransformIdentity : public ParallelTransform {
public:
  /// Constructor
  ParallelTransformIdentity(Mesh& mesh_in = *mesh)
    : ParallelTransform(mesh_in) {}
  /*!
   * Merges the yup and ydown() fields of f, so that
   * f.yup() = f.ydown() = f
   */ 
  void calcYUpDown(Field3D &f) override {f.mergeYupYdown();}
  
  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */ 
  const Field3D toFieldAligned(const Field3D &f, const REGION UNUSED(region)) override {
    return f;
  }
  
  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  const Field3D fromFieldAligned(const Field3D &f, const REGION UNUSED(region)) override {
    return f;
  }

  CoordinateSystem getCoordinateSystem3D() const override {
    return CoordinateSystem::FieldAligned;
  }

  CoordinateSystem getCoordinateSystem2D() const override {
    return CoordinateSystem::Axisymmetric;
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
  ShiftedMetric() = delete;
  /// Constructor
  ShiftedMetric(Mesh& mesh_in = *mesh);
  
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
  const Field3D toFieldAligned(const Field3D &f, const REGION region=RGN_NOX) override;

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  const Field3D fromFieldAligned(const Field3D &f, const REGION region=RGN_NOX) override;

  CoordinateSystem getCoordinateSystem3D() const override {
    return CoordinateSystem::Orthogonal;
  }

  CoordinateSystem getCoordinateSystem2D() const override {
    return CoordinateSystem::Axisymmetric;
  }

  bool canToFromFieldAligned() override{
    return true;
  }

private:
  /// This is the shift in toroidal angle (z) which takes a point from
  /// X-Z orthogonal to field-aligned along Y.
  Field2D zShift;

  int nmodes;

  Tensor<dcomplex> toAlignedPhs;   ///< Cache of phase shifts for transforming from X-Z
                                   ///orthogonal coordinates to field-aligned coordinates
  Tensor<dcomplex> fromAlignedPhs; ///< Cache of phase shifts for transforming from
                                   ///field-aligned coordinates to X-Z orthogonal
                                   ///coordinates

  Tensor<dcomplex> yupPhs;   ///< Cache of phase shifts for calculating yup fields
  Tensor<dcomplex> ydownPhs; ///< Cache of phase shifts for calculating ydown fields

  /*!
   * Shift a 2D field in Z. 
   * Since 2D fields are constant in Z, this has no effect
   */
  const Field2D shiftZ(const Field2D &f, const Field2D &UNUSED(zangle), const REGION UNUSED(region)=RGN_NOX){return f;};

  /*!
   * Shift a 3D field \p f in Z by the given \p zangle
   *
   * @param[in] f  The field to shift
   * @param[in] zangle   Toroidal angle (z)
   *
   */ 
  const Field3D shiftZ(const Field3D &f, const Field2D &zangle, const REGION region=RGN_NOX);

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
   * @param[in] in  A 1D array of length \p len
   * @param[in] len  Length of the in and out arrays
   * @param[in] zangle  The angle (z coordinate) to shift by
   * @param[out] out  A 1D array of length \p len, already allocated
   */
  void shiftZ(const BoutReal *in, int len, BoutReal zangle,  BoutReal *out);

  /*!
   * Shift a given 1D array, assumed to be in Z, by the given \p zangle
   *
   * @param[in] in  A 1D array of length thismesh.LocalNz
   * @param[in] phs Phase shift, assumed to have length (thismesh.LocalNz/2 + 1) i.e. the number of modes
   * @param[out] out  A 1D array of length thismesh.LocalNz, already allocated
   */
  void shiftZ(const BoutReal* in, const dcomplex* phs, BoutReal* out);
};


#endif // __PARALLELTRANSFORM_H__
