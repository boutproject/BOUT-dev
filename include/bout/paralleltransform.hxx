/*
 * Abstract class, representing transformation to calculate
 * values along Y
 */

#ifndef __PARALLELTRANSFORM_H__
#define __PARALLELTRANSFORM_H__

#include <field3d.hxx>
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
  /// Given a 3D field, calculate and set the Y up down fields
  virtual void calcYUpDown(Field3D &f) = 0;
  
  /// Convert a 3D field into field-aligned coordinates
  /// so that the y index is along the magnetic field
  virtual const Field3D toFieldAligned(const Field3D &f) = 0;
  
  /// Convert back from field-aligned coordinates
  /// into standard form
  virtual const Field3D fromFieldAligned(const Field3D &f) = 0;
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
  void calcYUpDown(Field3D &f) {f.mergeYupYdown();}
  
  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */ 
  const Field3D toFieldAligned(const Field3D &f) {
    return f;
  }
  
  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  const Field3D fromFieldAligned(const Field3D &f) {
    return f;
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
  void calcYUpDown(Field3D &f);
  
  /*!
   * Uses FFTs and a phase shift to align the grid points
   * with the y coordinate (along magnetic field usually). 
   * 
   * Note that the returned field will no longer be orthogonal
   * in X-Z, and the metric tensor will need to be changed 
   * if X derivatives are used.
   */
  const Field3D toFieldAligned(const Field3D &f);

  /*!
   * Converts a field back to X-Z orthogonal coordinates
   * from field aligned coordinates.
   */
  const Field3D fromFieldAligned(const Field3D &f);

  typedef std::vector<std::vector<std::vector<dcomplex>>> arr3Dvec;
private:
  ShiftedMetric();
  
  Mesh &mesh;

  Field2D zShift;
  std::vector<dcomplex> cmplx;
  std::vector<dcomplex> cmplxLoc;


  arr3Dvec toAlignedPhs;
  arr3Dvec fromAlignedPhs;

  arr3Dvec yupPhs;
  arr3Dvec ydownPhs;

  const Field2D shiftZ(const Field2D f, const Field2D UNUSED(zangle)){return f;};
  const Field3D shiftZ(const Field3D f, const Field2D zangle);
  const Field3D shiftZ(const Field3D f, const arr3Dvec &phs);
  void shiftZ(const BoutReal *in, int len, BoutReal zangle,  BoutReal *out);
  void shiftZ(const BoutReal *in, const std::vector<dcomplex> &phs, BoutReal *out);
};


#endif // __PARALLELTRANSFORM_H__
