/*
 * Abstract class, representing transformation to calculate
 * values along Y
 */

#ifndef __PARALLELTRANSFORM_H__
#define __PARALLELTRANSFORM_H__

#include <field3d.hxx>
#include <boutexception.hxx>
#include <dcomplex.hxx>

class Mesh;

#include <vector>

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
 */
class ParallelTransformIdentity : public ParallelTransform {
public:
  void calcYUpDown(Field3D &f) { }
  
  const Field3D toFieldAligned(const Field3D &f) {
    return f;
  }
  
  const Field3D fromFieldAligned(const Field3D &f) {
    return f;
  }
};

/*!
 * Shifted metric method
 */
class ShiftedMetric : public ParallelTransform {
public:
  ShiftedMetric(Mesh *m);
  
  void calcYUpDown(Field3D &f);
  
  const Field3D toFieldAligned(const Field3D &f);

  const Field3D fromFieldAligned(const Field3D &f);
private:
  ShiftedMetric();
  
  Mesh *mesh;

  Field2D zShift;
  std::vector<dcomplex> cmplx;
  
  const Field3D shiftZ(const Field3D f, const Field2D zangle);
  void shiftZ(const BoutReal *in, int len, BoutReal zangle,  BoutReal *out);
};


#endif // __PARALLELTRANSFORM_H__
