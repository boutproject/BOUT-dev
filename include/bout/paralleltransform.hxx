/*
 * Abstract class, representing transformation to calculate
 * values along Y
 */

#ifndef __PARALLELTRANSFORM_H__
#define __PARALLELTRANSFORM_H__

#include "bout_types.hxx"

class Field3D;
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
  ParallelTransform(Mesh& mesh_in) : mesh(mesh_in) {}
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
  virtual const Field3D toFieldAligned(const Field3D &f, const REGION region = RGN_ALL) = 0;
  
  /// Convert back from field-aligned coordinates
  /// into standard form
  virtual const Field3D fromFieldAligned(const Field3D &f, const REGION region = RGN_ALL) = 0;

  virtual bool canToFromFieldAligned() = 0;

protected:
  /// This method should be called in the constructor to check that if the grid
  /// has a 'coordinates_type' variable, it has the correct value
  virtual void checkInputGrid() = 0;

  Mesh &mesh; ///< The mesh this paralleltransform is part of
};

#endif // __PARALLELTRANSFORM_H__
