#ifndef __IDENTITY_TRANSFORM_HXX__
#define __IDENTITY_TRANSFORM_HXX__

#include <bout/paralleltransform.hxx>

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
  void calcYUpDown(Field3D& f) override { f.mergeYupYdown(); }

  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  const Field3D toFieldAligned(const Field3D& f) override { return f; }

  /*!
   * The field is already aligned in Y, so this
   * does nothing
   */
  const Field3D fromFieldAligned(const Field3D& f) override { return f; }

  bool canToFromFieldAligned() override { return true; }
};

#endif // __IDENTITY_TRANSFORM_HXX__
