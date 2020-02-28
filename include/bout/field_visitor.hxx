/*!
 * Visitor class for fields
 *
 * 
 */

class FieldVisitor;

#ifndef __FIELD_VISITOR_H__
#define __FIELD_VISITOR_H__

class Field2D;
class Field3D;
class FieldPerp;
class Vector2D;
class Vector3D;

class FieldVisitor {
public:
  virtual void accept(Field2D &f) = 0;
  virtual void accept(FieldPerp &f) = 0;
  virtual void accept(Field3D &f) = 0;
  virtual void accept(Vector2D &f) = 0;
  virtual void accept(Vector3D &f) = 0;
  virtual ~FieldVisitor() = default;
};

#endif // __FIELD_VISITOR_H__
