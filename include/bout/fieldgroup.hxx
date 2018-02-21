#ifndef __FIELDGROUP_H__
#define __FIELDGROUP_H__

#include "bout/field_data.hxx"
#include <bout/field3d.hxx>

#include <bout/vector2d.hxx>
#include <bout/vector3d.hxx>

#include <vector>

#include <algorithm>

/// Group together fields
///
/// Note: The FieldData class is used as a base class,
/// which is inherited by Field2D, Field3D, Vector2D and Vector3D
/// however Vector2D and Vector3D are stored by reference to their
/// components (x,y,z) as Field2D or Field3D objects.
class FieldGroup {
 public:
  FieldGroup() {}

  FieldGroup(const FieldGroup &other) : fvec(other.fvec), f3vec(other.f3vec) {}

  FieldGroup(FieldData &f) {fvec.push_back(&f); }

  /// Constructor with a single Field3D
  FieldGroup(Field3D &f) {fvec.push_back(&f); f3vec.push_back(&f); }

  /// Constructor with a single Vector2D
  /// This is needed so that fvec only contains Field2D or Field3D
  FieldGroup(Vector2D &v) {
    fvec.push_back(&v.x);
    fvec.push_back(&v.y);
    fvec.push_back(&v.z);
  }

  FieldGroup(Vector3D &v) {
    fvec.push_back(&v.x);
    fvec.push_back(&v.y);
    fvec.push_back(&v.z);
    f3vec.push_back(&v.x);
    f3vec.push_back(&v.y);
    f3vec.push_back(&v.z);
  }

  /*
   * Variadic constructor. Allows an arbitrary number of
   * FieldData arguments
   *
   * The explicit keyword prevents FieldGroup being constructed with arbitrary
   * types. In particular arguments to add() cannot be implicitly converted
   * to FieldGroup, leading to an infinite loop.
   */
  template <typename... Ts>
  explicit FieldGroup(Ts&... ts) { add(ts...); }


  /*!
   * Copy contents of another FieldGroup into this group.
   */
  void add(const FieldGroup &other) {
    fvec.insert(fvec.end(), other.fvec.begin(), other.fvec.end() );
    f3vec.insert(f3vec.end(), other.f3vec.begin(), other.f3vec.end() );
  }

  /*!
   * Add another FieldGroup's contents
   */
  FieldGroup& operator+=(const FieldGroup &other) {
    add(other); return *this;
  }

  /*!
   * Add a FieldData to the group.
   *
   * A pointer to this field will be stored internally,
   * so the lifetime of this variable should be longer
   * than the lifetime of this group.
   */
  void add(FieldData &f) {
    fvec.push_back(&f);
  }

  /*!
   * Add a 3D field, which goes into both vectors.
   *
   * A pointer to this field will be stored internally,
   * so the lifetime of this variable should be longer
   * than the lifetime of this group.
   */
  void add(Field3D &f) {
    fvec.push_back(&f);
    f3vec.push_back(&f);
  }

  void add(Vector2D &v) {
    fvec.push_back(&v.x);
    fvec.push_back(&v.y);
    fvec.push_back(&v.z);
  }

  void add(Vector3D &v) {
    fvec.push_back(&v.x);
    fvec.push_back(&v.y);
    fvec.push_back(&v.z);
    f3vec.push_back(&v.x);
    f3vec.push_back(&v.y);
    f3vec.push_back(&v.z);
  }

  /*!
   * add( FieldData ... )
   *
   * Add fields to this group. This is a variadic template
   * which allows Field3D objects to be treated as a special
   * case. An arbitrary number of fields can be added.
   */
  template <typename... Ts>
  void add(FieldData& t, Ts&... ts) {
    add(t);     // Add the first using functions above
    add(ts...); // Add the rest
  }

  template <typename... Ts>
  void add(Field3D& t, Ts&... ts) {
    add(t);     // Add the first using functions above
    add(ts...); // Add the rest
  }

  template <typename... Ts>
  void add(Vector3D& t, Ts&... ts) {
    add(t);     // Add the first using functions above
    add(ts...); // Add the rest
  }

  template <typename... Ts>
  void add(Vector2D& t, Ts&... ts) {
    add(t);     // Add the first using functions above
    add(ts...); // Add the rest
  }

  /*!
   * Return number of fields
   */
  int size() const {
    return fvec.size();
  }

  /*!
   * Return number of Field3Ds
   */
  int size_field3d() const {
    return f3vec.size();
  }

  /*!
   * Test whether this group is empty
   */
  bool empty() const {
    return fvec.empty();
  }

  /*!
   * Remove all fields from this group
   */
  void clear() {fvec.clear(); f3vec.clear(); }

  /*
   * Iteration over all fields
   */
  typedef std::vector<FieldData*>::iterator iterator;
  iterator begin() {
    return fvec.begin();
  }
  iterator end() {
    return fvec.end();
  }

  /*
   * Const iteration over all fields
   */
  typedef std::vector<FieldData*>::const_iterator const_iterator;
  const_iterator begin() const {
    return fvec.begin();
  }
  const_iterator end() const {
    return fvec.end();
  }

  const std::vector<FieldData*>& get() const {
    return fvec;
  }

  /*
   * Iteration over 3D fields
   */
  const std::vector<Field3D*>& field3d() const {
    return f3vec;
  }

  /*
   * Ensure that each field appears only once
   */
  void makeUnique();
 private:
  std::vector<FieldData*> fvec;  // Vector of fields
  std::vector<Field3D*>   f3vec; // Vector of 3D fields
};


FieldGroup operator+(const FieldGroup &lhs, const FieldGroup &rhs);

#endif // __FIELDGROUP_H__
