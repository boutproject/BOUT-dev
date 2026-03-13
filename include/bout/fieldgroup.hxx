#ifndef BOUT_FIELDGROUP_H
#define BOUT_FIELDGROUP_H

#include <bout/traits.hxx>
#include <bout/vector2d.hxx>
#include <bout/vector3d.hxx>

#include <vector>

class Field2D;
class Field3D;
class FieldPerp;
class Field;

/// Group together fields for easier communication
///
/// Note: The `Field` class is used as a base class,
/// which is inherited by `Field2D`, `Field3D`, `FieldPerp`;
/// however `Vector2D` and `Vector3D` are stored by reference to their
/// components ``(x, y, z)`` as `Field2D` or `Field3D` objects.
class FieldGroup {
public:
  FieldGroup() = default;
  FieldGroup(const FieldGroup& other) = default;
  FieldGroup(FieldGroup&& other) = default;
  FieldGroup& operator=(const FieldGroup& other) = default;
  FieldGroup& operator=(FieldGroup&& other) = default;
  ~FieldGroup() = default;

  FieldGroup(Field& f) { fvec.push_back(&f); }
  FieldGroup(Field3D& f) {
    fvec.push_back(&f);
    f3vec.push_back(&f);
  }

  /// Constructor with a single Vector2D \p v
  ///
  /// This is needed so that fvec only contains Field2D or Field3D
  FieldGroup(Vector2D& v) {
    fvec.push_back(&v.x);
    fvec.push_back(&v.y);
    fvec.push_back(&v.z);
  }

  /// Constructor with a single Vector3D \p v
  ///
  /// This is needed so that fvec only contains Field2D or Field3D
  FieldGroup(Vector3D& v) {
    fvec.push_back(&v.x);
    fvec.push_back(&v.y);
    fvec.push_back(&v.z);
    f3vec.push_back(&v.x);
    f3vec.push_back(&v.y);
    f3vec.push_back(&v.z);
  }

  /// Variadic constructor. Allows an arbitrary number of
  /// Field arguments
  ///
  /// The explicit keyword prevents FieldGroup being constructed with arbitrary
  /// types. In particular arguments to add() cannot be implicitly converted
  /// to FieldGroup, leading to an infinite loop.
  template <typename... Ts>
  explicit FieldGroup(Ts&... ts) {
    add(ts...);
  }

  /// Copy contents of another FieldGroup \p other into this group.
  void add(const FieldGroup& other) {
    fvec.insert(fvec.end(), other.fvec.begin(), other.fvec.end());
    f3vec.insert(f3vec.end(), other.f3vec.begin(), other.f3vec.end());
  }

  /// Add the contents of \p other to this
  FieldGroup& operator+=(const FieldGroup& other) {
    add(other);
    return *this;
  }

  /// Add a Field \p f to the group.
  ///
  /// A pointer to this field will be stored internally,
  /// so the lifetime of this variable should be longer
  /// than the lifetime of this group.
  void add(Field& f) { fvec.push_back(&f); }

  // Add a 3D field \p f, which goes into both vectors.
  //
  // A pointer to this field will be stored internally,
  // so the lifetime of this variable should be longer
  // than the lifetime of this group.
  void add(Field3D& f) {
    fvec.push_back(&f);
    f3vec.push_back(&f);
  }

  /// Add a Vector2D \p v to the group.
  ///
  /// Pointers to this vector's components will be stored internally,
  /// so the lifetime of this variable should be longer than the
  /// lifetime of this group.
  void add(Vector2D& v) {
    fvec.push_back(&v.x);
    fvec.push_back(&v.y);
    fvec.push_back(&v.z);
  }

  /// Add a Vector3D \p v to the group.
  ///
  /// Pointers to this vector's components will be stored internally,
  /// so the lifetime of this variable should be longer than the
  /// lifetime of this group.
  void add(Vector3D& v) {
    fvec.push_back(&v.x);
    fvec.push_back(&v.y);
    fvec.push_back(&v.z);
    f3vec.push_back(&v.x);
    f3vec.push_back(&v.y);
    f3vec.push_back(&v.z);
  }

  /// Add multiple fields to this group
  template <typename... Ts>
  void add(Field& t, Ts&... ts) {
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

  /// Return number of fields
  int size() const { return static_cast<int>(fvec.size()); }

  /// Return number of Field3Ds
  int size_field3d() const { return static_cast<int>(f3vec.size()); }

  /// Test whether this group is empty
  bool empty() const { return fvec.empty(); }

  /// Remove all fields from this group
  void clear() {
    fvec.clear();
    f3vec.clear();
  }

  /// Iteration over all fields
  auto begin() { return fvec.begin(); }
  auto end() { return fvec.end(); }

  /// Const iteration over all fields
  auto begin() const { return fvec.cbegin(); }
  auto end() const { return fvec.cend(); }

  const std::vector<Field*>& get() const { return fvec; }

  /// Iteration over 3D fields
  const std::vector<Field3D*>& field3d() const { return f3vec; }

  /// Ensure that each field appears only once
  void makeUnique();

private:
  std::vector<Field*> fvec;    // Vector of fields
  std::vector<Field3D*> f3vec; // Vector of 3D fields
};

/// Combine two FieldGroups
FieldGroup operator+(const FieldGroup& lhs, const FieldGroup& rhs);

#endif // BOUT_FIELDGROUP_H
