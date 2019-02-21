#ifndef __NEW_EMPTY_FIELD_H__
#define __NEW_EMPTY_FIELD_H__

#include "field.hxx"
#include "fieldperp.hxx"

/// Return an empty shell field of some type derived from Field, with metadata
/// copied but empty data array
template<typename T>
inline T emptyFrom(const T& f) {
  static_assert(std::is_base_of<Field, T>::value, "emptyFrom only works on Fields");
  return T(f.getMesh(), f.getLocation(), f.getDirectionX(), f.getDirectionY(), f.getDirectionZ()).allocate();
}

// Specialize newEmptyField templates for FieldPerp
template<>
inline FieldPerp emptyFrom<FieldPerp>(const FieldPerp& f) {
  return FieldPerp(f.getMesh(), f.getLocation(), f.getIndex(), f.getDirectionX(), f.getDirectionY(), f.getDirectionZ()).allocate();
}

#endif // __NEW_EMPTY_FIELD_H__
