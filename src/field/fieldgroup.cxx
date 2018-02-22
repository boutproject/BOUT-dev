
#include <bout/fieldgroup.hxx>

FieldGroup operator+(const FieldGroup &lhs, const FieldGroup &rhs) {
  return FieldGroup(lhs) += rhs;
}

void FieldGroup::makeUnique() {
  // Need to sort vector before making unique
  std::sort(fvec.begin(), fvec.end());

  // Remove duplicate entries (doesn't resize vector though)
  auto last = std::unique(fvec.begin(), fvec.end());

  // Resizes vector to remove memory no longer required
  fvec.erase(last, fvec.end());

  // Now do the same for the vector of Field3Ds
  std::sort(f3vec.begin(), f3vec.end());
  auto last_f3 = std::unique(f3vec.begin(), f3vec.end());
  f3vec.erase(last_f3, f3vec.end());
}
