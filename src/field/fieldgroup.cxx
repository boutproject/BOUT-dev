#include <bout/fieldgroup.hxx>
#include <algorithm>

FieldGroup operator+(const FieldGroup& lhs, const FieldGroup& rhs) {
  return FieldGroup(lhs) += rhs;
}

void FieldGroup::makeUnique() {
  // Need to sort vector before making unique
  std::ranges::sort(fvec);

  // Remove duplicate entries (doesn't resize vector though)
  auto fvec_dupes = std::ranges::unique(fvec);

  // Resizes vector to remove memory no longer required
  fvec.erase(fvec_dupes.begin(), fvec_dupes.end());

  // Now do the same for the vector of Field3Ds
  std::ranges::sort(f3vec);
  auto f3vec_dupes = std::ranges::unique(f3vec);
  f3vec.erase(f3vec_dupes.begin(), f3vec_dupes.end());
}
