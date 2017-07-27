
#include <bout/fieldgroup.hxx>

FieldGroup operator+(const FieldGroup &lhs, const FieldGroup &rhs) {
  return FieldGroup(lhs) += rhs;
}

void FieldGroup::makeUnique(){
  //Need to sort vector before making unique
  std::sort(fvec.begin(), fvec.end());
  //Remove duplicate entries (doesn't resize vector though)
  //auto last = std::unique(fvec.begin(), fvec.end());
  //Nicer to do the above but can't use auto until we have c++11 by default
  vector<FieldData*>::iterator last = std::unique(fvec.begin(), fvec.end());
  //Resizes vector to remove memory no longer required
  fvec.erase(last, fvec.end());
}
