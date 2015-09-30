
#include <bout/fieldgroup.hxx>

FieldGroup operator+(const FieldGroup &lhs, const FieldGroup &rhs) {
  return FieldGroup(lhs) += rhs;
};
