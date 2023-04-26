#include <bout/bout.hxx>
#include <bout/derivs.hxx>

#include "header.hxx"

int AdvDiff::rhs(BoutReal UNUSED(t)) {
  // Run communications
  V.getMesh()->communicate(V);

  ddt(V) = DDX(V);

  return 0;
}
