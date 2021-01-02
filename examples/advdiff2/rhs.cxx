#include <bout.hxx>
#include <derivs.hxx>

#include "header.hxx"

int AdvDiff::rhs(BoutReal UNUSED(t)) {
  // Run communications
  V.getMesh()->communicate(V);

  ddt(V) = DDX(V);

  return 0;
}
