#include <bout.hxx>
#include <derivs.hxx>

#include "header.hxx"

int AdvDiff::rhs(BoutReal UNUSED(t)) {
  // Run communications
  mesh->communicate(V);

  ddt(V) = DDX(V);

  return 0;
}
