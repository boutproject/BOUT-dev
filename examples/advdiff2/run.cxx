#include <bout.hxx>
#include <derivs.hxx>
#include "globals.hxx"

int physics_run(BoutReal UNUSED(t)) {
  // Run communications
  V.getMesh()->communicate(V);

  //ddt(V) = D2DX2(V) + 0.5*DDX(V) + D2DY2(V);
  ddt(V) = DDX(V);

  return 0;
}
