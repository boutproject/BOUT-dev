#include <bout/bout.hxx>
#include <bout/physicsmodel.hxx>

#include "header.hxx"

int AdvDiff::init(bool restarting) {
  // 2D initial profiles
  Field2D V0;

  // Read initial conditions
  mesh->get(V0, "V0");
  mesh->getCoordinates()->setDx(mesh->get("dx"));
  mesh->getCoordinates()->setDy(mesh->get("dy"));

  // read options

  // Set evolving variables
  SOLVE_FOR(V)

  if (!restarting) {
    // Set variables to these values (+ the initial perturbation)
    // NOTE: This must be after the calls to bout_solve
    V += V0;
  }

  return 0;
}

BOUTMAIN(AdvDiff)
