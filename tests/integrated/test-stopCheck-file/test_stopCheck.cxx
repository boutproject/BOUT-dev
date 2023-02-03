/*
 * Check stop tests
 * 
 */

#include "bout/unused.hxx"
#include <bout/bout.hxx>
#include <bout/physicsmodel.hxx>

class Test_stopcheck : public PhysicsModel {
  Field3D N;

protected:
  int init(bool UNUSED(restarting)) override;
  int rhs(BoutReal UNUSED(t)) override;
};

int Test_stopcheck::init(bool UNUSED(restarting)) {
  solver->add(N, "N");
  return 0;
}

int Test_stopcheck::rhs(BoutReal UNUSED(t)) {
  bout::globals::mesh->communicate(N);
  ddt(N) = 0.;
  return 0;
}

BOUTMAIN(Test_stopcheck)
