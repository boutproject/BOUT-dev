/*
 * Check stop tests
 * 
 */

#include <bout.hxx>
#include <bout/physicsmodel.hxx>
#include "unused.hxx"

class Test_stopcheck : public PhysicsModel {
protected:
  int init(bool UNUSED(restarting)) override;
  int rhs(BoutReal UNUSED(t)) override;
};


Field3D N;

int Test_stopcheck::init(bool UNUSED(restarting)) {
  solver->add(N,"N");
  return 0;
}

int Test_stopcheck::rhs(BoutReal UNUSED(t)) {
  bout::globals::mesh->communicate(N);
  ddt(N) = 0.;
  return 0;
}


BOUTMAIN(Test_stopcheck)
