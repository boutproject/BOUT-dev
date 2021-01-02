#ifndef INCLUDE_GUARD_H
#define INCLUDE_GUARD_H

#include <bout/physicsmodel.hxx>

class AdvDiff : public PhysicsModel {
  // Evolving variables
  Field3D V;

protected:
  int init(bool restarting) override;
  int rhs(BoutReal) override;
};

#endif
