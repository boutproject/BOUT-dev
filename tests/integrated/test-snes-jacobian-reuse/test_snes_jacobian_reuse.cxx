#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/unused.hxx>

/// Nonlinear relaxation to the steady state f = g = 1.
/// The Jacobian depends on the state, and the adaptive timestep
/// changes from step to step, so re-using the Jacobian exercises
/// both the timestep correction and the recalculation criteria.
class JacobianReuse : public PhysicsModel {
public:
  Field3D f, g;

  int init(bool UNUSED(restarting)) override {
    solver->add(f, "f");
    solver->add(g, "g");

    f = 0.1;
    g = 2.0;

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    ddt(f) = 1.0 - f * f * f;
    ddt(g) = f - g * g;
    return 0;
  }
};

BOUTMAIN(JacobianReuse);
