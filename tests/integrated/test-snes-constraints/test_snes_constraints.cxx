#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/unused.hxx>

/// Solves a Differential-Algebraic Equation (DAE) system
///
/// du/dt = -u + phi
/// phi - constraint_factor * u = 0
///
/// Eliminating phi gives the continuous reduced equation
///   du/dt = (-1 + constraint_factor) * u
///
/// with exact solution
///   u = u0 * exp( (-1 + constraint_factor) * t)
///
/// This integrated test advances a single implicit timestep with the SNES
/// solver, so the checked value is the one-step backward-Euler result rather
/// than the continuous exact solution.
class TestSnesConstraints : public PhysicsModel {
public:
  Field3D u;
  Field3D phi;

  constexpr BoutReal constraint_factor = 0.5;

  int init(bool UNUSED(restarting)) override {
    solver->add(u, "u");
    solver->constraint(phi, residual(phi), "phi");

    u = 1.0;
    phi = constraint_factor;

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    ddt(u) = -u + phi;
    residual(phi) = phi - constraint_factor * u;

    return 0;
  }
};

BOUTMAIN(TestSnesConstraints);
