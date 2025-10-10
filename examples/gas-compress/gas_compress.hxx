#include <bout/physicsmodel.hxx>

class GasCompress : public PhysicsModel {
protected:
  int init(bool restarting) override;
  int rhs(BoutReal time) override;

private:
  // Evolving variables
  Field3D N, P; // Density, Pressure
  Vector3D V;   // velocity

  // 2D initial profiles
  Field2D N0, P0;
  Vector2D V0;

  // parameters
  BoutReal gamma_ratio; // Ratio of specific heats
  BoutReal nu;          // Viscosity
  bool include_viscosity;
  bool sub_initial; // Subtract initial force balance from momentum equation

  Vector2D g; // Acceleration
};
