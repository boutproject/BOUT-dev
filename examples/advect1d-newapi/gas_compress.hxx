#include <bout/physicsmodel.hxx>

class GasCompress : public PhysicsModel {
protected:
  int init(bool restarting);
  int rhs(BoutReal t);
private:
  // Evolving variables 
  Field3D N, P; // Density, Pressure
  Vector3D V;   // velocity
  
  // parameters
  BoutReal gamma_ratio;   // Ratio of specific heats
  BoutReal nu;      // Viscosity
  bool include_viscosity;
  
  Vector2D g; // Acceleration
};
