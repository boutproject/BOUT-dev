
#include <bout/physicsmodel.hxx>

#include <bout/fv_ops.hxx>

/// A 1D fluid equation (in y) 
/// 
/// Evolves density, pressure and momentum
///
class Fluid : public PhysicsModel {
protected:

  int init(bool UNUSED(restart)) override {
    auto& opt = Options::root()["fluid"];

    gamma = opt["gamma"]
      .doc("Adiabatic index (ratio of specific heats)")
      .withDefault(5. / 3);

    SOLVE_FOR(n, p, nv);

    return 0;
  }
  
  int rhs(BoutReal UNUSED(time)) override {

    mesh->communicate(n, p, nv);

    // Calculate velocity from momentum
    Field3D v = nv / n;
    
    // Calculate sound speed
    Field3D cs = sqrt(gamma * p / n);
    
    // Density equation
    ddt(n) =
      - FV::Div_par(n, v, cs) // Limiter, flux conserving
      ;
    
    // Pressure equation
    ddt(p) =
      - FV::Div_par(p, v, cs) // Limiter, flux conserving
      - (gamma - 1.0) * p * Div_par(v) // Standard central differencing
      ;
    
    // Momentum equation
    ddt(nv) =
      - FV::Div_par(nv, v, cs) // Limiter, flux conserving
      - Grad_par(p)     // Central differencing
      ;

    return 0;
  }
  
private:
  Field3D n, p, nv;  ///< Density, pressure, parallel momentum

  BoutReal gamma; ///< Adiabatic index
};

BOUTMAIN(Fluid);
