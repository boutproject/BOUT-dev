/*
 * Non-linear diffusion
 *
 * This test case solves the diffusion equation along y:
 *
 * df/dt = Div_par( D Grad_par( f ) )
 *
 * using a nonlinear coefficient D, which depends on f
 * to a specified power alpha:
 *
 * D = f ^ alpha
 *
 * This is a model for parallel heat conduction in plasmas,
 * where the heat conduction coefficient depends on T^{5/2}
 *
 */

#include <bout/physicsmodel.hxx>
#include <bout/fv_ops.hxx>

class DiffusionNL : public PhysicsModel {
protected:
  int init(bool restarting) {
    // Get the input parameter alpha
    Options *opt = Options::getRoot();
    OPTION(opt, alpha, 2.5);
    
    // Specify that the operator is split
    // into convective and diffusive parts
    setSplitOperator();

    SOLVE_FOR(f);

    return 0;
  }
  int convective(BoutReal time) {
    ddt(f) = 0.0;
    return 0;
  }
  int diffusive(BoutReal time, bool linear) {
    if(!linear) {
      // Update diffusion coefficient
      D = f ^ alpha;
    }

    ddt(f) = FV::Div_par_K_Grad_par(D, f);

    return 0;
  }
private:
  Field3D f;      // Evolving quantity
  BoutReal alpha; // Input parameter, D = f^alpha
  Field3D D;      // The diffusion coefficient
};

BOUTMAIN(DiffusionNL);
