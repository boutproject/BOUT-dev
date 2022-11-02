/*!
 * This test implements 1D advection
 *
 */

#include <bout/physicsmodel.hxx>

class Advection : public PhysicsModel {
private:
  Field3D f; ///< The evolving field
  
protected:
  /// Initialise, specify evolving variables
  int init(bool) {
    // Solve for a single variable (f)
    SOLVE_FOR(f);
    return 0;
  }

  /// Calculate time derivatives
  /// 
  /// df/dt = 1 * df/dy
  ///
  /// Implemented using Vpar_Grad_par so 
  /// the method can be changed in the input
  int rhs(BoutReal) {
    mesh->communicate(f);
    
    ddt(f) = -Vpar_Grad_par(Field2D(1.0), f);
    
    return 0;
  }
};

BOUTMAIN(Advection);

