/**********************************************************************
 *
 * 1D modified advection-diffusion equation with operator splitting 
 * for use with the arkode solver
 *
 * Nick Walkden, 02/03/2015, nick.walkden@ccfe.ac.uk
 *
 * *******************************************************************/

#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

class IMEXexample : public PhysicsModel {
private:
  Field3D U;   // Evolving variable and auxilliary variable
  Field3D Vx, Dz; // Velocity, diffusion

protected:
  int init(bool) {
    setSplitOperator(); // Split into convective and diffusive
  
    // Get options
    auto& options = Options::root()["imex"];
    Vx = options["Vx"].doc("Velocity in X").withDefault(Field3D(100.0));
    Dz = options["Dz"].doc("Diffusion in Z").withDefault(Field3D(1.0));

    SOLVE_FOR(U);

    return 0;
  }

  int convective(BoutReal) {
    // Need communication
    U.getMesh()->communicate(U);
 
    // Passive advection
    ddt(U) = -VDDX(Vx,U);
   
    return 0;
  }
  
  int diffusive(BoutReal) {
    U.getMesh()->communicate(U);

    // Diffusion
    ddt(U) = Dz*D2DZ2(U);	
    
    return 0;
  }
};

BOUTMAIN(IMEXexample);
