/*
 * 1D temperature conduction problem
 * 
 */

#include <bout/physicsmodel.hxx>

class Conduction : public PhysicsModel {
private:
  
  Field3D T; // Evolving temperature equation only
  
  BoutReal chi; // Parallel conduction coefficient

protected:

  // This is called once at the start
  int init(bool UNUSED(restarting)) override {
    
    // Get the options
    auto& options = Options::root()["conduction"];
    
    // Read from BOUT.inp, setting default to 1.0
    // The doc() provides some documentation in BOUT.settings
    chi = options["chi"].doc("Conduction coefficient").withDefault(1.0);

    // Tell BOUT++ to solve T
    SOLVE_FOR(T);
    
    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    T.getMesh()->communicate(T); // Communicate guard cells
    
    ddt(T) = Div_par_K_Grad_par(chi, T); // Parallel diffusion Div_{||}( chi * Grad_{||}(T) )
  
    return 0;
  }
};

BOUTMAIN(Conduction);
