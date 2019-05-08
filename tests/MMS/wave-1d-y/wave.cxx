#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>
#include <unused.hxx>

class Wave1D : public PhysicsModel {
private:
  Field3D f, g; // Evolving variables
  
protected:
  int init(bool UNUSED(restarting)) {

    g.setLocation(CELL_YLOW); // g staggered 
    
    // Tell BOUT++ to solve f and g
    bout_solve(f, "f");
    bout_solve(g, "g");

    return 0;
  }
  
  int rhs(BoutReal UNUSED(t)) {
    mesh->communicate(f,g); // Communicate guard cells
    
    // Central differencing
    ddt(f) = DDY(g, CELL_CENTRE);
    ddt(g) = DDY(f, CELL_YLOW);
    
    return 0;
  }
};

BOUTMAIN(Wave1D); // Create a main() function

