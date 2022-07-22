#include <bout/physicsmodel.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <cmath>
#include <bout/constants.hxx>
#include <unused.hxx>

BoutReal Lx, Ly, Lz; // Size of the domain

class Wave1D : public PhysicsModel {
private:
  Field3D f, g; // Evolving variables

  Coordinates *coord;
protected:
  int init(bool) override {
    // Coordinate system
    coord = mesh->getCoordinates();
    
    // Get the options
    Options &meshoptions = Options::root()["mesh"];
    
    Lx = meshoptions["Lx"].withDefault(1.0);
    Ly = meshoptions["Ly"].withDefault(1.0);
    
    // this assumes equidistant grid
    int nguard = mesh->xstart;
    coord->dx = Lx / (mesh->GlobalNx - 2 * nguard);
    coord->dy = Ly / (mesh->GlobalNy - 2 * nguard);

    SAVE_ONCE(Lx, Ly);
    
    //set mesh
    coord->g11 = 1.0;
    coord->g22 = 1.0;
    coord->g33 = 1.0;
    coord->g12 = 0.0;
    coord->g13 = 0.0;
    coord->g23 = 0.0;
    
    coord->g_11 = 1.0;
    coord->g_22 = 1.0;
    coord->g_33 = 1.0;
    coord->g_12 = 0.0;
    coord->g_13 = 0.0;
    coord->g_23 = 0.0;
    coord->geometry();

    g.setLocation(CELL_XLOW); // g staggered to the left of f
    
    // Tell BOUT++ to solve f and g
    SOLVE_FOR(f, g);

    return 0;
  }

  int rhs(BoutReal t) override {
    mesh->communicate(f,g); // Communicate guard cells
    
    //update time-dependent boundary conditions
    f.applyBoundary(t);
    g.applyBoundary(t);

    // Central differencing
    ddt(f) = DDX(g, CELL_CENTRE);// + 20*SQ(coord->dx)*D2DX2(f);
    ddt(g) = DDX(f, CELL_XLOW);// + 20*SQ(coord->dx)*D2DX2(g);
    
    return 0;
  }
};


///////////////////////////////////////////////////////

BOUTMAIN(Wave1D); // Create a main() function
