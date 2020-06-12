// *******  Simulates 3D diffusion  *******

#include <bout.hxx>
#include <bout/physicsmodel.hxx>
// This gives the Laplace(f) options
#include <difops.hxx>
// Gives PI and TWOPI
#include <bout/constants.hxx>

class Diffusion_3d : public PhysicsModel {
protected:
  int init(bool UNUSED(restarting)) override;
  int rhs(BoutReal UNUSED(t)) override;
};


using bout::globals::mesh;

// Global variable initialization
// ############################################################################
Field3D n;               // Evolved variable
BoutReal D_par, D_perp;  // The diffusion constants
BoutReal Lx, Ly;         // The spatial domain size
bool use_grid;    // If the spatial size should be loaded from the grid
// ############################################################################

// Initialization of the physics
// ############################################################################
int Diffusion_3d::init(bool UNUSED(restarting)) {

    // Get the option (before any sections) in the BOUT.inp file
    Options *options   = Options::getRoot();

    // Get the diffusion constants
    // ************************************************************************
    // Get the section of the variables from [cst] specified in BOUT.inp
    // or in the command-line arguments
    Options *constants = options->getSection("cst");
    // Storing the variables with the following syntax
    // section_name->get("variable_name_in_input", variable_name_in_cxx,
    //                   default_value)
    constants->get("D_par",  D_par,  1.0);
    constants->get("D_perp", D_perp, 1.0);
    // ************************************************************************

    // Get domain dimensions
    // ************************************************************************
    Options *flags = options->getSection("flags");
    // Get the option
    flags->get("use_grid", use_grid, false);
    if(use_grid){
        // Loading variables from the grid file so that they can be saved into the
        // .dmp file (other variables such as dx, ny etc. are stored
        // automatically)
        GRID_LOAD2(Lx, Ly);
    }
    else{
        // Load from BOUT.inp
        Options *geometry = options->getSection("geom");
        geometry->get("Lx", Lx, 1.0);
        geometry->get("Ly", Ly, 1.0);
        // Calculate the internal number of points
        int internal_x_points = mesh->GlobalNx - 2*mesh->xstart;
        int internal_y_points = mesh->GlobalNy - 2*mesh->ystart;
        // Calculate dx and dy
        // dx = Lx/line_segments_in_x
        // On a line with equidistant points there is one less line
        // segment than points from the first to the last point.
        // The boundary lies (1/2)*dx away from the last point As there
        // are 2 boundaries there will effectively add one more line
        // segment in the domain. Hence
        mesh->getCoordinates()->dx = Lx/(internal_x_points);
        mesh->getCoordinates()->dy = Ly/(internal_y_points);
    }
    // ************************************************************************

    // Specify what values should be stored in the .dmp file
    // ************************************************************************
    // Save these variables once
    SAVE_ONCE4(Lx, Ly, D_par, D_perp);

    // Tell BOUT++ to solve for n
    SOLVE_FOR(n);
    //*************************************************************************
    return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int Diffusion_3d::rhs(BoutReal UNUSED(t)) {
    mesh->communicate(n); // Communicate guard cells

    // Density diffusion
    ddt(n) = D_par*Laplace_par(n) + D_perp*Laplace_perp(n);
    return 0;
}
// ############################################################################


BOUTMAIN(Diffusion_3d)
