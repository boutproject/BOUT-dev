// *******  Simulates 3D diffusion  *******

#include <bout.hxx>
#include <boutmain.hxx>
// This gives the Laplace(f) options
#include <difops.hxx>
// Gives PI and TWOPI
#include <bout/constants.hxx>

// Global variable initialization
// ############################################################################
Field3D n;               // The plasma density
BoutReal D_par, D_perp;  // The diffusion constants and the box dimensions
BoutReal Lx, Ly;         // The spatial domain size
bool Lx_Ly_from_grid;    // If the spatial size should be loaded from the grid
// ############################################################################

// Initialization of the physics
// ############################################################################
int physics_init(bool restarting) {

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
    flags->get("Lx_Ly_from_grid", Lx_Ly_from_grid, false);
    if(Lx_Ly_from_grid){
        // Loading variables from the grid file so that they can be saved into the
        // .dmp file
        GRID_LOAD2(Lx, Ly);
    }
    else{
        // Load from BOUT.inp
        Options *geometry = options->getSection("geom");
        geometry->get("Lx", Lx, TWOPI);
        geometry->get("Ly", Ly, TWOPI);
    }
    // ************************************************************************

    // Specify what values should be stored in the .dmp file
    // ************************************************************************
    // Save these variables once
    SAVE_ONCE2(Lx, Ly);

    // Tell BOUT++ to solve for n
    SOLVE_FOR(n);
    //*************************************************************************
    return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int physics_run(BoutReal t){
    mesh->communicate(n); // Communicate guard cells

    // Density diffusion
    ddt(n) = D_par *Laplace_par(n) + D_perp*Laplace_perp(n);
    return 0;
}
// ############################################################################
