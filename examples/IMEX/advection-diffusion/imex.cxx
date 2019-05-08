/**********************************************************************
 *
 * 1D modified advection-diffusion equation with operator splitting 
 * for use with the arkode solver
 *
 * Nick Walkden, 02/03/2015, nick.walkden@ccfe.ac.uk
 *
 * *******************************************************************/


#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>


int diffusive(BoutReal time);

Field3D U,Cx,Cz;   // Evolving variable and auxilliary variable

BoutReal cx,cz; //Advection velocity, diffusion rate

int physics_init(bool restarting) {

  // Give the solver two RHS functions
  // First function is explicit, second is implicit
  solver->setSplitOperator(physics_run,diffusive);
  
  // Get options
  auto globalOptions = Options::root();
  auto options = globalOptions["imex"];
  cz = options["cz"].withDefault(100.0);
  cx = options["cx"].withDefault(1.0);

  SOLVE_FOR(U);

  Cx = cx;
  Cz = cz;

  return 0;
}

int physics_run(BoutReal time) {

  // Need communication
  mesh->communicate(U);
 
  //Slow Passive advection
  ddt(U) = -VDDX(Cx,U);
   
  return 0;
}

int diffusive(BoutReal time) {

  mesh->communicate(U);
  //Fast passive advection
  
  ddt(U) = -VDDZ(Cz,U);	
  
  return 0;
}


