/*!************************************************************************
 * Main function
 *
 * BOUT++ can be used as a library, or by providing callback functions
 **************************************************************************/

#ifndef __BOUTMAIN_H__
#define __BOUTMAIN_H__

/*
 * This is a sloppy way to include main but for now it is a stop gap for converting all the examples
 */

#include <bout.hxx>
#include <options.hxx>
#include <msg_stack.hxx>
#include <boutexception.hxx>

#ifndef GLOBALORIGIN
#define GLOBAL extern
#else
#define GLOBAL
#endif

/*!
 * Solver object
 *
 * This is a global object here, to allow the physics models to be written as
 * standalone functions.
 */
Solver *solver;

/*!************************************************************************
 * Add variables to be solved
 **************************************************************************/

/*!
 * Add a variable to the solver. This will store a reference
 * to the given field, and associate that field 
 * with the given name.
 *
 * This function should probably be removed, replaced with either
 * solver->add(var, name)
 * 
 * or the macro
 * 
 * SOLVE_FOR(var)
 */
void bout_solve(Field2D &var, const char *name) {
  solver->add(var, name);
}

void bout_solve(Field3D &var, const char *name) {
  solver->add(var, name);
}

void bout_solve(Vector2D &var, const char *name) {
  solver->add(var, name);
}

void bout_solve(Vector3D &var, const char *name) {
  solver->add(var, name);
}

/// Macro to replace solver->add, passing variable name
#define SOLVE_FOR(var) solver->add(var, #var)
#define SOLVE_FOR2(var1, var2) { \
  solver->add(var1, #var1);       \
  solver->add(var2, #var2);}
#define SOLVE_FOR3(var1, var2, var3) { \
  solver->add(var1, #var1);             \
  solver->add(var2, #var2);             \
  solver->add(var3, #var3);}
#define SOLVE_FOR4(var1, var2, var3, var4) { \
  solver->add(var1, #var1);             \
  solver->add(var2, #var2);             \
  solver->add(var3, #var3);             \
  solver->add(var4, #var4);}
#define SOLVE_FOR5(var1, var2, var3, var4, var5) { \
  solver->add(var1, #var1);             \
  solver->add(var2, #var2);             \
  solver->add(var3, #var3);             \
  solver->add(var4, #var4);             \
  solver->add(var5, #var5);}
#define SOLVE_FOR6(var1, var2, var3, var4, var5, var6) { \
  solver->add(var1, #var1);             \
  solver->add(var2, #var2);             \
  solver->add(var3, #var3);             \
  solver->add(var4, #var4);             \
  solver->add(var5, #var5);             \
  solver->add(var6, #var6);}

/*!************************************************************************
 * Add constraints
 **************************************************************************/

/*!
 * Add a contraint variable
 */
bool bout_constrain(Field3D &var, Field3D &F_var, const char *name) {
  if (!solver->constraints()) return false; // Doesn't support constraints

  // Add to solver
  solver->constraint(var, F_var, name);

  return true;
}

// Physics functions
int physics_init(bool restarting);
int physics_run(BoutReal t);

/*!************************************************************************
 * Main function
 **************************************************************************/

int main(int argc, char **argv) {
  /// Start BOUT++
  int init_err = BoutInitialise(argc, argv);
  if (init_err < 0) {
    // User printed help message
    return 0;
  } else if (init_err > 0) {
    // Other errors
    return init_err;
  }

  try {
    /// Create the solver
    solver = Solver::create();

    /// Get the restart option
    bool restart;
    OPTION(Options::getRoot(), restart, false);

    /// Call the physics initialisation code
    output.write("Initialising physics module\n");
    int msg_point = msg_stack.push("Initialising physics module");
    if (physics_init(restart)) {
      output.write("Failed to initialise physics. Aborting\n");
      delete solver;
      MPI_Abort(BoutComm::get(), 1);
      BoutFinalise();
      return 1;
    }
    msg_stack.pop(msg_point);
    
    solver->outputVars(dump); // Add evolving variables to the output file
    
    bout_run(solver, physics_run);
    
    delete solver; // Delete the solver
  }catch (BoutException &e) {
    output << "Error encountered\n";
    output << e.what() << endl;
    // Shut down all processes by calling Abort
    MPI_Abort(BoutComm::get(), 1);
  }

  BoutFinalise();

  return(0);
}

#endif // __BOUTMAIN_H__
