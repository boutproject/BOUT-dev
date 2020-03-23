/*!************************************************************
 * Main function for legacy interface
 *
 * The user provides two standalone functions, physics_init 
 * and physics_run. This header defines a class LegacyModel,
 * which inherits from PhysicsModel. This class is used by
 * Solver, and forwards calls to the user standalone functions.
 *
 * It is recommended that users define a model by interiting
 * from PhysicsModel directly, rather than using this header.
 *
 * Note: This header file implements a main() function, so
 * should only be included once, and cannot be combined with
 * a user-defined main().
 **************************************************************/

#if defined(__BOUTMAIN_H__) and not defined(BOUT_NO_USING_NAMESPACE_BOUTGLOBALS)
// Include using statement by default in user code.
// Macro allows us to include bout.hxx or physicsmodel.hxx without the using
// statement in library code.
using namespace bout::globals;
#endif

#ifndef __BOUTMAIN_H__
#define __BOUTMAIN_H__

#include <bout.hxx>
#include <options.hxx>
#include <msg_stack.hxx>
#include <boutexception.hxx>
#include <bout/physicsmodel.hxx>

#ifndef GLOBALORIGIN
#define GLOBAL extern
#else
#define GLOBAL
#endif

// Physics functions. The user implements
// these two standalone functions, which are
// called by the PhysicsModel class.

/// Initialise the model. Called once at the start
///
/// @param[in] restarting  True if the simulation is restarting
/// @returns zero on success, non-zero error code
int physics_init(bool restarting);

/// Calculate the time derivative
///
/// @param[in] t  Simulation time
/// @returns  zero on success, non-zero error code
int physics_run(BoutReal t);

/// Need a global Solver pointer,
/// which is the same as the PhysicsModel solver
Solver *solver;

/// Class interface to Solvers, which emulates
/// the older standalone function interface
class LegacyModel : public PhysicsModel {
protected:
  /// Initialise
  int init(bool restarting) override {
    // Set the global solver pointer
    ::solver = this->solver;
    
    // Call the standalone function
    return physics_init(restarting);
  }

  /// Calculate time derivatives
  int rhs(BoutReal t) override {
    return physics_run(t);
  }
  
};

/// Global functions used by some legacy models
template<class T>
void bout_solve(T &var, const char *name) {
  solver->add(var, name);
}

/*!
 * Add a contraint variable
 */
bool bout_constrain(Field3D &var, Field3D &F_var, const char *name) {
  if (!solver->constraints()) return false; // Doesn't support constraints

  // Add to solver
  solver->constraint(var, F_var, name);

  return true;
}

/// Main function
BOUTMAIN(LegacyModel);

#endif // __BOUTMAIN_H__
