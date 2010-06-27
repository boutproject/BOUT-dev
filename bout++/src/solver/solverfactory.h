
class SolverFactory;

#ifndef __SOLVER_FACTORY_H__
#define __SOLVER_FACTORY_H__

#include "solver.h"

#include <string.h>

class SolverFactory {
 public:
  /// Return a pointer to the only instance
  static SolverFactory* getInstance();
  
  SolverType getDefaultSolverType();
  
  Solver* createSolver();
  Solver* createSolver(SolverType &);
  
 private:
  SolverFactory() {} // Prevent instantiation of this class
  static SolverFactory* instance; ///< The only instance of this class (Singleton)

};

#endif // __SOLVER_FACTORY_H__

