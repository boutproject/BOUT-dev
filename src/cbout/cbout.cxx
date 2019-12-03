#include <cbout/bout.h>

#include <bout.hxx>
#include <field3d.hxx>
#include <options.hxx>

#include <vector>
#include <string>

//////////////////////////////////////////////////////////////
// Error handling

#define CWRAP_START() \
  try {

#define CWRAP_END() \
  } catch (const BoutException& e) { \
    output.write(e.what());          \
    return 1;                        \
  }                                  \
  return BOUT_SUCCESS

//////////////////////////////////////////////////////////////
// Initialisation

extern "C" int bout_initialise(int argc, char** argv) {
  return BoutInitialise(argc, argv);
}

extern "C" int bout_finalise() {
  return BoutFinalise(true);
}

//////////////////////////////////////////////////////////////
// Fields

extern "C" int Field3D_create(Field3D** field) {
  CWRAP_START();
  *field = new Field3D();
  CWRAP_END();
}

extern "C" int Field3D_create_zerofrom(Field3D** field, Field3D* like) {
  CWRAP_START();
  *field = new Field3D(zeroFrom(*like));
  CWRAP_END();
}

extern "C" int Field3D_create_emptyfrom(Field3D** field, Field3D* like) {
  CWRAP_START();
  *field = new Field3D(emptyFrom(*like));
  CWRAP_END();
}

// Destroying

extern "C" int Field3D_delete(Field3D* field) {
  CWRAP_START();
  delete field;
  CWRAP_END();
}

// manipulating

extern "C" int Field3D_axpy(Field3D* result, BoutReal a, Field3D* x, Field3D* y) {
  CWRAP_START();
  (*result) = a * (*x) + (*y);
  CWRAP_END();
}

extern "C" int Field3D_scale(Field3D* x, BoutReal alpha) {
  CWRAP_START();
  (*x) = (*x) * alpha;
  CWRAP_END();
}

extern "C" int Field3D_getarray(Field3D* x, BoutReal **data) {
  CWRAP_START();
  *data = (*x)(0,0); // Pointer to start of internal data
  CWRAP_END();
}

//////////////////////////////////////////////////////////////
// Options
extern "C" Options* Options_root() {
  return Options::getRoot();
}

extern "C" Options* Options_get(Options* option, const char* name) {
  return &(*option)[name];
}

extern "C" int Options_get_int(int* result, Options* option) {
  CWRAP_START();
  *result = option->as<int>();
  CWRAP_END();
}

extern "C" int Options_get_BoutReal(BoutReal* result, Options* option) {
  CWRAP_START();
  *result = option->as<BoutReal>();
  CWRAP_END();
}

extern "C" int Options_get_Field3D(Field3D** result, Options* option) {
  CWRAP_START();
  **result = option->as<Field3D>(**result);
  CWRAP_END();
}

extern "C" int Options_set_int(Options* option, int value) {
  CWRAP_START();
  (*option) = value;
  CWRAP_END();
}

extern "C" int Options_set_BoutReal(Options* option, BoutReal value) {
  CWRAP_START();
  (*option) = value;
  CWRAP_END();
}

extern "C" int Options_set_Field3D(Options* option, Field3D* value) {
  CWRAP_START();
  (*option) = *value;
  CWRAP_END();
}

//////////////////////////////////////////////////////////////
// Physics Models

namespace {
  class CAPIModel : public PhysicsModel {
  public:
    void addEvolvingVar(const std::string &name) {
      evolving_names.push_back(name);
    }
    void setRHS(rhs_function function) {
      user_function = function;
    }
  protected:
    /// Initialise
    int init(bool restarting) override {
      evolving_vars.resize(evolving_names.size());
      for(std::size_t i = 0; i < evolving_names.size(); ++i) {
        solver->add(evolving_vars[i], evolving_names[i]);
      }
      return 0;
    }
    /// Calculate time derivatives
    int rhs(BoutReal t) override {
      // Put current state into an Options tree
      Options state;
      for(std::size_t i = 0; i < evolving_names.size(); ++i) {
        state[evolving_names[i]] = evolving_vars[i];
      }
      
      Options timederivs; // User will put time derivatives in this 

      // Run the user function
      if (int error = user_function(&timederivs, &state, t)) {
        return error;
      }

      // Copy time derivatives into evolving variables 
      for(std::size_t i = 0; i < evolving_names.size(); ++i) {
        ddt(evolving_vars[i]) = timederivs[evolving_names[i]].as<Field3D>(evolving_vars[i]);
      }
      return 0;
    }
  private:
    std::vector<std::string> evolving_names;
    std::vector<Field3D> evolving_vars;
    rhs_function user_function;
  };
}

extern "C" int PhysicsModel_create(PhysicsModel** model) {
  CWRAP_START();
  *model = new CAPIModel();
  CWRAP_END();
}

extern "C" int PhysicsModel_evolve(PhysicsModel* model, const char* name) {
  CWRAP_START();
  static_cast<CAPIModel*>(model)->addEvolvingVar(name);
  CWRAP_END();
}

extern "C" int PhysicsModel_set_rhs(PhysicsModel* model, rhs_function function) {
  CWRAP_START();
  static_cast<CAPIModel*>(model)->setRHS(function);
  CWRAP_END();
}

extern "C" int PhysicsModel_delete(PhysicsModel* model) {
  CWRAP_START();
  delete model;
  CWRAP_END();
}

//////////////////////////////////////////////////////////////
// Solver

extern "C" int Solver_create(Solver** solver) {
  CWRAP_START();
  // Note: This release converts a unique_ptr into a raw pointer.
  // Perhaps better to store the unique_ptr in a map and then return handle
  *solver = Solver::create().release();
  (*solver)->addMonitor(new BoutMonitor(), Solver::BACK);
  CWRAP_END();
}

extern "C" int Solver_add_model(Solver *solver, PhysicsModel *model) {
  CWRAP_START();
  solver->setModel(model);
  // This currently needs to be done after adding model
  solver->outputVars(bout::globals::dump);
  CWRAP_END();
}

extern "C" int Solver_solve(Solver *solver) {
  CWRAP_START();
  solver->solve();
  CWRAP_END();
}

extern "C" int Solver_delete(Solver *solver) {
  CWRAP_START();
  delete solver;
  CWRAP_END();
}
