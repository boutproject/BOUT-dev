#ifndef CBOUT_H_
#define CBOUT_H_

#ifdef __cplusplus
extern "C" {
#endif

  typedef double BoutReal;

  #define BOUT_SUCCESS 0
  
  /*******************/
  /* Initialisation  */

  int bout_initialise(int argc, char** argv);
  int bout_finalise();
  
  /*******************/
  /* Field3D objects */
  typedef struct Field3D Field3D;
  
  /* Creation */
  int Field3D_create(Field3D** field);
  Field3D* Field3D_new_zerofrom(Field3D* field);
  Field3D* Field3D_new_emptyfrom(Field3D* field);

  /* Destroying */
  void Field3D_delete(Field3D* field);
  
  /* Manipulating */

  /* AXPY:
     result = a * x + y 
     Can be optimised for e.g. a = 0.0, 1.0 and -1.0
  */
  int Field3D_axpy(Field3D* result, BoutReal a, Field3D* x, Field3D* y);

  /*
    Scale (multiply) a field by a given factor alpha
    x -> x * alpha
   */
  int Field3D_scale(Field3D* x, BoutReal alpha);
  
  /* Element access */

  int Field3D_getarray(Field3D* x, BoutReal **data);
  
  
  /* Get values using index */
  /*
    BoutReal Field3D_index3d(Field3D *field, int x, int y, int z);
    BoutReal Field3D_index1d(Field3D *field, int ind);
  */

  /* Set values using index */
  /*
  BoutReal Field3D_set_index3d(Field3D *field, int x, int y, int z, BoutReal value);
  BoutReal Field3D_set_index1d(Field3D *field, int ind, BoutReal value);
  */

  /*******************/
  /* Options         */

  typedef struct Options Options;

  Options* Options_root();
  Options* Options_get(Options* option, const char* name);
  
  int Options_get_int(int* result, Options* option);
  int Options_get_BoutReal(BoutReal* result, Options* option);
  int Options_get_Field3D(Field3D** result, Options* option);
  
  int Options_set_int(Options* option, int value);
  int Options_set_BoutReal(Options* option, BoutReal value);
  int Options_set_Field3D(Options* option, Field3D* value);
  
  /*******************/
  /* Physics Models  */
  
  typedef struct PhysicsModel PhysicsModel;

  /* */
  int PhysicsModel_create(PhysicsModel** model);

  /* Add a variable to be evolved */
  int PhysicsModel_evolve(PhysicsModel* model, const char* name);

  /* Function type for user-defined RHS functions */
  typedef int (*rhs_function)(Options* ddt, Options* state, BoutReal time);
  
  int PhysicsModel_set_rhs(PhysicsModel* model, rhs_function function);

  int PhysicsModel_delete(PhysicsModel* model);
  
  /*******************/
  /* Solvers         */

  typedef struct Solver Solver;

  int Solver_create(Solver** solver);
  int Solver_add_model(Solver *solver, PhysicsModel *model);
  int Solver_solve(Solver *solver);
  int Solver_delete(Solver *solver);
  
#ifdef __cplusplus
}
#endif
#endif /* CBOUT_H_ */
