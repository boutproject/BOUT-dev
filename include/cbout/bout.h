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
  Field3D* Field3D_new();
  Field3D* Field3D_zerofrom(Field3D* field);
  Field3D* Field3D_emptyfrom(Field3D* field);

  /* Destroying */
  void Field3D_delete(Field3D* field);
  
  /* Manipulating */

  /* AXPY:
     result = a * x + y 
     Can be optimised for e.g. a = 0.0, 1.0 and -1.0
  */
  void Field3D_axpy(Field3D* result, BoutReal a, Field3D* x, Field3D* y);

  /*
    Scale (multiply) a field by a given factor alpha
    x -> x * alpha
   */
  void Field3D_scale(Field3D* x, BoutReal alpha);
  
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
#ifdef __cplusplus
}
#endif
#endif /* CBOUT_H_ */
