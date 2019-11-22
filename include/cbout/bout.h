#ifndef CBOUT_H_
#define CBOUT_H_

#ifdef __cplusplus
extern "C" {
#endif

  typedef double BoutReal;
  
  /*******************/
  /* Field3D objects */
  typedef struct Field3D Field3D;

  /* Creation */
  Field3D* Field3D_new();
  Field3D* Field3D_zerofrom(Field3D* field);
  Field3D* Field3D_emptyfrom(Field3D* field);

  /* Destroying */
  void Field3D_delete(Field3D* field);

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
