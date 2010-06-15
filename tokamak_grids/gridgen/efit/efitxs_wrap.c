/*
C wrappers for for fortran routines READG, WRITEG etc
Fortran interface to IDL is platform specific, it
is easier to link via C wrappers
*/

#include <stdio.h>  
#include <stdlib.h>
#include "export.h" /* IDL export definitions */

#if UNDERSCORE
 #define READG readg_
 #define WRITEG writeg_
 #define READA reada_
 #define WRITEA writea_
 #define GET_NXY get_nxy_
 #define GET_PSI get_psi_
 #define PUT_PSI put_psi_
 #define GET_RZSEPS get_rzseps_
 #define PUT_RZSEPS put_rzseps_
 #define PUT_RZSEPS put_rzseps_
 #define WRITEDATA writedata_
#else
 #define READG readg
 #define WRITEG writeg
 #define READA reada
 #define WRITEA writea
 #define GET_NXY get_nxy
 #define GET_PSI get_psi
 #define PUT_PSI put_psi
 #define GET_RZSEPS get_rzseps
 #define PUT_RZSEPS put_rzseps
 #define WRITEDATA writedata
#endif



int read_data_g(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void READG();
  READG();

  return(1);
}  


int save_data_g(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void WRITEG();    
  WRITEG();

  return(1);
}  


int get_dims_g(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void GET_NXY(int nx, int ny);    

  if (argc != 2){
    fprintf(stderr,"get_nxy_c: Incorrect number of arguments\n");
    return(0);
  }

  GET_NXY((int) argv[0], (int) argv[1]);

  return(1);
}  



int get_data_g(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void GET_PSI(int nxefit, int nyefit, double** fold, 
		      double* xdim, double* zdim, double* rcentr, 
		      double* rgrid1, double* zmid, double* rmagx, double *zmagx, 
		      double* simagx, double* sibdry, double* bcentr);    

  if (argc != 13){
    fprintf(stderr,"get_psi_c: Incorrect number of arguments\n");
    return(0);
  }


  GET_PSI(
          (int) argv[0], (int) argv[1],
          (double**) argv[2],(double*) argv[3],(double*) argv[4],(double*) argv[5],
          (double*) argv[6],(double*) argv[7],(double*) argv[8],(double*) argv[9],
          (double*) argv[10],(double*) argv[11],(double*) argv[12]
          );

  return (1);
}  



int put_data_g(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void PUT_PSI(int nxefit, int nyefit, double** fold, 
		      double* xdim, double* zdim, double* rcentr, 
		      double* rgrid1, double* zmid, double* rmagx, double *zmagx, 
		      double* simagx, double* sibdry, double* bcentr);

  if (argc != 13){
    fprintf(stderr,"put_psi_c: Incorrect number of arguments\n");
    return(0);
  }


  PUT_PSI(
          (int) argv[0], (int) argv[1],
          (double**) argv[2],(double*) argv[3],(double*) argv[4],(double*) argv[5],
          (double*) argv[6],(double*) argv[7],(double*) argv[8],(double*) argv[9],
          (double*) argv[10],(double*) argv[11],(double*) argv[12]
          );

  return (1);
}  




int read_data_a(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void READA();    
  READA();

  return(1);
}  


int save_data_a(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void WRITEA();    
  WRITEA();

  return(1);
}  



int get_data_a(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void GET_RZSEPS(double *rseps1, double *zseps1, double *rseps2, double *zseps2);    

  if (argc != 4){
    fprintf(stderr,"get_rzseps_c: Incorrect number of arguments\n");
    return(0);
  }


  GET_RZSEPS((double *) argv[0], (double *) argv[1], (double *) argv[2], (double *) argv[3]);

  return (1);
}  


int put_data_a(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void PUT_RZSEPS(double *rseps1, double *zseps1, double *rseps2, double *zseps2);    

  if (argc != 4){
    fprintf(stderr,"put_rzseps_c: Incorrect number of arguments\n");
    return(0);
  }

  PUT_RZSEPS((double*) argv[0], (double*) argv[1], (double*) argv[2], (double*) argv[3]);

  return (1);
}  


int write_gridue(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void WRITEDATA(int nxm, int nym, int ixpt1, int ixpt2, int iysptrx1,
			double*** rm, double*** zm, double*** psi, 
			double*** br, double*** bz, double*** bpol, 
			double*** bphi, double*** b);    

  if (argc != 13){
    fprintf(stderr,"write_gridue: Incorrect number of arguments\n");
    return(0);
  }

  WRITEDATA(
	    (int) argv[0], (int) argv[1], (int) argv[2], (int) argv[3], (int) argv[4],
	    (double***) argv[5], (double***) argv[6], (double***) argv[7], 
	    (double***) argv[8], (double***) argv[9], (double***) argv[10], 
	    (double***) argv[11], (double***) argv[12]
	    );

  return (1);
}  
