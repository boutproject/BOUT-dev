/*****************************************************************
 * Interface functions to read PDB files from IDL
 * Note: not a shared library since has global variables
 *****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pdb.h"    /* PDBlib definitions */
#include "idl_export.h" /* IDL export definitions */
#include "utils.h"


/* GLOBAL VARIABLES */
PDBfile *pdb_file;

int pdb_nlist;			/* number of variables in a file */
char **pdb_var_names;	/* list of variable names */

#define MAX_PDB_STR_LEN 256
char buffer[MAX_PDB_STR_LEN];

/**** OPEN & CLOSE A FILE ****/

int pdb_open(int argc, void *argv[])
{
  IDL_STRING *filename;

  if(argc < 1) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }
  filename = (IDL_STRING*) argv[0];

  if(filename->s == (char*) NULL) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "filename is null");
    return(0);
  }

  if(argc > 1) {
    /* Open file for writing */
    if((pdb_file = PD_open(filename->s, "a")) == NULL) {
      IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
      return(0);
    }
  }else {
    /* Open file for reading */
    if ((pdb_file = PD_open(filename->s, "r")) == NULL) {
      IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
      return(0);
    }
  }
  return(1);
}

int pdb_close(int argc, void *argv[])
{
  PD_close(pdb_file);
  return(1);
}

/**** GET LIST OF VARIABLES IN FILE ****/

/* get number of variables in the file.
   argv[0] - string containing PDB filename
   argv[1] - pointer to integer for answer
*/

void uppercase(char *str, int n)
{
  int i;
  for(i=0;i!=n;i++) {
    if((str[i] >= 'a') && (str[i] <= 'z')) { /* lower case */
      str[i] += 'A' - 'a';
    }
  }
}

int idl_legalvar(int argc, void *argv[])
{
  /* Converts a given string to a legal IDL variable name */
  IDL_STRING *name;
  char* str;
  int n, i;

  if(argc != 1) {
    return(0);
  }
  name = (IDL_STRING*) argv[0];
  str = name->s;
  n = (int) name->slen;
  if(n == 0)
    return(0);

  /* convert to upper case */
  uppercase(str, n);

  /* Check first letter is a letter */  

  if((str[0] < 'A') || (str[0] > 'Z')) {
    str[0] = 'A';
  }
  for(i=1;i!=n;i++) {
    if( ((str[i] < 'A') || (str[i] > 'Z')) && ((str[i] < '0') || (str[i] > '9')) && (str[i] != '_') && (str[i] != '$')) {
      /* not a legal character */
      str[i] = '_';
    }
  }
  return(1);
}

int pdb_list_vars(int argc, void *argv[])
{
  pdb_var_names = PD_ls(pdb_file, NULL, NULL, &pdb_nlist);
  if(pdb_var_names == (char**) NULL) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(pdb_nlist);
}

/* get the name of the nth variable in the list */
char *pdb_list_varname(int argc, void *argv[])
{
  short *n;
  if(argc != 1) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return((char*) NULL);
  }
  n = (short*) argv[0];
  if((*n < 0) || (*n >= pdb_nlist)) {
    sprintf(buffer,  "n must be within range 0 to %d", pdb_nlist-1);
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, buffer);
    return((char*) NULL);
  }
  return(pdb_var_names[*n]);
}

/* free the data */
int pdb_list_free(int argc, void *argv[])
{
  return(1);
  /* currently causes segmentation fault
  for(i=0;i<pdb_nlist;i++) {
     free(pdb_var_names[i]);
  }
  free(pdb_var_names);
  return(1);
  */
}

/**** RETURN PROPERTIES OF A VARIABLE ****/

/* query properties of variable and return type as string */
char* pdb_query_type(int argc, void **argv[])
{
  IDL_STRING *name;
  syment *ep;
  static char type[32];
  char *typ;
  char *token, memb[32];
  
  if(argc != 1) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return((char*) NULL);
  }
  
  name = (IDL_STRING*) argv[0];
  
  strcpy(memb, name->s);
  token = strtok(memb, ".([");

  ep = PD_query_entry(pdb_file, token, NULL);
  if(ep == NULL) {
    sprintf(buffer, "Error querying variable %s: %s\n", name->s, PD_get_error());
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, buffer);
    return((char*) NULL);
  }
  typ = PD_entry_type(ep);
  if(typ == (char*) NULL) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "Error: Null type string");
    return((char*) NULL);
  }
  strncpy(type, typ, 31);	/* copy into static array */
  type[31] = 0;	/* make sure it's null terminated */
  
  return(type);
}

int pdb_query_ndim(int argc, void *argv[])
{
  IDL_STRING *name;
  syment *ep;
  int nd;
  char *token, memb[32];
  dimdes *dims;
  

  if(argc != 1) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  
  strcpy(memb, name->s);
  token = strtok(memb, ".([");

  ep = PD_query_entry(pdb_file, token, NULL);
  if(ep == (syment*) NULL) {
    sprintf(buffer, "Error querying variable %s: %s\n", name->s, PD_get_error());
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, buffer);
    return(0);
  }
  dims = PD_entry_dimensions(ep);
  nd = 0;
  if(dims != (dimdes*) NULL) {
    do {
      dims = dims->next;
      nd++; 
    }while(dims != (dimdes*) NULL);
  }
  return(nd);
}

/* get dimensions of array. has 2*ndim + 1 elements - (min, max) pairs & total nr*/
int pdb_query_dims(int argc, void *argv[])
{
	IDL_STRING *name;
	syment *ep;
	char *token, memb[32];
	IDL_LONG *idl_dims;
	int i;
	long total;
	dimdes* dims;

	if(argc != 2) {
	  IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
	  return(0);
	}

	name = (IDL_STRING*) argv[0];
	idl_dims = (IDL_LONG*) argv[1];
 
	strcpy(memb, name->s);
	token = strtok(memb, ".([");

	ep = PD_query_entry(pdb_file, token, NULL);
	if(ep == (syment*) NULL) {
	  sprintf(buffer, "Error querying variable %s: %s\n", name->s, PD_get_error());
	  IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, buffer);
	  return(0);
	}
	
	dims = PD_entry_dimensions(ep);
	if(dims != (dimdes*) NULL) {
	   i = 0;
	   total = 1;
	   do {
	     idl_dims[i] = (IDL_LONG) dims->index_min;
	     i++;
	     idl_dims[i] = (IDL_LONG) dims->index_max;
	     i++;
	     total *= dims->index_max - dims->index_min + 1;
	     dims = dims->next;
	   }while(dims != (dimdes*) NULL);
	   idl_dims[i] = (IDL_LONG) total;
	}else {
	  idl_dims[0] = 1;
	}

	return(1);
}

/**************** read data from PDB files **********************/

char* pdb_read_string(int argc, void *argv[])
{
  IDL_STRING *name;
  IDL_LONG *length;
  long ind1d[3];
  int len, i;
  
  if(argc != 2) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "pdb_read_string: incorrect number of arguments");
    return((char*) NULL);
  }
  name = (IDL_STRING*) argv[0];
  length = (IDL_LONG*) argv[1];
  
  len = (int) (*length);
  if(len > MAX_PDB_STR_LEN-1) {
    len = MAX_PDB_STR_LEN-1;
  }
  
  ind1d[0] = 0L;
  ind1d[1] = len-1;
  ind1d[2] = 1L;

  /*clean up the buffer*/
  for (i=0; i<MAX_PDB_STR_LEN;i++) buffer[i]=0;

  if (PD_read_as_alt(pdb_file, name->s, "char", buffer, ind1d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return((char*) NULL);
  }
  buffer[MAX_PDB_STR_LEN-1] = 0;
  return(buffer);
}

int pdb_read_0d_int(int argc, void *argv[])
{
  int var;
  IDL_STRING *name;
  IDL_LONG ivar;
  
  if(argc != 1) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }
  
  name = (IDL_STRING*) argv[0];
  var = 0;
  if (PD_read_as(pdb_file, name->s, "integer", &var) == FALSE)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
  ivar = (IDL_LONG) var;
  return(ivar);
}


int pdb_read_1d_int(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind1d[3];
  IDL_LONG *inds;
  int *array;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (int*) argv[2];

  ind1d[0] = inds[0];
  ind1d[1] = inds[1];
  ind1d[2] = inds[2];

  if (PD_read_as_alt(pdb_file, name->s, "int", array, ind1d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(1);
}


float pdb_read_0d_float(int argc, void *argv[])
{
  float var;
  IDL_STRING *name;
  
  if(argc != 1) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];

  var = 0;
  if (PD_read_as(pdb_file, name->s, "float", &var) == FALSE)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
  return(var);
}

int pdb_read_1d_float(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind1d[3];
  IDL_LONG *inds;
  float *array;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (float*) argv[2];

  ind1d[0] = inds[0];
  ind1d[1] = inds[1];
  ind1d[2] = inds[2];

  if (PD_read_as_alt(pdb_file, name->s, "float", array, ind1d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(1);
}

int pdb_read_2d_float(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind2d[6], xsize, ysize, i, j;
  IDL_LONG *inds;
  float *array;
  float **arr2d;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (float*) argv[2];

  
  for(i=0;i!=6;i++) {
    ind2d[i] = (long) inds[i];
  }
  
  /* allocate memory for the 2d array */
  xsize = ind2d[1] - ind2d[0] + 1;
  ysize = ind2d[4] - ind2d[3] + 1;

  arr2d = fmatrix(xsize, ysize);
  
  if (PD_read_as_alt(pdb_file, name->s, "float", *arr2d, ind2d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  
  /* now translate into IDL's format array[x][y] = array[x + xsize*y]*/
  
  for(i=0;i<xsize;i++) {
    for(j=0;j<ysize;j++) {
      array[i + xsize * j] = arr2d[i][j];
    }
  }
  free_fmatrix(arr2d);
  return(1);
}

int pdb_read_3d_float(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind3d[9], xsize, ysize, zsize, i, j, k;
  IDL_LONG *inds;
  float *array;
  float ***arr3d;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (float*) argv[2];

  
  for(i=0;i!=9;i++) {
    ind3d[i] = (long) inds[i];
  }
  
  /* allocate memory for the 3d array */
  xsize = ind3d[1] - ind3d[0] + 1;
  ysize = ind3d[4] - ind3d[3] + 1;
  zsize = ind3d[7] - ind3d[6] + 1;
  arr3d = f3tensor(xsize, ysize, zsize);
  
  if (PD_read_as_alt(pdb_file, name->s, "float", **arr3d, ind3d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  
  /* now translate into IDL's format array[x][y][z] = array[x + xsize*y + xsize*ysize*z]*/
  
  for(i=0;i<xsize;i++) {
    for(j=0;j<ysize;j++) {
      for(k=0;k<zsize;k++) {
	array[i + xsize*j + xsize*ysize*k] = arr3d[i][j][k];
      }
    }
  }
  free_f3tensor(arr3d);
  return(1);
}

int pdb_read_4d_float(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind4d[12], xsize, ysize, zsize, tsize, i, j, k, t;
  IDL_LONG *inds;
  float *array;
  float ****arr4d;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (float*) argv[2];

  for(i=0;i!=12;i++) {
    ind4d[i] = (long) inds[i];
  }
  
  /* allocate memory for the 4d array */
  xsize = ind4d[1] - ind4d[0] + 1;
  ysize = ind4d[4] - ind4d[3] + 1;
  zsize = ind4d[7] - ind4d[6] + 1;
  tsize = ind4d[10] - ind4d[9] + 1;
  arr4d = f4tensor(xsize, ysize, zsize, tsize);
  
  if (PD_read_as_alt(pdb_file, name->s, "float", ***arr4d, ind4d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }

  /* now translate into IDL's format array[x][y][z][t] = array[x + xsize*y + xsize*ysize*z + xsize*ysize*zsize*t]*/
  for(i=0;i<xsize;i++) {
    for(j=0;j<ysize;j++) {
      for(k=0;k<zsize;k++) {
	for(t=0;t<tsize;t++) {
	  array[i + xsize*j + xsize*ysize*k + xsize*ysize*zsize*t] = arr4d[i][j][k][t];
	}
      }
    }
  }
  free_f4tensor(arr4d);
  return(1);
}

/********************** WRITE DATA TO PDB FILES ************************/

int pdb_write_string(int argc, void *argv[])
{
  IDL_STRING *name;
  IDL_STRING *var;
  long ind1d[3];
  int len;
  
  if(argc != 2) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "pdb_write_string: incorrect number of arguments");
    return(0);
  }
  name = (IDL_STRING*) argv[0];
  var = (IDL_STRING*) argv[1];
  
  len = (int) var->slen;
  
  ind1d[0] = 0L;
  ind1d[1] = len-1;
  ind1d[2] = 1L;

  if (PD_write_alt(pdb_file, name->s, "char", var->s, 1, ind1d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(1);
}

int pdb_write_0d_int(int argc, void *argv[])
{
  int var;
  IDL_STRING *name;
  IDL_LONG *value;
  
  if(argc != 2) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }
  
  name = (IDL_STRING*) argv[0];
  value = (IDL_LONG*) argv[1];
  var = (int) (*value);

  if (PD_write(pdb_file, name->s, "integer", &var) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(1);
}

int pdb_write_1d_int(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind1d[3];
  IDL_LONG *inds;
  int *array;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (int*) argv[2];

  ind1d[0] = (long) inds[0];
  ind1d[1] = (long) inds[1];
  ind1d[2] = (long) inds[2];

  if (PD_write_alt(pdb_file, name->s, "integer", array, 1, ind1d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(1);
}

int pdb_write_0d_float(int argc, void *argv[])
{
  float *var;
  IDL_STRING *name;
  
  if(argc != 2) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  var = (float*) argv[1];

  if (PD_write(pdb_file, name->s, "float", var) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(1);
}

int pdb_write_1d_float(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind1d[3];
  IDL_LONG *inds;
  float *array;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (float*) argv[2];

  ind1d[0] = inds[0];
  ind1d[1] = inds[1];
  ind1d[2] = inds[2];

  if (PD_write_alt(pdb_file, name->s, "float", array, 1, ind1d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(1);
}

int pdb_write_2d_float(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind2d[6], xsize, ysize, i, j;
  IDL_LONG *inds;
  float *array;
  float **arr2d;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (float*) argv[2];

  
  for(i=0;i!=6;i++) {
    ind2d[i] = (long) inds[i];
  }
  
  /* allocate memory for the 2d array */
  xsize = ind2d[1] - ind2d[0] + 1;
  ysize = ind2d[4] - ind2d[3] + 1;

  arr2d = fmatrix(xsize, ysize);
  
  /* translate from IDL's format array[x][y] = array[x+xsize*y] */
  for(i=0;i<xsize;i++) {
    for(j=0;j<ysize;j++) {
      arr2d[i][j] = array[i + xsize*j];
    }
  }

  if (PD_write_alt(pdb_file, name->s, "float", *arr2d, 2, ind2d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  
  free_fmatrix(arr2d);
  return(1);
}

int pdb_write_3d_float(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind3d[9], xsize, ysize, zsize, i, j, k;
  IDL_LONG *inds;
  float *array;
  float ***arr3d;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (float*) argv[2];

  
  for(i=0;i!=9;i++) {
    ind3d[i] = (long) inds[i];
  }
  
  /* allocate memory for the 3d array */
  xsize = ind3d[1] - ind3d[0] + 1;
  ysize = ind3d[4] - ind3d[3] + 1;
  zsize = ind3d[7] - ind3d[6] + 1;
  arr3d = f3tensor(xsize, ysize, zsize);
  
  /* now translate from IDL's format array[x][y][z] = array[x + xsize*y + xsize*ysize*z]*/
  
  for(i=0;i<xsize;i++) {
    for(j=0;j<ysize;j++) {
      for(k=0;k<zsize;k++) {
	arr3d[i][j][k] = array[i + xsize*j + xsize*ysize*k];
      }
    }
  }

  if (PD_write_alt(pdb_file, name->s, "float", **arr3d, 3, ind3d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  
  
  free_f3tensor(arr3d);
  return(1);
}

int pdb_write_4d_float(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind4d[12], xsize, ysize, zsize, tsize, i, j, k, t;
  IDL_LONG *inds;
  float *array;
  float ****arr4d;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (float*) argv[2];

  for(i=0;i!=12;i++) {
    ind4d[i] = (long) inds[i];
  }
  
  /* allocate memory for the 4d array */
  xsize = ind4d[1] - ind4d[0] + 1;
  ysize = ind4d[4] - ind4d[3] + 1;
  zsize = ind4d[7] - ind4d[6] + 1;
  tsize = ind4d[10] - ind4d[9] + 1;
  arr4d = f4tensor(xsize, ysize, zsize, tsize);
  
  /* now translate from IDL's format array[x][y][z][t] = array[x + xsize*y + xsize*ysize*z + xsize*ysize*zsize*t]*/
  for(i=0;i<xsize;i++) {
    for(j=0;j<ysize;j++) {
      for(k=0;k<zsize;k++) {
	for(t=0;t<tsize;t++) {
	  arr4d[i][j][k][t] = array[i + xsize*j + xsize*ysize*k + xsize*ysize*zsize*t];
	}
      }
    }
  }

  if (PD_write_alt(pdb_file, name->s, "float", ***arr4d, 4, ind4d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }

  free_f4tensor(arr4d);
  return(1);
}

int pdb_write_0d_double(int argc, void *argv[])
{
  double *var;
  IDL_STRING *name;
  
  if(argc != 2) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  var = (double*) argv[1];

  if (PD_write(pdb_file, name->s, "double", var) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(1);
}

int pdb_write_1d_double(int argc, void *argv[])
{
  IDL_STRING *name;
  long ind1d[3];
  IDL_LONG *inds;
  double *array;
  
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }

  name = (IDL_STRING*) argv[0];
  inds = (IDL_LONG*) argv[1];
  array = (double*) argv[2];

  ind1d[0] = inds[0];
  ind1d[1] = inds[1];
  ind1d[2] = inds[2];

  if (PD_write_alt(pdb_file, name->s, "double", array, 1, ind1d) == FALSE) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, PD_get_error());
    return(0);
  }
  return(1);
}
