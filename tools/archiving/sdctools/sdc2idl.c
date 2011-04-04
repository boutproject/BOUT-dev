/**********************************************
 * A library to read and write SDC files 
 * from IDL
 **********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "export.h" /* IDL export definitions */
#include "sdclib.h"

#define FILE_MAGIC "BCD 1.0"

FILE *fp;
SDCfile *f;
char *name;
int nd;
long *ind;
int size;

/*** READING ***/

int sdc2idl_open(int argc, void *argv[])
{
  char buffer[256], c;
  int n;

  IDL_STRING *filename;

  if(argc != 1) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(1);
  }

  filename = (IDL_STRING*) argv[0];

  if(filename->s == (char*) NULL) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "filename is null");
    return(1);
  }
  
  /* Open file */
  
  if((fp = fopen(filename->s, "rb")) == (FILE*) NULL) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "Could not open file");
    return(1);
  }
  
  /* Read the magic string */
  
  n = strlen(FILE_MAGIC);
  READ_VAR(buffer, 1, n, fp, 1);
  if(strncmp(buffer, FILE_MAGIC, n) != 0) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "Not a valid file format");
    fclose(fp);
    return(1);
  }
  
  /* Read sizes of variables */
  READ_VAR(&c, 1, 1, fp, 1);
  if(c != sizeof(int)) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "Size of int doesn't match");
    fclose(fp);
    return(1);
  }
  READ_VAR(&c, 1, 1, fp, 1);
  if(c != sizeof(float)) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "Size of float doesn't match");
    fclose(fp);
    return(1);
  }
  READ_VAR(&c, 1, 1, fp, 1);
  if(c != sizeof(long)) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "Size of long doesn't match");
    fclose(fp);
    return(1);
  }

  printf("Sizes: %d, %d, %d\n", (int) sizeof(int), (int) sizeof(float), (int) sizeof(long));

  READ_STRING(name, fp, 1);

  READ_VAR(&nd, sizeof(int), 1, fp, 1);
  ind = (long*) malloc(sizeof(long)*3*(nd-1));
  READ_VAR(ind, sizeof(long), 3*(nd-1), fp, 1);
  READ_VAR(&size, sizeof(int), 1, fp, 1);

  if((f = sdc_open(fp)) == SDC_NULL) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "Invalid SDC data header");
    fclose(fp);
    return(1);
  }

  printf("Size: %d\n", size);

  return(0);
}

char* sdc2idl_get_name(int argc, void *argv[])
{
  return(name);
}

int sdc2idl_get_ndim(int argc, void *argv[])
{
  return(nd);
}

int sdc2idl_get_dims(int argc, void *argv[])
{
  IDL_LONG *dims;
  int i;

  if(argc != 1) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(1);
  }

  dims = (IDL_LONG*) argv[0];
  
  for(i=0;i<(nd-1);i++) {
    dims[2*i] = (IDL_LONG) ind[3*i];
    dims[2*i+1] = (IDL_LONG) ind[3*i+1];
  }
  dims[2*(nd-1)] = (IDL_LONG) 0;
  dims[2*(nd-1)+1] = (IDL_LONG) f->nt-1;
  return(0);
}

int sdc2idl_read(int argc, void *argv[])
{
  float *array, *a;
  int t, mint, maxt;

  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, "incorrect number of arguments");
    return(0);
  }
  
  mint = (int)  *( (IDL_LONG*) argv[0] );
  maxt = (int)  *( (IDL_LONG*) argv[1] );
  array = (float*) argv[2];
  
  a = array;

  for(t=mint;t<=maxt;t++) {
    sdc_read(f, t, a);
    a += size;
  }
  return(0);
}

int sdc2idl_close()
{
  sdc_close(f);
  fclose(fp);
  
  return(0);
}

/*** WRITING ***/

int sdc2idl_newfile(int argc, void *argv[])
{
  
  return(0);
}


