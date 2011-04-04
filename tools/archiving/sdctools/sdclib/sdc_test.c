
#include <stdio.h>
#include <math.h>

#include "sdclib.h"

#define PI 3.14159265

int main()
{
  FILE *fp;
  SDCfile *f;
  int i;
  float val[2];
  int n;

  n = 100;
  

  /* Open file */
  fp = fopen("test.sdc", "wb");

  /* Create new SDC data block */
  f = sdc_newfile(fp, 2, 10, 50, 1);

  for(i=0;i<n;i++) {
    val[0] = sin(((float) i)*2.0*PI/25.0);
    val[1] = cos(((float) i)*2.0*PI/60.0);
    
    sdc_write(f, val);
  }
  
  printf("Raw data takes %d bytes\n", (int) (n*2*sizeof(float)));

  sdc_close(f);

  fclose(fp);
  
  /* Re-open file */
  
  fp = fopen("test.sdc", "rb");
  
  f = sdc_open(fp);
  
  for(i=0;i<n;i++) {
    sdc_read(f, i, val);
    printf("%d: (%f,%f) (%f,%f)\n", i, 
	   val[0], sin(((float) i)*2.0*PI/25.0),
	   val[1], cos(((float) i)*2.0*PI/60.0));
  }

  sdc_close(f);
  fclose(fp);
  
  return(0);
}

