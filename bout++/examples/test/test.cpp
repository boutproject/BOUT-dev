#include "bout.h"
#include "fft.h"

int physics_init()
{
  // Testing FFT routines
  
  int n = 64;
  real *d = new real[n];
  real *d2 = new real[n];
  dcomplex *c = new dcomplex[n];
  
  srand(1234);
  for(int i=0;i<n;i++)
    d[i] = 2.*((double) rand()) / RAND_MAX - 1.0;
  
  real ts = MPI_Wtime();
  rfft(d, n, c);
  irfft(c, n, d2);
  real t1 = MPI_Wtime() - ts;

  for(int i=0;i<n;i++)
    output << i << ": " << d[i] << " -> " << c[i] << " -> " << d2[i] << endl;

  // Do timing
  ts = MPI_Wtime();
  for(int i=0;i<100;i++) {
    rfft(d, n, c);
    irfft(c, n, d2);
  }
  real t2 = (MPI_Wtime() - ts) / 100.;
  output << "TIMING: " << t1 << " : " << t2 << endl;

  exit(0);
}

int physics_run(real t)
{
  return 0;
}
