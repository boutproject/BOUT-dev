#include <bout.hxx>
#include <boutmain.hxx>
#include <fft.hxx>

int physics_init(bool restarting)
{
  // Testing FFT routines
  
  int n = 64;
  BoutReal *d = new BoutReal[n];
  BoutReal *d2 = new BoutReal[n];
  dcomplex *c = new dcomplex[n];
  
  srand(1234);
  for(int i=0;i<n;i++)
    d[i] = 2.*((double) rand()) / RAND_MAX - 1.0;
  
  BoutReal ts = MPI_Wtime();
  rfft(d, n, c);
  irfft(c, n, d2);
  BoutReal t1 = MPI_Wtime() - ts;

  for(int i=0;i<n;i++)
    output << i << ": " << d[i] << " -> " << c[i] << " -> " << d2[i] << endl;

  // Do timing
  ts = MPI_Wtime();
  for(int i=0;i<100;i++) {
    rfft(d, n, c);
    irfft(c, n, d2);
  }
  BoutReal t2 = (MPI_Wtime() - ts) / 100.;
  output << "TIMING: " << t1 << " : " << t2 << endl;

  exit(0);
}

int physics_run(BoutReal t)
{
  return 0;
}
