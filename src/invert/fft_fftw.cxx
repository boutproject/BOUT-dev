/*!************************************************************************
 * \file fft_fftw.cxx
 *
 * \brief FFT routines using external libraries
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 * Update by Yining Qin in July 26, 2021 (2021 Nersc Hackathon) to support  cuFFT
 **************************************************************************/

#include "bout/build_config.hxx"

#include <globals.hxx>
#include <options.hxx>
#include <fft.hxx>
#include <unused.hxx>

#if BOUT_HAS_FFTW
#include <bout/constants.hxx>
#include <bout/openmpwrap.hxx>

// test cuFFTW
//#include <cufftw.h>

#include <fftw3.h>

#include <cmath>

#if BOUT_USE_OPENMP 
#include <omp.h>
#endif
#else
#include <boutexception.hxx>
#endif

//---cuFFT define start----
#include "device_launch_parameters.h"
#include <cufft.h>

#define NX 3335 // Number of valid data
#define N 5335 // Length of data after padding 0
#define BATCH 1
#define BLOCK_SIZE 1024

using std::cout;
using std::endl;

//---cuFFT define end----


namespace bout {
namespace fft {

/*
 * Function: Determine if two cufftComplex arrays are equal
 * Input: idA input head of A pointer A
 * Input: idB output head of B pointer
 * Input: the number of elements in the size array
 * return: true or false
 */

bool cufftComplexIsEqual(cufftComplex *idA, cufftComplex *idB, const int size)
{
    for (int i = 0; i < size; i++)
    {
       if (abs(idA[i].x - idB[i].x) > 0.000001 || abs(idA[i].y - idB[i].y) > 0.000001)
           return false;
    }
    return true;
}

/* Function: Implements scaling of the cufftComplex array, which is multiplied by a number
 * Input: the head pointer of the idata input array
 * Output: head pointer of the odata output array
 * Input: the number of elements in the size array
 * Input: scale scale
 */
static __global__ void cufftComplexScale(cufftComplex *ida, cufftComplex *oda, const int size, float scale)
{
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
 
    if (threadID < size)
    {
        oda[threadID].x = ida[threadID].x * scale;
        oda[threadID].y = ida[threadID].y * scale;
    }
}

#define BATCH 1
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool  abort=true)
{
   if (code != cudaSuccess) 
   {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);            
        if (abort) exit(code);
   }
}


// using cuFFT for rfft, real to complex 

void rfft(MAYBE_UNUSED(const BoutReal *in), MAYBE_UNUSED(int length), MAYBE_UNUSED(dcomplex *out)) {

	cufftReal *hostInputData = (cufftReal*)malloc(length*sizeof(cufftReal));
	// copy data from input data to host
	for (int j=0; j<length; j++) hostInputData[j] = in[j];

	// --- Device side input data allocation and initialization
	cufftReal *deviceInputData; 
	gpuErrchk(cudaMalloc((void**)&deviceInputData, length * sizeof(cufftReal)));
	cudaMemcpy(deviceInputData, hostInputData, length * sizeof(cufftReal), cudaMemcpyHostToDevice);

	//--- Host side output data allocation
 	cufftComplex *hostOutputData = (cufftComplex*)malloc((length / 2 + 1) * BATCH * sizeof(cufftComplex));

	//--- Device side output data allocation
 	cufftComplex *deviceOutputData;
 	gpuErrchk(cudaMalloc((void**)&deviceOutputData, (length / 2 + 1) * sizeof(cufftComplex)));

	cufftResult cufftStatus;
	cufftHandle handle;
	// cufft 1d, real to complex, r2c
	cufftStatus = cufftPlan1d(&handle, length, CUFFT_R2C, BATCH);
	if (cufftStatus != cudaSuccess) { printf("cufftPlan1d failed!"); }       

	cufftStatus = cufftExecR2C(handle,  deviceInputData, deviceOutputData);
	if (cufftStatus != cudaSuccess) { printf("cufftExecR2C failed!"); }  
	gpuErrchk(cudaMemcpy(hostOutputData, deviceOutputData, (length / 2 + 1) * sizeof(cufftComplex), cudaMemcpyDeviceToHost));

  
  	//Normalising factor - this part is from original rfft code;
  	const BoutReal fac = 1.0 / length;
  	const int nmodes = (length/2) + 1;

  	// Store the output in out, and normalize
  	for(int j = 0; j < nmodes; j++){
    		out[j] = dcomplex( hostOutputData[j].x, hostOutputData[j].y)*fac ; // ? Normalise
	}
	
	cufftDestroy(handle);
	gpuErrchk(cudaFree(deviceOutputData));
	gpuErrchk(cudaFree(deviceInputData));

}


// using cuFFT for irfft, complex to real
//
void irfft(MAYBE_UNUSED(const dcomplex *in), MAYBE_UNUSED(int length), MAYBE_UNUSED(BoutReal *out)) {

	cufftComplex *hostInputData = (cufftComplex*)malloc(length*sizeof(cufftComplex));
	// copy data from input data to host
	for (int j=0; j<length; j++){
		 hostInputData[j].x = in[j].real();
		 hostInputData[j].y = in[j].imag();
		}

	// --- Device side input data allocation and initialization
	cufftComplex *deviceInputData; 
	gpuErrchk(cudaMalloc((void**)&deviceInputData, length * sizeof(cufftComplex)));
	cudaMemcpy(deviceInputData, hostInputData, length * sizeof(cufftComplex), cudaMemcpyHostToDevice);

	//--- Host side output data allocation
	 cufftReal *hostOutputData = (cufftReal*)malloc((length) * BATCH * sizeof(cufftReal));

	//--- Device side output data allocation
	 cufftReal *deviceOutputData;
	 gpuErrchk(cudaMalloc((void**)&deviceOutputData, length * sizeof(cufftReal)));

	cufftResult cufftStatus;
	cufftHandle handle;
	// cufft id,  complex to real, c2r
	cufftStatus = cufftPlan1d(&handle, length, CUFFT_C2R, BATCH);
	if (cufftStatus != cudaSuccess) { printf("cufftPlan1d failed!"); }       

	cufftStatus = cufftExecC2R(handle,  deviceInputData, deviceOutputData);
	if (cufftStatus != cudaSuccess) { printf("cufftExecC2R failed!"); }  
	gpuErrchk(cudaMemcpy(hostOutputData, deviceOutputData, (length ) * sizeof(cufftReal), cudaMemcpyDeviceToHost));

	 
	// Store the cufft output of to the out
	  for(int j=0;j<length;j++)
	    out[j] = hostOutputData[j];


	cufftDestroy(handle);
	gpuErrchk(cudaFree(deviceOutputData));
	gpuErrchk(cudaFree(deviceInputData));


}

/// Have we set fft_measure?
bool fft_initialised{false};
/// Should FFTW find an optimised plan by measuring various plans?
bool fft_measure{false};

void fft_init(Options* options) {
  if (fft_initialised) {
    return;
  }
  if (options == nullptr) {
    options = Options::getRoot()->getSection("fft");
  }
  fft_init((*options)["fft_measure"]
               .doc("Perform speed measurements to optimise settings?")
               .withDefault(false));
}

void fft_init(bool fft_measure) {
  bout::fft::fft_measure = fft_measure;
  fft_initialised = true;
}

/***********************************************************
 * Real FFTs
 ***********************************************************/

#if ! BOUT_USE_OPENMP
// Serial code
void b_rfft(MAYBE_UNUSED(const BoutReal *in), MAYBE_UNUSED(int length), MAYBE_UNUSED(dcomplex *out)) {
#if !BOUT_HAS_FFTW
  throw BoutException("This instance of BOUT++ has been compiled without fftw support.");

#else

  // static variables initialized once
  static double *fin;
  static fftw_complex *fout;
  static fftw_plan p;
  static int n = 0;

 
// If the current length mismatches n
  if(length != n) {
    // If n has been used before
    if(n > 0) {
      // Free the previously initialized data
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }

    fft_init();

    // Initialize the input for the fourier transformation
    fin = static_cast<double*>(fftw_malloc(sizeof(double) * length));
    // Initialize the output of the fourier transformation
    /* NOTE: Only the non-redundant output is given
     *       I.e the offset and the positive frequencies (so no mirroring
     *       around the Nyquist frequency)
     */
    fout = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * (length/2 + 1)));

    unsigned int flags = FFTW_ESTIMATE;
    if (fft_measure) {
      flags = FFTW_MEASURE;
    }

    /* fftw call
     * Plan a real-input/complex-output discrete Fourier transform (DFT)
     * in 1 dimensions. Returns a fftw_plan (containing pointers etc.)
     * r2c = real to complex
     */
    p = fftw_plan_dft_r2c_1d(length, fin, fout, flags);

    n = length;
  }

  // Put the input to fin
  for(int i=0;i<n;i++)
    fin[i] = in[i];

  // fftw call executing the fft
  fftw_execute(p);
  
  //Normalising factor
  const BoutReal fac = 1.0 / n;
  const int nmodes = (n/2) + 1;

  // Store the output in out, and normalize
  for(int i=0;i<nmodes;i++)
    out[i] = dcomplex(fout[i][0], fout[i][1]) * fac; // Normalise
 

#endif
}

void b_irfft(MAYBE_UNUSED(const dcomplex *in), MAYBE_UNUSED(int length), MAYBE_UNUSED(BoutReal *out)) {
#if !BOUT_HAS_FFTW
  throw BoutException("This instance of BOUT++ has been compiled without fftw support.");
#else
  // static variables initialized once
  static fftw_complex *fin;
  static double *fout;
  static fftw_plan p;
  static int n = 0;

  // If the current length mismatches n
  if(length != n) {
    // If n has been used before
    if(n > 0) {
      // Free the previously initialized data
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }

    fft_init();

    // Initilaize the input for the inverse fourier transformation
    /* NOTE: Only the non-redundant input is given
     *       I.e the offset and the positive frequencies (so no mirroring
     *       around the Nyquist frequency)
     */
    fin = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * (length/2 + 1)));
    // Initialize the output of the fourier transformation
    fout = static_cast<double*>(fftw_malloc(sizeof(double) * length));

    unsigned int flags = FFTW_ESTIMATE;
    if (fft_measure) {
      flags = FFTW_MEASURE;
    }

    /* fftw call
     * Plan a complex-input/real-output discrete Fourier transform (DFT)
     * in 1 dimensions. Returns a fftw_plan (containing pointers etc.)
     * c2r = complex to real
     */
    p = fftw_plan_dft_c2r_1d(length, fin, fout, flags);

    n = length;
  }

  // Store the real and imaginary parts in the proper way
  const int nmodes = (n/2) + 1;
  for(int i=0;i<nmodes;i++) {
    fin[i][0] = in[i].real();
    fin[i][1] = in[i].imag();
  }

  // fftw call executing the fft
  fftw_execute(p);

  // Store the output of the fftw to the out
  for(int i=0;i<n;i++)
    out[i] = fout[i];
#endif
}
#else
// Parallel thread-safe version of rfft and irfft
void b_rfft(MAYBE_UNUSED(const BoutReal *in), MAYBE_UNUSED(int length), MAYBE_UNUSED(dcomplex *out)) {
#if !BOUT_HAS_FFTW
  throw BoutException("This instance of BOUT++ has been compiled without fftw support.");
#else
  // ensure we are not nested
  ASSERT1(omp_get_active_level() < 2);
  static double *finall;
  static fftw_complex *foutall;
  static fftw_plan *p;
  static int size = 0, nthreads = 0;

 
 int th_id = omp_get_thread_num();
  int n_th = omp_get_num_threads(); // Number of threads

  // Sort out memory. Also, FFTW planning routines not thread safe
  if ((size != length) || (nthreads < n_th)) {
    // We make the check to see if the problem size has changed twice
    // intentionally. The first check ensures we don't pay the cost of
    // obtaining a lock for the critical section if we don't need to do
    // any work here. The second check is required to make sure that
    // only one thread does the actual setup when required. Note we can't
    // use a `single` block here as that requires all threads to reach the
    // block (implicit barrier) which may not be true in all cases (e.g.
    // if there are 8 threads but only 4 call the fft routine).
    BOUT_OMP(critical(rfft))
    if ((size != length) || (nthreads < n_th)) {
      if(size > 0) {
        // Free all memory
        for(int i=0;i<nthreads;i++)
          fftw_destroy_plan(p[i]);
        delete[] p;
        fftw_free(finall);
        fftw_free(foutall);
      }

      fft_init();

      finall = static_cast<double *>(fftw_malloc(sizeof(double) * length * n_th));
      foutall = static_cast<fftw_complex *>(
          fftw_malloc(sizeof(fftw_complex) * (length / 2 + 1) * n_th));
      p = new fftw_plan[n_th]; //Never freed

      unsigned int flags = FFTW_ESTIMATE;
      if (fft_measure) {
        flags = FFTW_MEASURE;
      }

      for(int i=0;i<n_th;i++)
        // fftw call
        // Plan a real-input/complex-output discrete Fourier transform (DFT)
        // in 1 dimensions. Returns a fftw_plan (containing pointers etc.)
        p[i] = fftw_plan_dft_r2c_1d(length, finall+i*length,
                                    foutall+i*(length/2 + 1), flags);
      size = length;
      nthreads = n_th;
    }
  }

  // Get working arrays for this thread
  double *fin = finall + th_id * size;
  fftw_complex *fout = foutall + th_id * (size / 2 + 1);

  for (int i = 0; i < size; i++)
    fin[i] = in[i];

  // fftw call executing the fft
  fftw_execute(p[th_id]);

  //Normalising factor
  const BoutReal fac = 1.0 / static_cast<BoutReal>(size);
  const int nmodes = (size / 2) + 1;

  for(int i=0;i<nmodes;i++)
    out[i] = dcomplex(fout[i][0], fout[i][1]) * fac; // Normalise
#endif
}

void b_irfft(MAYBE_UNUSED(const dcomplex *in), MAYBE_UNUSED(int length), MAYBE_UNUSED(BoutReal *out)) {
#if !BOUT_HAS_FFTW
  throw BoutException("This instance of BOUT++ has been compiled without fftw support.");
#else
  static fftw_complex *finall;
  static double *foutall;
  static fftw_plan *p;
  static int size = 0, nthreads = 0;

  int th_id = omp_get_thread_num();
  int n_th = omp_get_num_threads(); // Number of threads

  // Sort out memory. Also, FFTW planning routines not thread safe
  if ((size != length) || (nthreads < n_th)) {
    // We make the check to see if the problem size has changed twice
    // intentionally. The first check ensures we don't pay the cost of
    // obtaining a lock for the critical section if we don't need to do
    // any work here. The second check is required to make sure that
    // only one thread does the actual setup when required. Note we can't
    // use a `single` block here as that requires all threads to reach the
    // block (implicit barrier) which may not be true in all cases (e.g.
    // if there are 8 threads but only 4 call the fft routine).
    BOUT_OMP(critical(irfft))
    if ((size != length) || (nthreads < n_th)) {
      if (size > 0) {
        // Free all memory
        for (int i = 0; i < nthreads; i++)
          fftw_destroy_plan(p[i]);
        delete[] p;
        fftw_free(finall);
        fftw_free(foutall);
      }

      fft_init();

      finall = static_cast<fftw_complex *>(
          fftw_malloc(sizeof(fftw_complex) * (length / 2 + 1) * n_th));
      foutall = static_cast<double *>(fftw_malloc(sizeof(double) * length * n_th));

      p = new fftw_plan[n_th]; // Never freed

      unsigned int flags = FFTW_ESTIMATE;
      if (fft_measure) {
        flags = FFTW_MEASURE;
      }

      for (int i = 0; i < n_th; i++)
        p[i] = fftw_plan_dft_c2r_1d(length, finall + i * (length / 2 + 1),
                                    foutall + i * length, flags);
      size = length;
      nthreads = n_th;
    }
  }

  // Get working arrays for this thread
  fftw_complex *fin = finall + th_id * (size / 2 + 1);
  double *fout = foutall + th_id * size;

  const int nmodes = (size / 2) + 1;

  for(int i=0;i<nmodes;i++) {
    fin[i][0] = in[i].real();
    fin[i][1] = in[i].imag();
  }

  // fftw call executing the fft
  fftw_execute(p[th_id]);

  for (int i = 0; i < size; i++)
    out[i] = fout[i];
#endif
}
#endif
//  Discrete sine transforms (B Shanahan)

void DST(MAYBE_UNUSED(const BoutReal *in), MAYBE_UNUSED(int length), MAYBE_UNUSED(dcomplex *out)) {
#if !BOUT_HAS_FFTW
  throw BoutException("This instance of BOUT++ has been compiled without fftw support.");
#else
  static double *fin;
  static fftw_complex *fout;
  static fftw_plan p;
  static int n = 0;

  ASSERT1(length > 0);

  if(length != n) {
    if(n > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
      }

    //  fft_init();

    // Could be optimized better
    fin = static_cast<double *>(fftw_malloc(sizeof(double) * 2 * length));
    fout = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * 2 * length));

    unsigned int flags = FFTW_ESTIMATE;
    if (fft_measure) {
      flags = FFTW_MEASURE;
    }

    // fftw call
    // Plan a real-input/complex-output discrete Fourier transform (DFT)
    // in 1 dimensions. Returns a fftw_plan (containing pointers etc.)
    p = fftw_plan_dft_r2c_1d(2*(length-1), fin, fout, flags);

    n = length;
  }
  for(int i=0;i<n;i++)
    fin[i] = in[i];

  fin[0] = 0.;
  fin[length-1]=0.;

  for (int j = 1; j < length-1; j++){
      fin[j] = in[j];
      fin[2*(length-1)-j] = - in[j];
  }

  // fftw call executing the fft
  fftw_execute(p);

  out[0]=0.0;
  out[length-1]=0.0;

  for(int i=1;i<length-1;i++)
    out[i] = -fout[i][1] / (static_cast<BoutReal>(length) - 1); // Normalise
#endif
}

void DST_rev(MAYBE_UNUSED(dcomplex *in), MAYBE_UNUSED(int length), MAYBE_UNUSED(BoutReal *out)) {
#if !BOUT_HAS_FFTW
  throw BoutException("This instance of BOUT++ has been compiled without fftw support.");
#else
  static fftw_complex *fin;
  static double *fout;
  static fftw_plan p;
  static int n = 0;

 printf("call DST_rev\n");

  ASSERT1(length > 0);

  if(length != n) {
    if(n > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }

    //fft_init();

    // Could be optimized better
    fin =
        static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * 2 * (length - 1)));
    fout = static_cast<double *>(fftw_malloc(sizeof(double) * 2 * (length - 1)));

    unsigned int flags = FFTW_ESTIMATE;
    if (fft_measure) {
      flags = FFTW_MEASURE;
    }

    p = fftw_plan_dft_c2r_1d(2*(length-1), fin, fout, flags);

    n = length;
  }

  for(int i=0;i<n;i++){
    fin[i][0] = in[i].real();
    fin[i][1] = in[i].imag();
  }

  fin[0][0] = 0.; fin[0][1] = 0.;
  fin[length-1][0] = 0.; fin[length-1][1] = 0.;

  for (int j = 1; j < length-1; j++){
    fin[j][0] = 0.; fin[j][1] = -in[j].real()/2.;
    fin[2*(length-1)-j][0] = 0.; fin[2*(length-1)-j][1] =  in[j].real()/2.;
  }

  // fftw call executing the fft
  fftw_execute(p);

  out[0]=0.0;
  out[length-1]=0.0;
  for(int i=1;i<length-1;i++)
    out[i] = fout[i];
#endif
}

Array<dcomplex> rfft(const Array<BoutReal>& in) {
  ASSERT1(!in.empty());

  int size{in.size()};
  Array<dcomplex> out{(size / 2) + 1};

  bout::fft::rfft(in.begin(), size, out.begin());
  return out;
}

Array<BoutReal> irfft(const Array<dcomplex>& in, int length) {
  ASSERT1(!in.empty());
  ASSERT1(in.size() == (length / 2) + 1);

  Array<BoutReal> out{length};

  bout::fft::irfft(in.begin(), length, out.begin());
  return out;
}

} // namespace fft
} // namespace bout

