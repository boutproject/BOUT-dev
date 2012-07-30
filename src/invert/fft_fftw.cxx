/*!************************************************************************
 * \file fft.cpp
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
 * 
 **************************************************************************/

#include <globals.hxx>
#include <options.hxx>
#include <fft.hxx>
#include <bout/constants.hxx>

#include <fftw3.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

bool fft_options = false;
bool fft_measure;

void fft_init()
{
  if(fft_options)
    return;
  //#pragma omp critical
  {
    Options *opt = Options::getRoot();
    opt = opt->getSection("fft");
    opt->get("fft_measure", fft_measure, false);
    fft_options = true;
  }
}

#ifndef _OPENMP
// Serial code
void cfft(dcomplex *cv, int length, int isign)
{
  static fftw_complex *in, *out;
  static fftw_plan pf, pb;
  static int n = 0;

  if(length != n) {
    if(n > 0) {
      fftw_destroy_plan(pf);
      fftw_destroy_plan(pb);
      fftw_free(in);
      fftw_free(out);
    } 
    fft_init();

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
    
    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;

    pf = fftw_plan_dft_1d(length, in, out, FFTW_FORWARD, flags);
    pb = fftw_plan_dft_1d(length, in, out, FFTW_BACKWARD, flags);
    
    n = length;
  }

  // Load input data
  for(int i=0;i<n;i++) {
    in[i][0] = cv[i].Real();
    in[i][1] = cv[i].Imag();
  }
  
  if(isign < 0) {
    // Forward transform
    fftw_execute(pf);
    for(int i=0;i<n;i++)
      cv[i] = dcomplex(out[i][0], out[i][1]) / ((double) n); // Normalise
  }else {
    // Backward
    fftw_execute(pb);
    for(int i=0;i<n;i++)
      cv[i] = dcomplex(out[i][0], out[i][1]);
  }
}
#else
// Parallel thread-safe version of cfft
void cfft(dcomplex *cv, int length, int isign)
{
  static fftw_complex *inall, *outall;
  static fftw_plan *pf, *pb;
  static int size = 0, nthreads;
  
  int th_id = omp_get_thread_num();
  #pragma omp critical
  {
    // Sort out memory. Also, FFTW planning routines not thread safe
    int n_th = omp_get_num_threads(); // Number of threads
    if((size != length) || (nthreads < n_th)) {
      if(size > 0) {
        // Free all memory
        for(int i=0;i<nthreads;i++) {
          fftw_destroy_plan(pf[i]);
          fftw_destroy_plan(pb[i]);
        }
        delete[] pf;
        delete[] pb;
        fftw_free(inall);
        fftw_free(outall);
      }
      
      inall = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length * n_th);
      outall = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length * n_th);
      
      pf = new fftw_plan[n_th];
      pb = new fftw_plan[n_th];
      
      unsigned int flags = FFTW_ESTIMATE;
      if(fft_measure)
        flags = FFTW_MEASURE;
      
      for(int i=0;i<n_th;i++) {
        pf[i] = fftw_plan_dft_1d(length, inall+i*length, outall+i*length, 
                                 FFTW_FORWARD, flags);
        pb[i] = fftw_plan_dft_1d(length, inall+i*length, outall+i*length,
                                 FFTW_BACKWARD, flags);
      }
      
      nthreads = n_th;
      size = length;
    }
  }
  // Get working arrays for this thread
  fftw_complex *in = inall+th_id*length;
  fftw_complex *out = outall+th_id*length;
  
  // Load input data
  for(int i=0;i<length;i++) {
    in[i][0] = cv[i].Real();
    in[i][1] = cv[i].Imag();
  }
  
  if(isign < 0) {
    // Forward transform
    fftw_execute(pf[th_id]);
    for(int i=0;i<length;i++)
      cv[i] = dcomplex(out[i][0], out[i][1]) / ((double) length); // Normalise
  }else {
    // Backward
    fftw_execute(pb[th_id]);
    for(int i=0;i<length;i++)
      cv[i] = dcomplex(out[i][0], out[i][1]);
  }
}
#endif


void ZFFT(dcomplex *cv, BoutReal zoffset, int isign, bool shift)
{
  int jz, ikz;
  BoutReal kwave;
  
  int ncz = mesh->ngz-1;
  if((isign > 0) && (mesh->ShiftXderivs) && shift) {
    // Reverse FFT
    for(jz=0;jz<ncz;jz++) {
      if (jz <= ncz/2) ikz=jz; else ikz=jz-ncz;
      kwave=ikz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
      
      // Multiply by EXP(ik*zoffset)
      cv[jz] *= dcomplex(cos(kwave*zoffset) , sin(kwave*zoffset));
    }
  }

  cfft(cv, ncz, isign);

  if((isign < 0) && (mesh->ShiftXderivs) && shift) {
    // Forward FFT
    for(jz=0;jz<ncz;jz++) {
      if (jz <= ncz/2) ikz=jz; else ikz=jz-ncz;
      kwave=ikz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
      
      // Multiply by EXP(-ik*zoffset)
      cv[jz] *= dcomplex(cos(kwave*zoffset) , -sin(kwave*zoffset));
    }
  }
}
/***********************************************************
 * Real FFTs
 ***********************************************************/

#ifndef _OPENMP
// Serial code
void rfft(BoutReal *in, int length, dcomplex *out) {
  static double *fin;
  static fftw_complex *fout;
  static fftw_plan p;
  static int n = 0;
  
  if(length != n) {
    if(n > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }
    
    fft_init();
    
    fin = (double*) fftw_malloc(sizeof(double) * length);
    fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (length/2 + 1));

    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;
    
    p = fftw_plan_dft_r2c_1d(length, fin, fout, flags);
    
    n = length;
  }
  
  for(int i=0;i<n;i++)
    fin[i] = in[i];
  
  fftw_execute(p);

  for(int i=0;i<(n/2)+1;i++)
    out[i] = dcomplex(fout[i][0], fout[i][1]) / ((double) n); // Normalise
}

void irfft(dcomplex *in, int length, BoutReal *out)
{
  static fftw_complex *fin;
  static double *fout;
  static fftw_plan p;
  static int n = 0;
  
  if(length != n) {
    if(n > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }
    
    fft_init();
    
    fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (length/2 + 1));
    fout = (double*) fftw_malloc(sizeof(double) * length);

    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;
    
    p = fftw_plan_dft_c2r_1d(length, fin, fout, flags);
    
    n = length;
  }
  
  for(int i=0;i<(n/2)+1;i++) {
    fin[i][0] = in[i].Real();
    fin[i][1] = in[i].Imag();
  }
  
  fftw_execute(p);

  for(int i=0;i<n;i++)
    out[i] = fout[i];
}

#else
// Parallel thread-safe version of rfft and irfft
void rfft(BoutReal *in, int length, dcomplex *out) {
  static double *finall;
  static fftw_complex *foutall;
  static fftw_plan *p;
  static int size = 0, nthreads;
  
  int th_id = omp_get_thread_num();
#pragma omp critical(rfft)
  {
    // Sort out memory. Also, FFTW planning routines not thread safe
    int n_th = omp_get_num_threads(); // Number of threads
    //output << "Num threads = " << n_th << " id = " << th_id << endl;
    
    if((size != length) || (nthreads < n_th)) {
      if(size > 0) {
        // Free all memory
        for(int i=0;i<nthreads;i++)
          fftw_destroy_plan(p[i]);
        delete[] p;
        fftw_free(finall);
        fftw_free(foutall);
      }
    
      fft_init();
      
      finall = (double*) fftw_malloc(sizeof(double) * length * n_th);
      foutall = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (length/2 + 1) * n_th);
      p = new fftw_plan[n_th];
      
      unsigned int flags = FFTW_ESTIMATE;
      if(fft_measure)
        flags = FFTW_MEASURE;
      
      for(int i=0;i<n_th;i++)
        p[i] = fftw_plan_dft_r2c_1d(length, finall+i*length, 
                                    foutall+i*(length/2 + 1), flags);
      size = length;
      nthreads = n_th;
    }
  }
  
  // Get working arrays for this thread
  double *fin = finall + th_id * length;
  fftw_complex *fout = foutall + th_id * (length/2 + 1);
  
  for(int i=0;i<length;i++)
    fin[i] = in[i];
  
  fftw_execute(p[th_id]);

  for(int i=0;i<(length/2)+1;i++)
    out[i] = dcomplex(fout[i][0], fout[i][1]) / ((double) length); // Normalise
}

void irfft(dcomplex *in, int length, BoutReal *out)
{
  static fftw_complex *finall;
  static double *foutall;
  static fftw_plan *p;
  static int size = 0, nthreads;
  
  int th_id = omp_get_thread_num();
  int n_th = omp_get_num_threads(); // Number of threads
#pragma omp critical(irfft)
  {
    // Sort out memory. Also, FFTW planning routines not thread safe
    
    if((size != length) || (nthreads < n_th)) {
      if(size > 0) {
        // Free all memory
        for(int i=0;i<nthreads;i++)
          fftw_destroy_plan(p[i]);
        delete[] p;
        fftw_free(finall);
        fftw_free(foutall);
      }
    
      fft_init();
      
      finall = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (length/2 + 1) * n_th);
      foutall = (double*) fftw_malloc(sizeof(double) * length * n_th);
      
      p = new fftw_plan[n_th];
      
      unsigned int flags = FFTW_ESTIMATE;
      if(fft_measure)
        flags = FFTW_MEASURE;
      
      for(int i=0;i<n_th;i++)
        p[i] = fftw_plan_dft_c2r_1d(length, finall+i*(length/2 + 1), 
                                    foutall+i*length, flags);
      size = length;
      nthreads = n_th;
    }
  }
  
  // Get working arrays for this thread
  fftw_complex *fin = finall + th_id * (length/2 + 1);
  double *fout = foutall + th_id * length;
  
  for(int i=0;i<(length/2)+1;i++) {
    fin[i][0] = in[i].Real();
    fin[i][1] = in[i].Imag();
  }
  
  fftw_execute(p[th_id]);

  for(int i=0;i<length;i++)
    out[i] = fout[i];
}
#endif

void ZFFT(BoutReal *in, BoutReal zoffset, dcomplex *cv, bool shift)
{
  int jz;
  BoutReal kwave;

  int ncz = mesh->ngz-1;

  rfft(in, ncz, cv);

  if((mesh->ShiftXderivs) && shift) {
    // Forward FFT
    for(jz=0;jz<=ncz/2;jz++) {
      kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
      
      // Multiply by EXP(-ik*zoffset)
      cv[jz] *= dcomplex(cos(kwave*zoffset) , -sin(kwave*zoffset));
    }
  }
}

void ZFFT_rev(dcomplex *cv, BoutReal zoffset, BoutReal *out, bool shift)
{
  int jz;
  BoutReal kwave;
  
  int ncz = mesh->ngz-1;
  
  if((mesh->ShiftXderivs) && shift) {
    for(jz=0;jz<=ncz/2;jz++) { // Only do positive frequencies
      kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
      
      // Multiply by EXP(ik*zoffset)
      cv[jz] *= dcomplex(cos(kwave*zoffset) , sin(kwave*zoffset));
    }
  }

  irfft(cv, ncz, out);
}
