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

#include "globals.h"
#include "fft.h"

#include <fftw3.h>
#include <math.h>

bool fft_options = false;
bool fft_measure;

void fft_init()
{
  if(fft_options)
    return;

  options.setSection("fft");
  options.get("fft_measure", fft_measure, false);
  fft_options = true;
}

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

void ZFFT(dcomplex *cv, real zoffset, int isign, bool shift)
{
  int jz, ikz;
  real kwave;
  
  if((isign > 0) && (ShiftXderivs) && shift) {
    // Reverse FFT
    for(jz=0;jz<ncz;jz++) {
      if (jz <= ncz/2) ikz=jz; else ikz=jz-ncz;
      kwave=ikz*2.0*PI/zlength; // wave number is 1/[rad]
      
      // Multiply by EXP(ik*zoffset)
      cv[jz] *= dcomplex(cos(kwave*zoffset) , sin(kwave*zoffset));
    }
  }

  cfft(cv, ncz, isign);

  if((isign < 0) && (ShiftXderivs) && shift) {
    // Forward FFT
    for(jz=0;jz<ncz;jz++) {
      if (jz <= ncz/2) ikz=jz; else ikz=jz-ncz;
      kwave=ikz*2.0*PI/zlength; // wave number is 1/[rad]
      
      // Multiply by EXP(-ik*zoffset)
      cv[jz] *= dcomplex(cos(kwave*zoffset) , -sin(kwave*zoffset));
    }
  }
}
/***********************************************************
 * Real FFTs
 ***********************************************************/

void rfft(real *in, int length, dcomplex *out)
{
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

void irfft(dcomplex *in, int length, real *out)
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

void ZFFT(real *in, real zoffset, dcomplex *cv, bool shift)
{
  int jz;
  real kwave;

  rfft(in, ncz, cv);

  if((ShiftXderivs) && shift) {
    // Forward FFT
    for(jz=0;jz<=ncz/2;jz++) {
      kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
      
      // Multiply by EXP(-ik*zoffset)
      cv[jz] *= dcomplex(cos(kwave*zoffset) , -sin(kwave*zoffset));
    }
  }
}

void ZFFT_rev(dcomplex *cv, real zoffset, real *out, bool shift)
{
  int jz;
  real kwave;
  
  if((ShiftXderivs) && shift) {
    for(jz=0;jz<=ncz/2;jz++) { // Only do positive frequencies
      kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
      
      // Multiply by EXP(ik*zoffset)
      cv[jz] *= dcomplex(cos(kwave*zoffset) , sin(kwave*zoffset));
    }
  }

  irfft(cv, ncz, out);
}
