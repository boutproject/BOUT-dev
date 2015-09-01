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

#include <utils.hxx>

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
    in[i][0] = cv[i].real();
    in[i][1] = cv[i].imag();
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
    in[i][0] = cv[i].real();
    in[i][1] = cv[i].imag();
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
    
    // fftw call
    // Plan a real-input/complex-output discrete Fourier transform (DFT)
    // in 1 dimensions. Returns a fftw_plan (containing pointers etc.)
    p = fftw_plan_dft_r2c_1d(length, fin, fout, flags);
    
    n = length;
  }
  
  for(int i=0;i<n;i++)
    fin[i] = in[i];
  
  // fftw call executing the fft
  fftw_execute(p);

  for(int i=0;i<(n/2)+1;i++)
    out[i] = dcomplex(fout[i][0], fout[i][1]) / ((double) n); // Normalise
}

//Hacked 2D version
void rfft(BoutReal **in, const int length1, const int length2, dcomplex **out, bool transpose) {
  static double *fin;
  static fftw_complex *fout;
  static fftw_plan p;
  static int n2 = 0;
  static int n1 = 0;
  static int nk = 0;
  static bool tpose=false;
  if(length1 != n1 || length2 != n2 || transpose != tpose) {
    if(n2 > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }
    
    fft_init();

    n1 = length1;
    n2 = length2;
    nk = (n2/2)+1;
    tpose = transpose;
    
    fin = (double*) fftw_malloc(sizeof(double) * n1 * n2);
    fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1 * nk);

    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;

    /*
      Arguments are
        rank : Dimensionality of transform (1d-->1)
	size : number of data elements (length2)
	howmany : number to do (length1)
	input : fin
	inembed : Can choose to FFT a sub-array, pass NULL to ignore
	istride : Stride to use in input array (1 -- or length1)
	idist : kth transorm at input+k*idist (this is length2 -- or 1)
	output : fout
	onembed : Can choose to FFT a sub-array, pass NULL to ignore
	ostride : Stride to use in output array (1 -- or length1)
	odist : kth transorm at output+k*odist (this is nk -- or 1)
	flags : controls transform a bit
    */
    int const * l2 = &length2;
    int istride=length1, idist=1;
    int ostride=length1, odist=1;
    if(tpose){
      istride=1 ; idist=n2;
      ostride=1 ; odist=nk;
    }
    p = fftw_plan_many_dft_r2c(1,l2,length1,fin,
			       NULL,istride,idist,fout,
			       NULL,ostride,odist,flags);
  }

  //Copy data
  if(tpose){
    for(int i=0;i<n1;i++)
      for(int j=0;j<n2;j++)
	fin[j+i*n2] = in[i][j];
  }else{
    for(int i=0;i<n1;i++)
      for(int j=0;j<n2;j++)
	fin[i+j*n1] = in[i][j];
  }
  //Do the transform
  fftw_execute(p);

  //Copy data back
  BoutReal fac=1.0/(double)n2;
  if(tpose){
    for(int i=0;i<n1;i++)
      for(int j=0;j<nk;j++)
	out[i][j] = dcomplex(fout[j+i*nk][0], fout[j+i*nk][1])*fac;// / ((double) n2); // Normalise
  }else{
    for(int i=0;i<n1;i++)
      for(int j=0;j<nk;j++)
	out[i][j] = dcomplex(fout[i+j*n1][0], fout[i+j*n1][1])*fac;// / ((double) n2); // Normalise
  };
}

//Hacked 3D version
void rfft(BoutReal ***in, const int length1, const int length2, const int length3, dcomplex ***out, bool transpose) {
  static double *fin;
  static fftw_complex *fout;
  static fftw_plan p;
  static bool first=true;
  static int n3 = 0;
  static int n2 = 0;
  static int n1 = 0;
  static int nk = 0;
  static int nmany = 0;
  static bool tpose=false;
  if(first || transpose != tpose) {

    if(n3 > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }
    
    fft_init();

    //Set static vars
    n1 = length1;
    n2 = length2;
    n3 = length3;
    nmany = n1*n2;
    nk = (n3/2)+1;
    tpose = transpose;
    first=false;

    fin = (double*) fftw_malloc(sizeof(double) * nmany * n3);
    fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nmany * nk);

    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;

    int const * sz = &n3;
    int istride=nmany, idist=1;
    int ostride=nmany, odist=1;
    if(tpose){
      istride=1 ; idist=n3;
      ostride=1 ; odist=nk;
    }
    p = fftw_plan_many_dft_r2c(1,sz,nmany,fin,
			       NULL,istride,idist,fout,
			       NULL,ostride,odist,flags);
  }

  //Copy data into fftw input array, really just a flattening of
  //fld.block->data, except we skip the repeat z point.
  int itot=0;
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	for(int k=0;k<n3;k++){
	  fin[k+itot*n3] = in[i][j][k];
	}
	itot++;
      }
    }
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	for(int k=0;k<n3;k++){
	  fin[itot+k*nmany] = in[i][j][k];
	}
	itot++;
      }
    }
  }

  //Do the transform
  fftw_execute(p);

  //Copy data out of fftw output array into output
  BoutReal fac=1.0/(double)n3;
  itot=0;
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	for(int k=0;k<nk;k++){
	  out[i][j][k] = dcomplex(fout[k+itot*nk][0], fout[k+itot*nk][1])*fac; // Normalise
	}
	itot++;
      }
    }
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	for(int k=0;k<nk;k++){
	  out[i][j][k] = dcomplex(fout[itot+k*nmany][0], fout[itot+k*nmany][1])*fac; // Normalise
	}
	itot++;
      }
    }
  };
}

//Hacked Field3D version
void rfft(Field3D &fld, dcomplex ***out, bool transpose) {
  static double *fin;
  static fftw_complex *fout;
  static fftw_plan p;
  static bool first=true;
  static int n3 = 0;
  static int n2 = 0;
  static int n1 = 0;
  static int nk = 0;
  static int nmany = 0;
  static bool tpose=false;
  if(first || transpose != tpose) {

    if(n3 > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }
    
    fft_init();

    //Set static vars
    n1 = mesh->ngx;
    n2 = mesh->ngy;
    n3 = mesh->ngz-1;
    nmany = n1*n2;
    nk = (n3/2)+1;
    tpose = transpose;
    first=false;

    fin = (double*) fftw_malloc(sizeof(double) * nmany * n3);
    fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nmany * nk);

    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;

    int const * sz = &n3;
    int istride=nmany, idist=1;
    int ostride=nmany, odist=1;
    if(tpose){
      istride=1 ; idist=n3;
      ostride=1 ; odist=nk;
    }
    p = fftw_plan_many_dft_r2c(1,sz,nmany,fin,
			       NULL,istride,idist,fout,
			       NULL,ostride,odist,flags);
  }

  BoutReal ***r=fld.getData();

  //Copy data into fftw input array, really just a flattening of
  //fld.block->data, except we skip the repeat z point.
  int itot=0;
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	for(int k=0;k<n3;k++){
	  fin[k+itot*n3] = r[i][j][k];
	}
	itot++;
      }
    }
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	for(int k=0;k<n3;k++){
	  fin[itot+k*nmany] = r[i][j][k];
	}
	itot++;
      }
    }
  }

  //Do the transform
  fftw_execute(p);

  //Copy data out of fftw output array into output
  BoutReal fac=1.0/(double)n3;
  itot=0;
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	for(int k=0;k<nk;k++){
	  out[i][j][k] = dcomplex(fout[k+itot*nk][0], fout[k+itot*nk][1])*fac; // Normalise
	}
	itot++;
      }
    }
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	for(int k=0;k<nk;k++){
	  out[i][j][k] = dcomplex(fout[itot+k*nmany][0], fout[itot+k*nmany][1])*fac; // Normalise
	}
	itot++;
      }
    }
  };
}

void rfft(Field3D &fld, bool transpose) {
  rfft(fld,fld.fft_coef,transpose);
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
    fin[i][0] = in[i].real();
    fin[i][1] = in[i].imag();
  }

  // fftw call executing the fft
  fftw_execute(p);

  for(int i=0;i<n;i++)
    out[i] = fout[i];
}


//Hacked 2D version
void irfft(dcomplex **in, const int length1, const int length2, BoutReal **out, bool transpose) {
  static double *fout;
  static fftw_complex *fin;
  static fftw_plan p;
  static int n2 = 0;
  static int n1 = 0;
  static int nk = 0;
  static bool tpose=false;
  
  if(length1 != n1 || length2 != n2 || transpose != tpose) {
    if(n2 > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }
    
    fft_init();

    n1 = length1;
    n2 = length2;
    nk = (n2/2)+1;
    tpose=transpose;
    
    fout = (double*) fftw_malloc(sizeof(double) * n1 * n2);
    fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1 * nk);

    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;

    /*
      Arguments are
        rank : Dimensionality of transform (1d-->1)
	size : number of data elements (nk)
	howmany : number to do (length1)
	input : fin
	inembed : Can choose to FFT a sub-array, pass NULL to ignore
	istride : Stride to use in input array (1 -- or length1)
	idist : kth transorm at input+k*idist (this is nk -- or 1)
	output : fout
	onembed : Can choose to FFT a sub-array, pass NULL to ignore
	ostride : Stride to use in output array (1 -- or length1)
	odist : kth transorm at output+k*odist (this is length2 -- or 1)
	flags : controls transform a bit
    */
    int const * l2 = &n2;
    int istride=length1, idist=1;
    int ostride=length1, odist=1;
    if(tpose){
      istride=1; idist=nk;
      ostride=1; odist=n2;
    }
    p = fftw_plan_many_dft_c2r(1,l2,length1,fin,
			       NULL,istride,idist,fout,
			       NULL,ostride,odist,flags);
  }

  //Copy data
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<nk;j++){
	fin[j+i*nk][0] = in[i][j].real();
	fin[j+i*nk][1] = in[i][j].imag();
      };
    };
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<nk;j++){
	fin[i+j*length1][0] = in[i][j].real();
	fin[i+j*length1][1] = in[i][j].imag();
      };
    };
  };
  
  //Do the transform
  fftw_execute(p);

  //Copy data back
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	out[i][j] = fout[j+i*n2];
      };
    };
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
	out[i][j] = fout[i+j*n1];
      };
    };
  };
}

//Hacked 3D version
void irfft(dcomplex ***in, const int length1, const int length2, const int length3, BoutReal ***out, bool transpose) {
  static double *fout;
  static fftw_complex *fin;
  static fftw_plan p;
  static bool first=true;
  static int n3 = 0;
  static int n2 = 0;
  static int n1 = 0;
  static int nk = 0;
  static int nmany = 0;
  static bool tpose=false;
  if(first || transpose != tpose) {

    if(n3 > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }
    
    fft_init();

    //Set static vars
    n1 = length1;
    n2 = length2;
    n3 = length3;
    nmany = n1*n2;
    nk = (n3/2)+1;
    tpose = transpose;
    first=false;

    fout = (double*) fftw_malloc(sizeof(double) * nmany * n3);
    fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nmany * nk);

    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;

    int const * sz = &n3;
    int istride=nmany, idist=1;
    int ostride=nmany, odist=1;
    if(tpose){
      istride=1 ; idist=nk;
      ostride=1 ; odist=n3;
    }
    p = fftw_plan_many_dft_c2r(1,sz,nmany,fin,
			       NULL,istride,idist,fout,
			       NULL,ostride,odist,flags);
  }

  //Copy data into fftw input array, really just a flattening of
  //fld.block->data, except we skip the repeat z point.
  int itot=0;
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
  	for(int k=0;k<nk;k++){
  	  fin[k+itot*nk][0] = in[i][j][k].real();
  	  fin[k+itot*nk][1] = in[i][j][k].imag();
  	}
  	itot++;
      }
    }
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
  	for(int k=0;k<nk;k++){
  	  fin[itot+k*nmany][0] = in[i][j][k].real();
  	  fin[itot+k*nmany][1] = in[i][j][k].imag();
  	}
  	itot++;
      }
    }
  }

  //Do the transform
  fftw_execute(p);

  //Copy data out of fftw output array into output
  itot=0;
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
  	for(int k=0;k<n3;k++){
  	  out[i][j][k] = fout[k+itot*n3];
  	}
  	itot++;
      }
    }
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
  	for(int k=0;k<n3;k++){
  	  out[i][j][k] = fout[itot+k*nmany];
  	}
  	itot++;
      }
    }
  };
}

//Hacked Field3D version
void irfft(Field3D &fld, dcomplex ***in, bool transpose) {
  static double *fout;
  static fftw_complex *fin;
  static fftw_plan p;
  static bool first=true;
  static int n3 = 0;
  static int n2 = 0;
  static int n1 = 0;
  static int nk = 0;
  static int nmany = 0;
  static bool tpose=false;
  if(first || transpose != tpose) {

    if(n3 > 0) {
      fftw_destroy_plan(p);
      fftw_free(fin);
      fftw_free(fout);
    }
    
    fft_init();

    //Set static vars
    n1 = mesh->ngx;
    n2 = mesh->ngy;
    n3 = mesh->ngz-1;
    nmany = n1*n2;
    nk = (n3/2)+1;
    tpose = transpose;
    first=false;

    fout = (double*) fftw_malloc(sizeof(double) * nmany * n3);
    fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nmany * nk);

    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;

    int const * sz = &n3;
    int istride=nmany, idist=1;
    int ostride=nmany, odist=1;
    if(tpose){
      istride=1 ; idist=nk;
      ostride=1 ; odist=n3;
    }
    p = fftw_plan_many_dft_c2r(1,sz,nmany,fin,
			       NULL,istride,idist,fout,
			       NULL,ostride,odist,flags);
  }

  BoutReal ***r=fld.getData();

  //Copy data into fftw input array, really just a flattening of
  //fld.block->data, except we skip the repeat z point.
  int itot=0;
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
  	for(int k=0;k<nk;k++){
  	  fin[k+itot*nk][0] = in[i][j][k].real();
  	  fin[k+itot*nk][1] = in[i][j][k].imag();
  	}
  	itot++;
      }
    }
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
  	for(int k=0;k<nk;k++){
  	  fin[itot+k*nmany][0] = in[i][j][k].real();
  	  fin[itot+k*nmany][1] = in[i][j][k].imag();
  	}
  	itot++;
      }
    }
  }

  //Do the transform
  fftw_execute(p);

  //Copy data out of fftw output array into output
  itot=0;
  if(tpose){
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
  	for(int k=0;k<n3;k++){
  	  r[i][j][k] = fout[k+itot*n3];
  	}
  	r[i][j][n3] = r[i][j][0];
  	itot++;
      }
    }
  }else{
    for(int i=0;i<n1;i++){
      for(int j=0;j<n2;j++){
  	for(int k=0;k<n3;k++){
  	  r[i][j][k] = fout[itot+k*nmany];
  	}
  	r[i][j][n3] = r[i][j][0];
  	itot++;
      }
    }
  };
}

void irfft(Field3D &fld, bool transpose) {
  irfft(fld,fld.fft_coef,transpose);
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
  double *fin = finall + th_id * length;
  fftw_complex *fout = foutall + th_id * (length/2 + 1);
  
  for(int i=0;i<length;i++)
    fin[i] = in[i];

  // fftw call executing the fft
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
    fin[i][0] = in[i].real();
    fin[i][1] = in[i].imag();
  }
  
  // fftw call executing the fft
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

  // Taking the fft of the real input in, and storing the output in cv
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

void ZFFT(BoutReal *in, int ix, int iy, dcomplex *cv, bool shift)
{
  int ncz = mesh->ngz-1;
  int nkz = 1+(ncz/2);
  static dcomplex ***phs=(dcomplex ***) NULL;
  if(phs == (dcomplex ***) NULL){
    phs=c3tensor(mesh->ngx,mesh->ngy,nkz);
    BoutReal kwave;
    BoutReal fac=2.0*PI/mesh->zlength;
    for(int jx=0;jx<mesh->ngx;jx++){
      for(int jy=0;jy<mesh->ngy;jy++){
	for(int jz=0;jz<nkz;jz++){
	  kwave=jz*fac; // wave number is 1/[rad]
	  dcomplex phase(cos(kwave*mesh->zShift[jx][jy]) , -sin(kwave*mesh->zShift[jx][jy]));
	  phs[jx][jy][jz]=phase;
	}
      }
    }
  };


  rfft(in, ncz, cv);

  if((mesh->ShiftXderivs) && shift) {
    // Forward FFT
    for(int jz=0;jz<nkz;jz++) {
      
      // Multiply by EXP(-ik*zoffset)
      cv[jz] *= phs[ix][iy][jz];
    }
  }
}

void ZFFT(BoutReal **in, int nx, int iy, dcomplex **cv, bool shift)
{
  int ncz = mesh->ngz-1;
  int nkz = 1+(ncz/2);
  static dcomplex ***phs=(dcomplex ***) NULL;
  if(phs == (dcomplex ***) NULL){
    phs=c3tensor(mesh->ngx,mesh->ngy,nkz);
    BoutReal kwave;
    BoutReal fac=2.0*PI/mesh->zlength;
    for(int jx=0;jx<mesh->ngx;jx++){
      for(int jy=0;jy<mesh->ngy;jy++){
	for(int jz=0;jz<nkz;jz++){
	  kwave=jz*fac; // wave number is 1/[rad]
	  dcomplex phase(cos(kwave*mesh->zShift[jx][jy]) , -sin(kwave*mesh->zShift[jx][jy]));
	  phs[jx][jy][jz]=phase;
	}
      }
    }
  };


  rfft(in, nx, ncz, cv, true);

  if((mesh->ShiftXderivs) && shift) {
    // Forward FFT
    for(int jx=0;jx<nx;jx++) {
      for(int jz=0;jz<nkz;jz++) {
	// Multiply by EXP(-ik*zoffset)
	cv[jx][jz] *= phs[jx][iy][jz];
      }
    }
  }
}

void ZFFT(BoutReal ***in, int nx, int ny, dcomplex ***cv, bool shift)
{
  int ncz = mesh->ngz-1;
  int nkz = 1+(ncz/2);
  static dcomplex ***phs=(dcomplex ***) NULL;
  if(phs == (dcomplex ***) NULL){
    phs=c3tensor(mesh->ngx,mesh->ngy,nkz);
    BoutReal kwave;
    BoutReal fac=2.0*PI/mesh->zlength;
    for(int jx=0;jx<mesh->ngx;jx++){
      for(int jy=0;jy<mesh->ngy;jy++){
	for(int jz=0;jz<nkz;jz++){
	  kwave=jz*fac; // wave number is 1/[rad]
	  dcomplex phase(cos(kwave*mesh->zShift[jx][jy]) , -sin(kwave*mesh->zShift[jx][jy]));
	  phs[jx][jy][jz]=phase;
	}
      }
    }
  };
  
  
  rfft(in, nx, ny, ncz, cv, true);
  
  if((mesh->ShiftXderivs) && shift) {
    // Forward FFT
    for(int jx=0;jx<nx;jx++) {
      for(int jy=0;jy<ny;jy++) {
	for(int jz=0;jz<nkz;jz++) {
	  // Multiply by EXP(-ik*zoffset)
	  cv[jx][jy][jz] *= phs[jx][jy][jz];
	}
      }
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

void ZFFT_rev(dcomplex *cv, int ix, int iy, BoutReal *out, bool shift)
{
  int ncz = mesh->ngz-1;
  int nkz = 1+(ncz/2);
  static dcomplex ***phs=(dcomplex ***) NULL;
  if(phs == (dcomplex ***) NULL){
    phs=c3tensor(mesh->ngx,mesh->ngy,nkz);
    BoutReal kwave;
    BoutReal fac=2.0*PI/mesh->zlength;
    for(int jx=0;jx<mesh->ngx;jx++){
      for(int jy=0;jy<mesh->ngy;jy++){
	for(int jz=0;jz<nkz;jz++){
	  kwave=jz*fac; // wave number is 1/[rad]
	  dcomplex phase(cos(kwave*mesh->zShift[jx][jy]) , sin(kwave*mesh->zShift[jx][jy]));
	  phs[jx][jy][jz]=phase;
	}
      }
    }
  };

  
  if((mesh->ShiftXderivs) && shift) {
    for(int jz=0;jz<=ncz/2;jz++) { // Only do positive frequencies
      // Multiply by EXP(ik*zoffset)
      cv[jz] *= phs[ix][iy][jz];
    }
  }

  irfft(cv, ncz, out);
}

void ZFFT_rev(dcomplex **cv, int nx, int iy, BoutReal **out, bool shift)
{
  int ncz = mesh->ngz-1;
  int nkz = 1+(ncz/2);
  static dcomplex ***phs=(dcomplex ***) NULL;
  if(phs == (dcomplex ***) NULL){
    phs=c3tensor(mesh->ngx,mesh->ngy,nkz);
    BoutReal kwave;
    BoutReal fac=2.0*PI/mesh->zlength;
    for(int jx=0;jx<mesh->ngx;jx++){
      for(int jy=0;jy<mesh->ngy;jy++){
	for(int jz=0;jz<nkz;jz++){
	  kwave=jz*fac; // wave number is 1/[rad]
	  dcomplex phase(cos(kwave*mesh->zShift[jx][jy]) , sin(kwave*mesh->zShift[jx][jy]));
	  phs[jx][jy][jz]=phase;
	}
      }
    }
  };

  
  if((mesh->ShiftXderivs) && shift) {
    for(int jx=0;jx<nx;jx++){
      for(int jz=0;jz<=ncz/2;jz++) { // Only do positive frequencies
	// Multiply by EXP(ik*zoffset)
	cv[jx][jz] *= phs[jx][iy][jz];
      }
    }
  }
  
  irfft(cv, nx, ncz, out, true);
}

void ZFFT_rev(dcomplex ***cv, int nx, int ny, BoutReal ***out, bool shift)
{
  int ncz = mesh->ngz-1;
  int nkz = 1+(ncz/2);
  static dcomplex ***phs=(dcomplex ***) NULL;
  if(phs == (dcomplex ***) NULL){
    phs=c3tensor(mesh->ngx,mesh->ngy,nkz);
    BoutReal kwave;
    BoutReal fac=2.0*PI/mesh->zlength;
    for(int jx=0;jx<mesh->ngx;jx++){
      for(int jy=0;jy<mesh->ngy;jy++){
	for(int jz=0;jz<nkz;jz++){
	  kwave=jz*fac; // wave number is 1/[rad]
	  dcomplex phase(cos(kwave*mesh->zShift[jx][jy]) , sin(kwave*mesh->zShift[jx][jy]));
	  phs[jx][jy][jz]=phase;
	}
      }
    }
  };

  
  if((mesh->ShiftXderivs) && shift) {
    for(int jx=0;jx<nx;jx++){
      for(int jy=0;jy<ny;jy++){
	for(int jz=0;jz<nkz;jz++) { // Only do positive frequencies
	  // Multiply by EXP(ik*zoffset)
	  cv[jx][jy][jz] *= phs[jx][jy][jz];
	}
      }
    }
  }

  irfft(cv, nx, ny, ncz, out, true);
}


//  Discrete sine transforms (B Shanahan)

void DST(BoutReal *in, int length, dcomplex *out) {
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

    //  fft_init();
   
    // Could be optimized better
    fin = (double*) fftw_malloc(sizeof(double) * 2 * length);
    fout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2 * length);


    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;

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
    out[i] = -fout[i][1]/ ((double) length-1); // Normalise 
}

void DST_rev(dcomplex *in, int length, BoutReal *out) {
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

    //fft_init();

    // Could be optimized better
    fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2 * (length-1));
    fout = (double*) fftw_malloc(sizeof(double)  * 2 * (length-1));

    unsigned int flags = FFTW_ESTIMATE;
    if(fft_measure)
      flags = FFTW_MEASURE;

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
}
